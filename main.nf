#!/usr/bin/env nextflow
tmpdir = config.tmpdir
reference = config.reference
cores = config.alignment_cores
variant_cores = config.variant_cores
genome = config.genome
SM_alignments_dir = config.SM_alignments_dir
analysis_dir = config.analysis_dir

// Define contigs here!
contig_list = ["I", "II", "III", "IV", "V", "X", "MtDNA"]
contigs = Channel.from(contig_list)

fq_concordance_script = file("fq_concordance.R")
process_concordance = file("process_concordance.R")
plot_pairwise_script = file("plot_pairwise.R")

/*
    Filtering configuration
*/
min_depth=3
qual=30
mq=40
dv_dp=0.5

println "Running Concordance on Wild Isolates"
println "Using Reference: ${genome}" 

// Construct strain and isotype lists
import groovy.json.JsonSlurper

def strain_set = []

if (test == 'true') {
    strain_json = "strain_set_test.json"
} else {
    strain_json = "strain_set.json"
}

// Strain
def strainFile = new File(strain_json)
def strainJSON = new JsonSlurper().parseText(strainFile.text)
strain_set_file = Channel.fromPath(strain_json)

strainJSON.each { SM, RG ->
    RG.each { k, v ->
        strain_set << [SM, k, v[0], v[1], v[2]]
    }
}

process setup_dirs {

    executor 'local'

    input:
        file strain_set_file

    """
        mkdir -p ${analysis_dir}
        cp ${strain_set_file} ${analysis_dir}/strain_set.json
    """
}

/*
    Fastq concordance
*/
process perform_alignment {

    cpus cores

    tag { fq_pair_id }

    input:
        set SM, RG, fq1, fq2, fq_pair_id from strain_set
    output:
        set val(fq_pair_id), file("${fq_pair_id}.bam"), file("${fq_pair_id}.bam.bai") into fq_bam_set
        set val(SM), file("${fq_pair_id}.bam"), file("${fq_pair_id}.bam.bai") into SM_aligned_bams

    
    """
        bwa mem -t ${cores} -R '${RG}' ${reference} ${fq1} ${fq2} | \\
        sambamba view --nthreads=${cores} --sam-input --format=bam --with-header /dev/stdin | \\
        sambamba sort --nthreads=${cores} --show-progress --tmpdir=${tmpdir} --out=${fq_pair_id}.bam /dev/stdin
        sambamba index --nthreads=${cores} ${fq_pair_id}.bam
    """
}

fq_bam_set.into { fq_cov_bam; fq_stats_bam; fq_idx_stats_bam }

/*
    Fastq coverage
*/
process coverage_fq {

    tag { fq_pair_id }

    input:
        set val(fq_pair_id), file("${fq_pair_id}.bam"), file("${fq_pair_id}.bam.bai") from fq_cov_bam
    output:
        file("${fq_pair_id}.coverage.tsv") into fq_coverage


    """
        bam coverage ${fq_pair_id}.bam > ${fq_pair_id}.coverage.tsv
    """
}


process coverage_fq_merge {

    publishDir analysis_dir + "/fq", mode: 'copy'

    input:
        file fq_set from fq_coverage.toSortedList()

    output:
        file("fq_coverage.full.tsv")
        file("fq_coverage.tsv")

    """
        echo -e 'fq\\tcontig\\tstart\\tend\\tproperty\\tvalue' > fq_coverage.full.tsv
        cat ${fq_set} >> fq_coverage.full.tsv

        cat <(echo -e 'fq\\tcoverage') <( cat fq_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6) > fq_coverage.tsv
    """
}

/*
    fq idx stats
*/

process fq_idx_stats {
    
    tag { fq_pair_id }

    input:
        set val(fq_pair_id), file("${fq_pair_id}.bam"), file("${fq_pair_id}.bam.bai") from fq_idx_stats_bam
    output:
        file fq_idxstats into fq_idxstats_set

    """
        samtools idxstats ${fq_pair_id}.bam | awk '{ print "${fq_pair_id}\\t" \$0 }' > fq_idxstats
    """
}

process fq_combine_idx_stats {

    publishDir analysis_dir + "/fq", mode: 'copy'

    input:
        file("?.stat.txt") from fq_idxstats_set.toSortedList()

    output:
        file("fq_bam_idxstats.tsv")

    """
        echo -e "SM\\treference\\treference_length\\tmapped_reads\\tunmapped_reads" > fq_bam_idxstats.tsv
        cat *.stat.txt >> fq_bam_idxstats.tsv
    """

}

/*
    fq bam stats
*/

process fq_bam_stats {

    tag { fq_pair_id }

    input:
        set val(fq_pair_id), file("${fq_pair_id}.bam"), file("${fq_pair_id}.bam.bai") from fq_stats_bam

    output:
        file 'bam_stat' into fq_bam_stat_files

    """
        cat <(samtools stats ${fq_pair_id}.bam | grep ^SN | cut -f 2- | awk '{ print "${fq_pair_id}\t" \$0 }' | sed 's/://g') > bam_stat
    """
}

process combine_fq_bam_stats {

    publishDir analysis_dir + "/fq", mode: 'copy'

    input:
        file("*.stat.txt") from fq_bam_stat_files.toSortedList()

    output:
        file("fq_bam_stats.tsv")

    """
        echo -e "fq_pair_id\\tvariable\\tvalue\\tcomment" > fq_bam_stats.tsv
        cat *.stat.txt >> fq_bam_stats.tsv
    """
}


/* 
  Merge - Generate SM Bam
*/

process merge_bam {

    cpus cores

    publishDir SM_alignments_dir + "/WI/strain", mode: 'copy', pattern: '*.bam*'

    tag { SM }

    input:
        set SM, bam, index from SM_aligned_bams.groupTuple()

    output:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") into bam_set
        file("${SM}.duplicates.txt") into duplicates_file
        
    """

    count=`echo ${bam.join(" ")} | tr ' ' '\\n' | wc -l`

    if [ "\${count}" -eq "1" ]; then
        ln -s ${bam.join(" ")} ${SM}.merged.bam
        ln -s ${bam.join(" ")}.bai ${SM}.merged.bam.bai
    else
        sambamba merge --nthreads=${cores} --show-progress ${SM}.merged.bam ${bam.sort().join(" ")}
        sambamba index --nthreads=${cores} ${SM}.merged.bam
    fi

    picard MarkDuplicates I=${SM}.merged.bam O=${SM}.bam M=${SM}.duplicates.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=false
    sambamba index --nthreads=${cores} ${SM}.bam
    """
}

bam_set.into { merged_bams_for_coverage; merged_bams_individual; merged_bams_union; bams_idxstats; bams_stats; fq_concordance_bam  }


/*
    SM_idx_stats
*/

process SM_idx_stats {
    
    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bams_idxstats
    output:
        file('bam_idxstats.txt') into bam_idxstats_set

    """
        samtools idxstats ${SM}.bam | awk '{ print "${SM}\\t" \$0 }' > bam_idxstats.txt
    """
}

process SM_combine_idx_stats {

    publishDir analysis_dir + "/SM", mode: 'copy'

    input:
        file("*.stat.txt") from bam_idxstats_set.toSortedList()

    output:
        file("SM_bam_idxstats.tsv")

    """
        echo -e "SM\\treference\\treference_length\\tmapped_reads\\tunmapped_reads" > SM_bam_idxstats.tsv
        cat *.stat.txt | sort >> SM_bam_idxstats.tsv
    """

}


/*
    SM bam stats
*/

process SM_bam_stats {

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bams_stats

    output:
        file('bam_stat.txt') into SM_bam_stat_files

    """
        cat <(samtools stats ${SM}.bam | grep ^SN | cut -f 2- | awk '{ print "${SM}\t" \$0 }' | sed 's/://g') > bam_stat.txt
    """
}

process combine_SM_bam_stats {

    publishDir analysis_dir + "/SM", mode: 'copy'

    input:
        file("?.stat.txt") from SM_bam_stat_files.toSortedList()

    output:
        file("SM_bam_stats.tsv")

    """
        echo -e "fq_pair_id\\tvariable\\tvalue\\tcomment" > SM_bam_stats.tsv
        cat *.stat.txt | sort >> SM_bam_stats.tsv
    """
}



process format_duplicates {

    publishDir analysis_dir + "/duplicates", mode: 'copy'

    input:
        val duplicates_set from duplicates_file.toSortedList()

    output:
        file("bam_duplicates.tsv")


    """
        echo -e 'filename\\tlibrary\\tunpaired_reads_examined\\tread_pairs_examined\\tsecondary_or_supplementary_rds\\tunmapped_reads\\tunpaired_read_duplicates\\tread_pair_duplicates\\tread_pair_optical_duplicates\\tpercent_duplication\\testimated_library_size' > bam_duplicates.tsv
        for i in ${duplicates_set.join(" ")}; do
            f=\$(basename \${i})
            cat \${i} | awk -v f=\${f/.duplicates.txt/} 'NR >= 8 && \$0 !~ "##.*" && \$0 != ""  { print f "\\t" \$0 } NR >= 8 && \$0 ~ "##.*" { exit }'  >> bam_duplicates.tsv
        done;
    """
}

/*
    Coverage Bam
*/
process coverage_SM {

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from merged_bams_for_coverage

    output:
        file("${SM}.coverage.tsv") into SM_coverage


    """
        bam coverage ${SM}.bam > ${SM}.coverage.tsv
    """
}


process coverage_SM_merge {

    publishDir analysis_dir + "/SM", mode: 'copy'


    input:
        val sm_set from SM_coverage.toSortedList()

    output:
        file("SM_coverage.full.tsv")
        file("SM_coverage.tsv") into SM_coverage_merged
        file("SM_coverage.tsv") into SM_coverage_network

    """
        echo -e 'bam\\tcontig\\tstart\\tend\\tproperty\\tvalue' > SM_coverage.full.tsv
        cat ${sm_set.join(" ")} >> SM_coverage.full.tsv

        # Generate condensed version
        cat <(echo -e 'strain\\tcoverage') <(cat SM_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6 | sort) > SM_coverage.tsv
    """

}

process call_variants_individual {

    cpus variant_cores

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from merged_bams_individual

    output:
        file("${SM}.individual.sites.tsv") into individual_sites

    """
    # Perform individual-level calling
    contigs="`samtools view -H ${SM}.bam | grep -Po 'SN:([^\\W]+)' | cut -c 4-40`"
    echo \${contigs} | tr ' ' '\\n' | xargs --verbose -I {} -P ${cores} sh -c "samtools mpileup --redo-BAQ -r {} --BCF --output-tags DP,AD,ADF,ADR,SP --fasta-ref ${reference} ${SM}.bam | bcftools call --skip-variants indels --variants-only --multiallelic-caller -O z  -  > ${SM}.{}.individual.vcf.gz"
    order=`echo \${contigs} | tr ' ' '\\n' | awk '{ print "${SM}." \$1 ".individual.vcf.gz" }'`
    
    # Output variant sites
    bcftools concat \${order} -O v | vk geno het-polarization - | bcftools view -O z > ${SM}.individual.vcf.gz
    bcftools index ${SM}.individual.vcf.gz
    rm \${order}

    bcftools view -M 2 -m 2 -O v ${SM}.individual.vcf.gz | \\
    bcftools filter --include 'DP > 3' | \\
    egrep '(^#|1/1)' | \\
    bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\n' > ${SM}.individual.sites.tsv
    """
}


process fq_concordance {

    cpus variant_cores

    tag { SM }

    input:
        set val(SM), file("input.bam"), file("input.bam.bai") from fq_concordance_bam

    output:
        file('out.tsv') into fq_concordance_out

    """
        # Split bam file into individual read groups; Ignore MtDNA
        contigs="`samtools view -H input.bam | grep -Po 'SN:([^\\W]+)' | cut -c 4-40 | grep -v 'MtDNA' | tr ' ' '\\n'`"
        samtools split -f '%!.%.' input.bam
        # DO NOT INDEX ORIGINAL BAM; ELIMINATES CACHE!
        bam_list="`ls -1 *.bam | grep -v 'input.bam'`"

        ls -1 *.bam | grep -v 'input.bam' | xargs --verbose -I {} -P ${variant_cores} sh -c "samtools index {}"

        # Generate a site list for the set of fastqs
        rg_list="`samtools view -H input.bam | grep '@RG.*ID:' | cut -f 2 | sed  's/ID://'`"
        # Perform individual-level calling
        for rg in \$rg_list; do
            echo \${contigs} | tr ' ' '\\n' | xargs --verbose -I {} -P ${variant_cores} sh -c "samtools mpileup --redo-BAQ -r {} --BCF --output-tags DP,AD,ADF,ADR,SP --fasta-ref ${reference} \${rg}.bam | bcftools call --skip-variants indels --variants-only --multiallelic-caller -O v | bcftools query -f '%CHROM\\t%POS\\n' >> {}.\${rg}.site_list.tsv"
        done;
        cat *.site_list.tsv  | sort --temporary-directory=${tmpdir} -k1,1 -k2,2n | uniq > site_list.srt.tsv
        bgzip site_list.srt.tsv -c > site_list.srt.tsv.gz && tabix -s1 -b2 -e2 site_list.srt.tsv.gz
        
        # Call a union set of variants
        for rg in \$rg_list; do
            echo \${contigs} | tr ' ' '\\n' | xargs --verbose -I {} -P ${variant_cores} sh -c "samtools mpileup --redo-BAQ -r {} --BCF --output-tags DP,AD,ADF,ADR,SP --fasta-ref ${reference} \${rg}.bam | bcftools call -T site_list.srt.tsv.gz --skip-variants indels --multiallelic-caller -O z > {}.\${rg}.vcf.gz"
            order=`echo \${contigs} | tr ' ' '\\n' | awk -v rg=\${rg} '{ print \$1 "." rg ".vcf.gz" }'`
            # Output variant sites
            bcftools concat \${order} -O v | \\
            vk geno het-polarization - | \\
            bcftools filter -O u --threads ${variant_cores} --set-GTs . --include "QUAL >= 10 || FORMAT/GT == '0/0'" |  \\
            bcftools filter -O u --threads ${variant_cores} --set-GTs . --include "FORMAT/DP > 3" | \\
            bcftools filter -O u --threads ${variant_cores} --set-GTs . --include "INFO/MQ > ${mq}" | \\
            bcftools filter -O u --threads ${variant_cores} --set-GTs . --include "(FORMAT/AD[1])/(FORMAT/DP) >= ${dv_dp} || FORMAT/GT == '0/0'" | \\
            bcftools query -f '%CHROM\\t%POS[\\t%GT\\t${SM}\\n]' | grep -v '0/1' | awk -v rg=\${rg} '{ print \$0 "\\t" rg }' > \${rg}.rg_gt.tsv
        done;
        cat *.rg_gt.tsv > rg_gt.tsv
        touch out.tsv
        Rscript --vanilla ${fq_concordance_script} 
    """
}

process combine_fq_concordance {

    publishDir analysis_dir + "/concordance", mode: 'copy'

    input:
        file("out*.tsv") from fq_concordance_out.toSortedList()

    output:
        file("fq_concordance.tsv")

    """
        cat <(echo 'a\tb\tconcordant_sites\ttotal_sites\tconcordance\tSM') out*.tsv > fq_concordance.tsv
    """


}


/*
    Merge individual sites
*/

process merge_variant_list {

    publishDir analysis_dir + "/sitelist", mode: 'copy'
    
    input:
        val sites from individual_sites.toSortedList()

    output:
        set file("sitelist.tsv.gz"), file("sitelist.tsv.gz.tbi") into gz_sitelist
        set file("sitelist.tsv.gz"), file("sitelist.tsv.gz.tbi") into fq_gz_sitelist
        file("sitelist.tsv") into sitelist
        file("sitelist.count.txt")


    """
        echo ${sites}
        cat ${sites.join(" ")} | sort -k1,1 -k2,2n | uniq > sitelist.tsv
        cat sitelist.tsv | wc -l > sitelist.count.txt
        bgzip sitelist.tsv -c > sitelist.tsv.gz && tabix -s1 -b2 -e2 sitelist.tsv.gz
    """
}


/* 
    Call variants using the merged site list
*/


union_vcf_channel = merged_bams_union.spread(gz_sitelist)


process call_variants_union {

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai"), file('sitelist.tsv.gz'), file('sitelist.tsv.gz.tbi') from union_vcf_channel

    output:
        file("${SM}.union.vcf.gz") into union_vcf_set

    """
        contigs="`samtools view -H ${SM}.bam | grep -Po 'SN:([^\\W]+)' | cut -c 4-40`"
        echo \${contigs} | tr ' ' '\\n' | xargs --verbose -I {} -P ${cores} sh -c "samtools mpileup --redo-BAQ -r {} --BCF --output-tags DP,AD,ADF,ADR,INFO/AD,SP --fasta-ref ${reference} ${SM}.bam | bcftools call -T sitelist.tsv.gz --skip-variants indels --multiallelic-caller -O z  -  > ${SM}.{}.union.vcf.gz"
        order=`echo \${contigs} | tr ' ' '\\n' | awk '{ print "${SM}." \$1 ".union.vcf.gz" }'`

        # Output variant sites
        bcftools concat \${order} -O v | \\
        vk geno het-polarization - | \\
        bcftools filter -O u --threads ${variant_cores} --set-GTs . --include "QUAL >= ${qual} || FORMAT/GT == '0/0'" |  \\
        bcftools filter -O u --threads ${variant_cores} --set-GTs . --include "FORMAT/DP > ${min_depth}" | \\
        bcftools filter -O u --threads ${variant_cores} --set-GTs . --include "INFO/MQ > ${mq}" | \\
        bcftools filter -O u --threads ${variant_cores} --set-GTs . --include "(FORMAT/AD[1])/(FORMAT/DP) >= ${dv_dp} || FORMAT/GT == '0/0'" | \\
        bcftools view -O z > ${SM}.union.vcf.gz
        bcftools index ${SM}.union.vcf.gz
        rm \${order}
    """

}



process generate_union_vcf_list {

    cpus 1 

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
       val vcf_set from union_vcf_set.toSortedList()

    output:
       file("union_vcfs.txt") into union_vcfs

    """
        echo ${vcf_set.join(" ")} | tr ' ' '\\n' > union_vcfs.txt
    """
}

union_vcfs_in = union_vcfs.spread(contigs)

process merge_union_vcf_chromosome {

    tag { chrom }

    input:
        set file(union_vcfs:"union_vcfs.txt"), val(chrom) from union_vcfs_in

    output:
        val(chrom) into contigs_list_in
        file("${chrom}.merged.raw.vcf.gz") into raw_vcf

    """
        bcftools merge --threads 10 --regions ${chrom} -O z -m all --file-list ${union_vcfs} > ${chrom}.merged.raw.vcf.gz
        bcftools index ${chrom}.merged.raw.vcf.gz
    """
}


// Generate a list of ordered files.
contig_raw_vcf = contig_list*.concat(".merged.raw.vcf.gz")

process concatenate_union_vcf {

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
        val merge_vcf from raw_vcf.toSortedList()

    output:
        set file("merged.raw.vcf.gz"), file("merged.raw.vcf.gz.csi") into raw_vcf_concatenated

    """
        for i in ${merge_vcf.join(" ")}; do
            ln  -s \${i} `basename \${i}`;
        done;
        chrom_set="";
        bcftools concat -O z ${contig_raw_vcf.join(" ")}  > merged.raw.vcf.gz
        bcftools index merged.raw.vcf.gz
    """
}

process filter_union_vcf {

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
        set file("merged.raw.vcf.gz"), file("merged.raw.vcf.gz.csi") from raw_vcf_concatenated

    output:
        set file("concordance.vcf.gz"), file("concordance.vcf.gz.csi") into filtered_vcf
        set file("concordance.vcf.gz"), file("concordance.vcf.gz.csi") into filtered_vcf_pairwise

    """
        bcftools view merged.raw.vcf.gz | \\
        vk filter ALT --max=0.99 - | \\
        vk filter MISSING --max=0.05 - | \\
        vk filter REF --min=1 - | \\
        vk filter ALT --min=1 - | \\
        bcftools view -O z - > concordance.vcf.gz
        bcftools index concordance.vcf.gz
    """
}

filtered_vcf.into { filtered_vcf_gtcheck; filtered_vcf_stat; filtered_vcf_phylo }


process calculate_gtcheck {

    publishDir analysis_dir + "/concordance", mode: 'copy'

    input:
        set file("concordance.vcf.gz"), file("concordance.vcf.gz.csi") from filtered_vcf_gtcheck

    output:
        file("gtcheck.tsv") into gtcheck
        file("gtcheck.tsv") into pairwise_compare_gtcheck

    """
        echo -e "discordance\\tsites\\tavg_min_depth\\ti\\tj" > gtcheck.tsv
        bcftools gtcheck -H -G 1 concordance.vcf.gz | egrep '^CN' | cut -f 2-6 >> gtcheck.tsv
    """

}


process stat_tsv {

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
        set file("concordance.vcf.gz"), file("concordance.vcf.gz.csi") from filtered_vcf_stat

    output:
        file("concordance.stats") into filtered_stats

    """
        bcftools stats --verbose concordance.vcf.gz > concordance.stats
    """

}

/*
    Perform concordance analysis
*/

concordance_script = Channel.fromPath("process_concordance.R")

process process_concordance_results {

    publishDir analysis_dir + "/concordance", mode: "copy"

    input:
        file 'gtcheck.tsv' from gtcheck
        file 'filtered.stats.txt' from filtered_stats
        file 'SM_coverage.tsv' from SM_coverage_merged

    output:
        file("concordance.svg")
        file("concordance.png")
        file("xconcordance.svg")
        file("xconcordance.png")
        file("isotype_groups.tsv") into isotype_groups
        file("isotype_count.txt")

    """
    # Run concordance analysis
    Rscript --vanilla ${process_concordance}
    """

}

process generate_isotype_groups {

    executor 'local'

    input:
        file("isotype_groups.tsv") from isotype_groups

    output:
        file("pairwise_groups.txt") into pairwise_groups

    """
    cat isotype_groups.tsv | awk '{ curr_strain = \$2; curr_group = \$1; if (group_prev == curr_group) { print prev_strain "," curr_strain "\t" \$1 "\t" \$3 } ; prev_strain = \$2; group_prev = \$1; }' > pairwise_groups.txt
    """

}

pairwise_groups_input = pairwise_groups.splitText( by:1 )

// Look for diverged regions among isotypes.
process pairwise_variant_compare {

    publishDir analysis_dir + "/pairwise", mode: 'copy'

    tag { pair }

    input:
        val(pair_group) from pairwise_groups_input
        set file("concordance.vcf.gz"), file("concordance.vcf.gz.csi") from filtered_vcf_pairwise.first() 

    output:
        file("${group}.${isotype}.${pair.replace(",","_")}.png")
        file("${group}.${isotype}.${pair.replace(",","_")}.tsv")

    script:
        pair_group = pair_group.trim().split("\t")
        pair = pair_group[0]
        group = pair_group[1]
        isotype = pair_group[2]

    """
        bcftools view -s ${pair} concordance.vcf.gz | grep '0/0' | grep '1/1' | cut -f 1,2 > out.tsv
        Rscript --vanilla ${plot_pairwise_script}
        mv out.png ${group}.${isotype}.${pair.replace(",","_")}.png
        mv out.tsv ${group}.${isotype}.${pair.replace(",","_")}.tsv
    """

}

filtered_vcf_phylo_contig = filtered_vcf_phylo.spread(["I", "II", "III", "IV", "V", "X", "MtDNA", "genome"])

/*
    Phylo analysis
*/

process phylo_analysis {

    publishDir analysis_dir + "/phylo", mode: "copy"

    tag { contig }

    input:
        set file("merged.filtered.vcf.gz"), file("merged.filtered.vcf.gz.csi"), val(contig) from filtered_vcf_phylo_contig

    output:
        set val(contig), file("${contig}.tree") into trees

    """
        # generate tree
        if [ "${contig}" == "genome" ]
        then
            vk phylo tree nj merged.filtered.vcf.gz > genome.tree
        else
            vk phylo tree nj merged.filtered.vcf.gz ${contig} > ${contig}.tree
        fi
    """
}

trees_phylo = trees.spread(Channel.fromPath("process_trees.R"))

process plot_trees {

    publishDir analysis_dir + "/phylo", mode: "copy"

    tag { contig }

    input:
        set val(contig), file("${contig}.tree"), file("process_trees.R") from trees_phylo

    output:
        file("${contig}.svg")
        file("${contig}.png")

    """
    Rscript --vanilla process_trees.R ${contig}
    """

}

