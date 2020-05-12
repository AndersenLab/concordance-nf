#!/usr/bin/env nextflow
/*
 * Authors:
 * - Daniel Cook <danielecook@gmail.com>
 * - Ye Wang <yewangfaith@gmail.com>
 * - Dan Lu <dan.lu@northwestern.edu>
 */

nextflow.preview.dsl=2
// For now, this pipeline requires NXF_VER >= 20.01.0

/*
    Params
*/

// Define contigs here!
CONTIG_LIST = ["I", "II", "III", "IV", "V", "X", "MtDNA"]
contig_list = Channel.from(CONTIG_LIST)

/* 
    ======
    Params
    ======
*/

date = new Date().format( 'yyyyMMdd' )
params.out = "concordance-${date}"
params.debug = false


//reference_handle = file("${params.reference}/*")

//fq_concordance_script = file("fq_concordance.R")

// Debug
if (params.debug == true) {
    println """

        *** Using debug mode ***

    """
    params.vcf = "${workflow.projectDir}/test_data/concordance_chrI.hard-filter.vcf.gz"
    params.bam_coverage = "${workflow.projectDir}/test_data/SM_coverage.tsv"

} else {
    params.vcf = "(required)"
    params.bam_coverage = "(required)"
}



/* 
    ==
    UX
    ==
*/

param_summary = '''

┌─┐┌─┐┌┐┌┌─┐┌─┐┬─┐┌┬┐┌─┐┌┐┌┌─┐┌─┐  ┌┐┌┌─┐
│  │ │││││  │ │├┬┘ ││├─┤││││  ├┤───│││├┤ 
└─┘└─┘┘└┘└─┘└─┘┴└──┴┘┴ ┴┘└┘└─┘└─┘  ┘└┘└  
                                                         
''' + """

    parameters          description                                                               Set/Default
    ==========          ===========                                                               =======

    --debug             Set to 'true' to test                                                     ${params.debug}
    --cores             Regular job cores                                                         ${params.cores}
    --out               Directory to output results                                               ${params.out}
    --vcf               Hard filtered vcf                                                         ${params.vcf}
    --bam_coverage      Strain level mqc_mosdepth-coverage-per-contig_1.txt from alignment-nf     ${params.bam_coverage}

    HELP: http://andersenlab.org/dry-guide/pipeline-concordance/

"""


workflow {

    hard_filtered_vcf = Channel.fromPath("${params.vcf}")
    vcf_index = Channel.fromPath("${params.vcf}.tbi")
    bam_coverage = Channel.fromPath("${params.bam_coverage}")


    hard_filtered_vcf.combine(vcf_index) | (calculate_gtcheck & heterozygosity_check & query_between_group_pairwise_gt & npr1_allele_check )

    calculate_gtcheck.out.combine(bam_coverage) | process_concordance_results

    process_concordance_results.out.isotype_groups_ch | generate_isotype_groups

    generate_isotype_groups.out.splitText( by:1 ).combine(hard_filtered_vcf).combine(vcf_index) | pairwise_variant_compare // this is for strains of same isotype

    bam_coverage | strain_pairwise_list

    strain_pairwise_list.out.splitText( by:1 ).combine(query_between_group_pairwise_gt.out) | between_group_pairwise 

    between_group_pairwise.out.cutoff_distribution | cutoff_distribution

    between_group_pairwise.out.between_group_pairwise_out.toSortedList() | merge_betweengroup_pairwise_output


    process_concordance_results.out.isotype_groups_ch.combine(merge_betweengroup_pairwise_output.out).combine(npr1_allele_check.out) | combine_pairwise_results

}



process calculate_gtcheck {

    publishDir "${params.out}/concordance", mode: 'copy'

    input:
        tuple file("concordance.vcf.gz"), file("concordance.vcf.gz.csi")

    output:
        file("gtcheck.tsv")

    """
        echo -e "discordance\\tsites\\tavg_min_depth\\ti\\tj" > gtcheck.tsv
        bcftools gtcheck -H -G 1 concordance.vcf.gz | egrep '^CN' | cut -f 2-6 >> gtcheck.tsv
    """
}
/*
process stat_tsv {

    publishDir "${params.out}/vcf", mode: 'copy'

    input:
        set file("concordance.vcf.gz"), file("concordance.vcf.gz.csi") from filtered_vcf_stat

    output:
        file("concordance.stats") into filtered_stats

    """
        bcftools stats --verbose concordance.vcf.gz > concordance.stats
    """
}
*/
process process_concordance_results {

    publishDir "${params.out}/concordance", mode: "copy"

    input:
        tuple file("gtcheck.tsv"), file("SM_coverage.tsv") // from gtcheck
//        file 'filtered.stats.txt' from filtered_stats
//        file 'SM_coverage.tsv' from for_concordance

    output:
        file("concordance.pdf")
        file("concordance.png")
        file("xconcordance.pdf")
        file("xconcordance.png")
        path "isotype_groups.tsv", emit: isotype_groups_ch
        file("isotype_count.txt")
        file("WI_metadata.tsv")
        file("problem_strains.tsv")

    """
    # Run concordance analysis
    Rscript --vanilla ${workflow.projectDir}/bin/process_concordance.R SM_coverage.tsv 
    """
}

//isotype_groups_ch.into { isotype_groups; for_combined_final}

process generate_isotype_groups {

    input:
        file("isotype_groups.tsv") //from isotype_groups

    output:
        file("pairwise_groups.txt") //into pairwise_groups

    """
    cat isotype_groups.tsv | awk '{ curr_strain = \$2; curr_group = \$1; if (group_prev == curr_group) { print prev_strain "," curr_strain "\t" \$1 "\t" \$3 } ; prev_strain = \$2; group_prev = \$1; }' > pairwise_groups.txt
    """

}

//pairwise_groups_input = pairwise_groups.splitText( by:1 )

// Look for diverged regions among isotypes.
process pairwise_variant_compare {

    publishDir "${params.out}/concordance/pairwise/within_group", mode: 'copy', overwrite: true

    tag { pair }

    input:
        tuple val(pair_group), file("concordance.vcf.gz"), file("concordance.vcf.gz.csi") //from filtered_vcf_pairwise

    output:
        file("${group}.${isotype}.${pair.replace(",","_")}.png")
        file("${group}.${isotype}.${pair.replace(",","_")}.tsv")

    script:
        pair_group = pair_group.trim().split("\t")
        pair = pair_group[0]
        group = pair_group[1]
        isotype = pair_group[2]

    """
        bcftools query -f '%CHROM\t%POS[\t%GT]\n' -s ${pair} concordance.vcf.gz > out.tsv
        Rscript --vanilla ${workflow.projectDir}/bin/plot_pairwise.R
        mv out.png ${group}.${isotype}.${pair.replace(",","_")}.png
        mv out.tsv ${group}.${isotype}.${pair.replace(",","_")}.tsv
    """
}

process heterozygosity_check {

    cpus params.cores

    publishDir "${params.out}/concordance", mode: "copy"

    input:
        set file("concordance.vcf.gz"), file("concordance.vcf.gz.csi") //from het_check_vcf

    output:
        file("heterozygosity.tsv")

    """
        bcftools query -l concordance.vcf.gz | xargs --verbose -I {} -P ${task.cpus} sh -c "bcftools query -f '[%SAMPLE\t%GT\n]' --samples={} concordance.vcf.gz | grep '0/1' | uniq -c >> heterozygosity.tsv"
    """

}

// The belows are new processes for futher checking

process strain_pairwise_list {

    executor 'local'

    publishDir "${params.out}/concordance/pairwise/between_strains", mode: "copy"

    input:
        file("SM_coverage.tsv")// from for_strain_list

    output:
        path("strain_pairwise_list.tsv")//into strain_pairwise

    """
        # generate strain level pairwise comparison list
        cat SM_coverage.tsv | cut -f1 | sed '1d' | grep -v "^N2\$" > raw_strain.tsv

        for i in `cat raw_strain.tsv` ; do
            echo \${i}-N2
        done | tr "-" "\t" > strain_pairwise_list.tsv
    """
}

//strain_pairwise.splitText( by:1 )
//               .set { new_strain_pairwise }


process query_between_group_pairwise_gt {

    publishDir "${params.out}/variation", mode: 'copy', overwrite: true

    cpus params.cores

    input:
        set file("concordance.vcf.gz"), file("concordance.vcf.gz.csi") //from strain_pairwise_vcf

    output:
        file("out_gt.tsv") //into gt_pairwise

    """
        # query the genotype for pairwise comparison
        bcftools query -f '%CHROM\\t%POS[\\t%GT]\\n' concordance.vcf.gz -H | sed 's/\\:GT//g' | sed 's/\\[[0-9]*\\]//g' | sed 's/\\#//g' | csvtk rename -t -f 1 -n CHROM > out_gt.tsv
    """
}

process between_group_pairwise {

    publishDir "${params.out}/concordance/pairwise/between_group", mode: 'copy', overwrite: true, pattern: '*.png'

    tag "${sp1}_${sp2}"

    input:
        tuple val(pair_group), file("out_gt.tsv") //from gt_pairwise

    output:
        path "${sp1}-${sp2}.tsv", emit: between_group_pairwise_out
        path "${sp1}-${sp2}-distribution.tsv", emit: cutoff_distribution
        file("${sp1}-${sp2}.disconcordance.png") optional true
        file("${sp1}-${sp2}.hist.png") optional true
    
    script:
        pair_group = pair_group.trim().split("\t")
        sp1 = pair_group[0]
        sp2 = pair_group[1]

    """            
        csvtk cut -t -f CHROM,POS,${sp1},${sp2} out_gt.tsv > ${sp1}-${sp2}.queried.tsv
        Rscript --vanilla ${workflow.projectDir}/bin/process_strain_pairwise.R ${sp1} ${sp2} ${sp1}-${sp2}.queried.tsv
        mv condition_results.tsv ${sp1}-${sp2}.tsv
        mv for_distribution.tsv ${sp1}-${sp2}-distribution.tsv
        rm ${sp1}-${sp2}.queried.tsv
    """
}

process npr1_allele_check {

    cpus params.cores

    publishDir "${params.out}/concordance", mode: 'copy', overwrite: true

    input:
        set file("concordance.vcf.gz"), file("concordance.vcf.gz.csi") //from npr1_allele

    output:
        file("npr1_allele_strain.tsv") //into npr1_out

    """
        echo -e 'problematic_strain\\tgt' > npr1_allele_strain.tsv
        bcftools view --threads ${params.cores} -t X:4768788 concordance.vcf.gz | bcftools query -f '[%SAMPLE\\t%GT\\n]' | awk '\$2 != "1/1"' >> npr1_allele_strain.tsv
    """
}

process merge_betweengroup_pairwise_output {

    publishDir "${params.out}/concordance", mode: 'copy', overwrite: true

    input:
        file(bg_pairwise) //from between_group_pairwise_out.toSortedList()

    output:
        file("merge_betweengroup_pairwise_output.tsv")// into combine_pairwise_results_ch


    """
        echo ${bg_pairwise}
        echo -e 'pairwise\\tconcordant_bin_gt_70\\tmax_discordant_bin_count_lt_3\\tmean_discordant_bin_count_lt_2.5\\tno_bin_lt_0.9\\tsuspected_introgress' > merge_betweengroup_pairwise_output.tsv
        cat ${bg_pairwise.join(" ")} | cut -f 1,2,3,4,5,6 >> merge_betweengroup_pairwise_output.tsv
    """
}

process cutoff_distribution {

    publishDir "${params.out}/concordance", mode: 'copy', overwrite: true

    input:
        file(cutoff_val) //from cutoff_distribution.toSortedList()

    output:
        file("cutoff_distribution.tsv")


    """
        echo ${cutoff_val}
        ls -al 1>&2
        echo -e 'pairwise\\tconcordant_bin_gt_70\\tmax_discordant_bin_count_lt_3\\tmean_discordant_bin_count_lt_2.5\\tmin_bin' > cutoff_distribution.tsv
        cat ${cutoff_val.join(" ")} | cut -f 1,2,3,4,5 >> cutoff_distribution.tsv
    """
}

process combine_pairwise_results {

    publishDir "${params.out}/concordance", mode: 'copy', overwrite: true

    input:
        tuple file("isotype_groups.tsv"), file("merge_betweengroup_pairwise_output.tsv"), file("npr1_allele_strain.tsv")

    output:
        file("new_isotype_groups.tsv")

    """
        Rscript --vanilla ${workflow.projectDir}/bin/merge_groups_info.R isotype_groups.tsv merge_betweengroup_pairwise_output.tsv npr1_allele_strain.tsv
    """
}


