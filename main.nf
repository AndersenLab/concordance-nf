#!/usr/bin/env nextflow
/*
 * Authors:
 * - Daniel Cook <danielecook@gmail.com>
 * - Ye Wang <yewangfaith@gmail.com>
 * - Dan Lu <dan.lu@northwestern.edu>
 */

nextflow.enable.dsl=2
// For now, this pipeline requires NXF_VER >= 23.0

/* 
    ======
    Params
    ======
*/

date = new Date().format( 'yyyyMMdd' )
params.out = "concordance-${date}"
params.species = "c_elegans"


// Debug
if (params.debug == true) {
    println """

        *** Using debug mode ***

    """
    params.vcf = "${workflow.projectDir}/test_data/concordance.test.vcf.gz"
    params.bam_coverage = "${workflow.projectDir}/test_data/SM_coverage.tsv"
    params.concordance_cutoff = 0.99
    params.species = "c_elegans"

} else {
    params.vcf = "(required)"
    params.bam_coverage = "(required)"
    params.concordance_cutoff = 0.9995
}

def log_summary() {

out = '''

┌─┐┌─┐┌┐┌┌─┐┌─┐┬─┐┌┬┐┌─┐┌┐┌┌─┐┌─┐  ┌┐┌┌─┐
│  │ │││││  │ │├┬┘ ││├─┤││││  ├┤───│││├┤ 
└─┘└─┘┘└┘└─┘└─┘┴└──┴┘┴ ┴┘└┘└─┘└─┘  ┘└┘└  
                                                         
''' + """

nextflow main.nf -profile quest --debug=true

nextflow main.nf -profile quest --vcf=a.vcf.gz --bam_coverage=mqc_mosdepth-coverage-per-contig_1.txt --info_sheet=WI_elegans.tsv

    parameters           description                                                               Set/Default
    ==========           ===========                                                               =======

    --debug              Set to 'true' to test                                                     ${params.debug}
    --out                Directory to output results                                               ${params.out}
    --vcf                Hard filtered vcf                                                         ${params.vcf}
    --bam_coverage       Table with "strain" and "coverage" as header                              ${params.bam_coverage}
    --species            'c_elegans' will check for npr1. All other values will skip this          ${params.species}
    --concordance_cutoff Cutoff of concordance value to count two strains as same isotype          ${params.concordance_cutoff}

    HELP: http://andersenlab.org/dry-guide/pipeline-concordance/
    ----------------------------------------------------------------------------------------------------------------------------
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]

"""
out
}


log.info(log_summary())

if (params.help) {
    exit 1
}



workflow {

    get_species_sheet()

    hard_filtered_vcf = Channel.fromPath("${params.vcf}")
    vcf_index = Channel.fromPath("${params.vcf}.tbi") // make sure the index format is consistent in process inputs
    bam_coverage = Channel.fromPath("${params.bam_coverage}")

    hard_filtered_vcf.combine(vcf_index) | (calculate_gtcheck)

// moved npr1 allele check to alignment
// if (params.species == "c_elegans") {

//     hard_filtered_vcf.combine(vcf_index) | (calculate_gtcheck & npr1_allele_check)

// } else {

//     hard_filtered_vcf.combine(vcf_index) | (calculate_gtcheck)

// }


    calculate_gtcheck.out.combine(bam_coverage)
        .combine(get_species_sheet.out) | process_concordance_results

    process_concordance_results.out.isotype_groups_ch | generate_isotype_groups

    generate_isotype_groups.out.splitText( by:1 ).combine(hard_filtered_vcf).combine(vcf_index) | within_group_pairwise // this is for strains of same isotype


//    hard_filtered_vcf.combine(vcf_index) | (query_between_group_pairwise_gt & strain_pairwise_list_N2)

//    strain_pairwise_list_N2.out.splitText( by:1 ).combine(query_between_group_pairwise_gt.out) | between_group_pairwise // this is for each pair of strains

//    between_group_pairwise.out.cutoff_distribution | cutoff_distribution

//    between_group_pairwise.out.between_group_pairwise_out.toSortedList() | merge_betweengroup_pairwise_output


//    process_concordance_results.out.isotype_groups_ch.combine(merge_betweengroup_pairwise_output.out).combine(npr1_allele_check.out) | combine_pairwise_results

}


process get_species_sheet {
    
    label 'xs'

    publishDir "${params.out}", mode: 'copy'

    output:
        file("*.tsv")

    """
    Rscript --vanilla ${workflow.projectDir}/bin/download_google_sheet.R ${params.species}
        
    """

}


process calculate_gtcheck {

    label 'xs'

    publishDir "${params.out}/concordance", mode: 'copy'

    input:
        tuple file("concordance.vcf.gz"), file("concordance.vcf.gz.tbi")

    output:
        file("gtcheck.tsv")

    """
        echo -e "discordance\\tsites\\tavg_min_depth\\ti\\tj" > gtcheck.tsv
        bcftools view -i 'TYPE="snp"' -O u concordance.vcf.gz | bcftools gtcheck -H -G 1 | egrep '^CN' | cut -f 2-6 >> gtcheck.tsv
    """
}



process process_concordance_results {

    label 'xs'

    publishDir "${params.out}/concordance", mode: "copy"

    input:
        tuple file("gtcheck.tsv"), file("SM_coverage.tsv"), file("WI_info_sheet.tsv")

    output:
        file("concordance.pdf")
        file("concordance.png")
        file("xconcordance.pdf")
        file("xconcordance.png")
        path "isotype_groups.tsv"
        path "isotype_groups_new.tsv", emit: isotype_groups_ch
        file("isotype_count.txt")
        file("WI_metadata.tsv")
        file("problem_strains.tsv")

    """
    # Run concordance analysis
    Rscript --vanilla ${workflow.projectDir}/bin/process_concordance.R SM_coverage.tsv WI_info_sheet.tsv ${params.concordance_cutoff}
    """
}


// update 20220202 - do we need ALL groups or just new ones?
process generate_isotype_groups {

    executor 'local'
    container null

    input:
        file("isotype_groups_new.tsv") //from isotype_groups

    output:
        file("pairwise_groups.txt") //into pairwise_groups

    """
    cat isotype_groups_new.tsv | awk '{ curr_strain = \$2; curr_group = \$1; if (group_prev == curr_group) { print prev_strain "," curr_strain "\t" \$1 "\t" \$3 } ; prev_strain = \$2; group_prev = \$1; }' > pairwise_groups.txt
    """

}


process within_group_pairwise {

    label 'xs'

    publishDir "${params.out}/concordance/pairwise/within_group", mode: 'copy', overwrite: true

    tag { pair }

    input:
        tuple val(pair_group), file("concordance.vcf.gz"), file("concordance.vcf.gz.tbi") //from filtered_vcf_pairwise

    output:
        file("*.png")

    script:
        pair_group = pair_group.trim().split("\t")
        pair = pair_group[0]
        group = pair_group[1]
        isotype = pair_group[2]

    """
        bcftools query -f '%CHROM\t%POS[\t%GT]\n' -s ${pair} concordance.vcf.gz > out.tsv
        Rscript --vanilla ${workflow.projectDir}/bin/plot_pairwise.R ${pair} ${group} ${isotype}

    """
}


// The belows are new processes by Ye for futher checking


process strain_pairwise_list_N2 {

    executor 'local'

//    publishDir "${params.out}/concordance/pairwise/between_strains", mode: "copy"

    input:
        tuple file("concordance.vcf.gz"), file("concordance.vcf.gz.tbi")

    output:
        path("strain_pairwise_list.tsv")//into strain_pairwise

    """
        # generate strain level pairwise comparison list
        bcftools query -l concordance.vcf.gz | grep -v "^N2\$" > raw_strain.tsv

        for i in `cat raw_strain.tsv` ; do
            echo \${i}-N2
        done | tr "-" "\t" > strain_pairwise_list.tsv
    """
}


process query_between_group_pairwise_gt {

    publishDir "${params.out}/variation", mode: 'copy', overwrite: true

    label 'md'

    input:
        tuple file("concordance.vcf.gz"), file("concordance.vcf.gz.tbi") 

    output:
        file("out_gt.tsv") 

    """
        # query the genotype for pairwise comparison
        bcftools query -f '%CHROM\\t%POS[\\t%GT]\\n' concordance.vcf.gz -H | sed 's/\\:GT//g' | sed 's/\\[[0-9]*\\]//g' | sed 's/\\#//g' | csvtk rename -t -f 1 -n CHROM > out_gt.tsv
    """
}

process between_group_pairwise {

    publishDir "${params.out}/concordance/pairwise/between_group", mode: 'copy', overwrite: true, pattern: '*.png'

    tag "${sp1}_${sp2}"

    label 'xs'

    input:
        tuple val(pair_group), file("out_gt.tsv") 

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

// Moved this process to alignment now
// process npr1_allele_check {

//     cpus params.cores

//     publishDir "${params.out}/concordance", mode: 'copy', overwrite: true

//     input:
//         tuple file("concordance.vcf.gz"), file("concordance.vcf.gz.tbi") //from npr1_allele

//     output:
//         file("npr1_allele_strain.tsv") //into npr1_out

//     """
//         echo -e 'problematic_strain\\tgt' > npr1_allele_strain.tsv
//         bcftools view --threads ${params.cores} -t X:4768788 concordance.vcf.gz | bcftools query -f '[%SAMPLE\\t%GT\\n]' | awk '\$2 != "1/1"' >> npr1_allele_strain.tsv
//     """
// }

process merge_betweengroup_pairwise_output {

    executor 'local'
    container null

    publishDir "${params.out}/concordance", mode: 'copy', overwrite: true

    input:
        file(bg_pairwise) 

    output:
        file("merge_betweengroup_pairwise_output.tsv")


    """
        echo ${bg_pairwise}
        echo -e 'pairwise\\tconcordant_bin_gt_70\\tmax_discordant_bin_count_lt_3\\tmean_discordant_bin_count_lt_2.5\\tno_bin_lt_0.9\\tsuspected_introgress' > merge_betweengroup_pairwise_output.tsv
        cat ${bg_pairwise.join(" ")} | cut -f 1,2,3,4,5,6 >> merge_betweengroup_pairwise_output.tsv
    """
}

process cutoff_distribution {

    executor 'local'
    container null

    publishDir "${params.out}/concordance", mode: 'copy', overwrite: true

    input:
        file(cutoff_val) 

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

    label 'xs'

    publishDir "${params.out}/concordance", mode: 'copy', overwrite: true

    input:
        tuple file("isotype_groups.tsv"), file("merge_betweengroup_pairwise_output.tsv"), file("npr1_allele_strain.tsv")

    output:
        file("new_isotype_groups.tsv")

    """
        Rscript --vanilla ${workflow.projectDir}/bin/merge_groups_info.R isotype_groups.tsv merge_betweengroup_pairwise_output.tsv npr1_allele_strain.tsv
    """
}


