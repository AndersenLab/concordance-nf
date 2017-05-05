# concordance-nf

Perform Fastq-Profiling in the current working directory.

### Usage
```
nextflow run main.nf -resume
```

### Configuration

The following variables should be set in either `~/.nextflow/config` OR `nextflow.config`

A temporary directory:
```
tmpdir = "/projects/b1042/AndersenLab/tmp`"
```

A genome name and reference file. The reference must be indexed with bwa.
```
genome = "WS245"
reference = "/projects/b1059/data/genomes/c_elegans/${genome}/${genome}.fa.gz"
```

__Cores/threads for different processes__
```
alignment_cores = 16
variant_cores = 6
compression_threads = 4
```

__Output directories__

The `analysis_dir` will output a structured directory of results that are detailed below. The `SM_alignments_dir` will output bams for each strain.

```
analysis_dir = "/projects/b1059/analysis/WI-concordance" # For output of results.
SM_alignments_dir = "/projects/b1059/data/alignments/WI/SM" # For sample level alignments.
```

Process level options can also be set:
```
process {
    module='gcc/5.1.0:R/3.3.1'
    $perform_alignment {
        cpus = 4
    }
    $call_variants_individual {
        cpus = 6
        memory = '8G'
    }
    $call_variants_union {
        cpus = 6
        memory = '8G'
    }
    $merge_union_vcf {
        cpus = 20
        memory = '50G'
    }
    $filter_union_vcf {
        cpus = 20
        memory = '50G'
    }

}
```

### Results

The pipeline results are output in the following structure.

```
├── strain_set.json # A list of strains used in the analysis (json format).
├── concordance
│   ├── concordance.png                      # Histogram of concordance, colored by isotype.
│   ├── concordance.svg
│   ├── fq_concordance.tsv -   -  -  -  -  - # Concordance numbers at a fastq-level
│   ├── gtcheck.tsv                          # output of bcftools gtcheck 
│   ├── isotype_groups.tsv -   -  -  -  -  - # Generated isotype groupings.
│   ├── problem_SM.tsv                       # Problematic groupings, if any.
│   ├── xconcordance.png -  -  -  -  -  -  - # Zoomed in view of concordance.
│   └── xconcordance.svg
├── duplicates
│   └── bam_duplicates.tsv                   # Summary of duplicate reads (determined by picard).
├── phylo
│   ├── genome.png                           # Genome phylogeny
│   ├── genome.svg
│   ├── genome.tree                          # Newick format of tree
│   ├── I.png # Chromosome I phylogeny
│   ├── I.svg
│   ├── I.tree
│   ├── ...
│   ├── MtDNA.png
│   ├── MtDNA.svg
│   ├── MtDNA.tree
├── sitelist
│   ├── sitelist.count.txt                   # Number of sites called.
│   ├── sitelist.tsv                         # List of sites used in genotyping.
│   ├── sitelist.tsv.gz                      # Compressed and indexed (with .tbi) list of sites.
│   └── sitelist.tsv.gz.tbi                
├── fq
│   ├── fq_bam_idxstats.tsv   -  -  -  -  -  # Stats generated with `samtools idxstats`
│   ├── fq_bam_stats.tsv                     # Stats generated with `samtools stats`
│   ├── fq_coverage.full.tsv  -  -  -  -  -  # Detailed coverage numbers
│   └── fq_coverage.tsv                      # Summary coverage numbers.
├── SM
│   ├── SM_bam_idxstats.tsv                  # Strain-level statistics; Same as fq descriptions above.
│   ├── SM_bam_stats.tsv
│   ├── SM_coverage.full.tsv
│   └── SM_coverage.tsv
└── vcf
    ├── concordance.vcf.gz        # Filtered VCF, filtered for true SNPs (no monomorphic sites)
    ├── concordance.vcf.gz.csi    # Concordance VCF Index
    ├── concordance.stats         # Stats from concordance vcf. Contains unumber of SNPs
    └── union_vcfs.txt            # List of VCFs that are combined for performing concordance.
```

### Details on files above

Below are details regarding some of the generated files above:

#### concordance/

__isotype_groups.tsv__

|   group | strain   | isotype   |   coverage |   unique_isotypes_per_group | strain_conflicts   |
|--------:|:---------|:----------|-----------:|----------------------------:|:-------------------|
|       1 | AB1      | AB1       |    69.4687 |                           1 | FALSE              |
|     112 | AB4      | CB4858    |   158.358  |                           1 | FALSE              |
|     112 | ECA251   | CB4858    |    73.5843 |                           1 | FALSE              |
|     112 | JU1960   | NA        |    55.0373 |                           1 | FALSE              |
|     159 | BRC20067 | BRC20067  |    33.5934 |                           1 | FALSE              |
|     159 | BRC20113 | NA        |    38.9916 |                           1 | FALSE              |

The `isotype_groups.tsv` classifies the strains into isotypes. There may exist issues, however.

