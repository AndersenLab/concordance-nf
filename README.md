# concordance-nf

Perform Fastq-Profiling in the current working directory.

### Usage

```
# cd to directory of fastqs
nextflow run Andersenlab/concordance-nf
```
### Configuration

The following variables should be set in either `~/.nextflow/config` OR `nextflow.config`

A temporary directory:
tmpdir = "/projects/b1042/AndersenLab/tmp`"
genome = "WS245"
reference = "/projects/b1059/data/genomes/c_elegans/${genome}/${genome}.fa.gz"
cores = 16 # For numb
compression_threads = 4
analysis_dir = "/projects/b1059/analysis/WI-concordance"
data_dir = "/projects/b1059/data"

email="Danielecook@gmail.com"

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
