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