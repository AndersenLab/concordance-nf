# concordance-nf

Perform Fastq-Profiling in the current working directory.

### Usage

The `concordance-nf` pipeline takes as input a JSON file. The JSON file is organized as a hierarchy as follows:

```
{
    "<strain_name>": {
        "@RG\tID:<sample_id>LB:<sequencing_library>SM:<strain_name>": [
            "/path/to/fq1.fq.gz",
            "/path/to/fq2.fq.gz",
            "<sample_id>"
        ]
}
```

The following variables are defined above:

* `<sample_id>` - A unique id for the sample. This should be unique for every sequencing run of a given sample.
* `<sequencing_library>` - The sequencing library for a sample. In general, the Nextera Barcodes used can be substituted for the sequencing library.
* `<sample_name>` - `SM` stands for sample, but here we are using it to reflect the strain name.

A complete example is below:
```
{
  "AB1": {
    "@RG\tID:original_wi_seq____BGI1-RET2-AB1-trim-\tLB:RET2\tSM:AB1": [
      "/Users/dancook/Documents/quest/concordance-nf/test_data/BGI1-RET2-AB1-trim-1P.fq.gz",
      "/Users/dancook/Documents/quest/concordance-nf/test_data/BGI1-RET2-AB1-trim-2P.fq.gz",
      "BGI1-RET2-AB1-trim-"
    ],
    "@RG\tID:original_wi_seq____BGI2-RET2-AB1-trim-\tLB:RET2\tSM:AB1": [
      "/Users/dancook/Documents/quest/concordance-nf/test_data/BGI2-RET2-AB1-trim-1P.fq.gz",
      "/Users/dancook/Documents/quest/concordance-nf/test_data/BGI2-RET2-AB1-trim-2P.fq.gz",
      "BGI2-RET2-AB1-trim-"
    ],
    "@RG\tID:original_wi_seq____BGI3-RET2b-AB1-trim-\tLB:RET2b\tSM:AB1": [
      "/Users/dancook/Documents/quest/concordance-nf/test_data/BGI3-RET2b-AB1-trim-1P.fq.gz",
      "/Users/dancook/Documents/quest/concordance-nf/test_data/BGI3-RET2b-AB1-trim-2P.fq.gz",
      "BGI3-RET2b-AB1-trim-"
    ],
    "@RG\tID:original_wi_seq____MMP-ML_2-AB1-trim-\tLB:ML_2\tSM:AB1": [
      "/Users/dancook/Documents/quest/concordance-nf/test_data/MMP-ML_2-AB1-trim-1P.fq.gz",
      "/Users/dancook/Documents/quest/concordance-nf/test_data/MMP-ML_2-AB1-trim-2P.fq.gz",
      "MMP-ML_2-AB1-trim-"
    ]
  },
  "AB2": {
    "@RG\tID:original_wi_seq____Rockman-R999.1-AB2-trim-\tLB:R999.1\tSM:AB2": [
      "/Users/dancook/Documents/quest/concordance-nf/test_data/R999.1-AB2-trim-1P.fq.gz",
      "/Users/dancook/Documents/quest/concordance-nf/test_data/R999.1-AB2-trim-2P.fq.gz",
      "Rockman-R999.1-AB2-trim-"
    ]
  }
```

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

### Results

The pipeline results w
