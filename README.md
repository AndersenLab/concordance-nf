![Build Docker (env/concordance.Dockerfile)](https://github.com/AndersenLab/concordance-nf/workflows/Build%20Docker%20(env/concordance.Dockerfile)/badge.svg)

# concordance-nf

Group isotypes based on strain-level vcf.

### Usage
```
nextflow main.nf -profile quest --debug=true

nextflow main.nf -profile quest --vcf=a.vcf.gz --bam_coverage=mqc_mosdepth-coverage-per-contig_1.txt
```

### Parameters
```
                     details                                                                   default
--debug              Set to 'true' to test                                                     false
--cores              Regular job cores                                                         see config
--out                Directory to output results                                               concordance-date
--vcf                Hard filtered vcf                                                         NA, required
--bam_coverage       Table with "strain" and "coverage" as header                              NA, required
--info_sheet         Strain sheet containing exisiting isotype assignment                      NA, required
--species            'ce' will check for npr1. All other values will skip this                 empty
--concordance_cutoff Cutoff of concordance value to count two strains as same isotype          0.9995

```

### Results and notes
Note all files are generated by default. 
```
├── concordance
│   │
│   ├── gtcheck.tsv              # The raw pairwise comparison of variants in strains, used by /bin/process_concordance.R to define isotypes.
│   │                            # IMPORTANT: one should manually run process_concordance.R to inspect isotype assignment and adjust concordance cutoff.
│   │ 
│   │   # Output of process_concordance.R with provided concordance cutoff. 
│   ├── isotype_count.txt        # Number of isotypes.
│   ├── isotype_groups.tsv       # Isotype assignment incorporating the existing isotypes.
│   ├── problem_strains.tsv      # Strains with potential issues to look into.
│   ├── WI_metadata.tsv          # A copy of the --info_sheet.
│   ├── concordance.png          # Histogram of concordance, colored by isotype.
│   ├── concordance.pdf
│   ├── xconcordance.png         # Zoomed in view of the plots.
│   ├── xconcordance.pdf
│   │ 
│   │   # Results of checks Ye added. 
│   ├── npr1_allele_strain.tsv   # Any strains with the N2 npr-1 allele (lab dervied) is most likely not a wild isolate. 
│   ├── cutoff_distribution.tsv  # These files are results taking into consideration of introgression among strains.
│   ├── merge_betweengroup_pairwise_output.tsv
│   ├── new_isotype_groups.tsv
│   │ 
│   │   # Pairwise disconcordance plots.
│   └── pairwise
│       ├── within_group         # Plot for selective pairs of strains from the same isotypes. Note that not all pairs are plotted.
│       │                        # So far any strain pairs with red bars (discordance above cutoff) has been assigned to different isotypes. 
│       └── between_group        # Plot between strain pairs that meets Ye's criteria of introgression.
│                                # If process strain_pairwise_list_N2 is used, this step only compare each strain with N2.
│                                # If process strain_pairwise_list is used, this step compare every possible pair of strains.
└── varition
    └── out_gt.tsv               # this is created as input for scripts in "concordance-nf/manually_pairwise" folder.
```

# Relevant Docker Images

* `andersenlab/concordance` ([link](https://hub.docker.com/r/andersenlab/concordance)): Docker image is created within this pipeline using GitHub actions. Whenever a change is made to `env/concordance.Dockerfile` or `.github/workflows/build_docker.yml` GitHub actions will create a new docker image and push if successful
