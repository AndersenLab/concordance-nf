There are some pairwises don't fit the cutoffs in this pipeline. But we need them to do futhur identification. As this process need manually pick the pairwise comparisons, so it's hard to incorporate it into the whole concordance-nf pipeline.

Once you got the pairwise comparisons, first, put them into `mannal_pairwise.txt`. Then update `out_gt.tsv` from the most recent runs. Finally, run `bash run_manual_pairwise.sh`. The plots will be in manual_pairwise_plots.

All the codes should be run under this directory. 

# exact genotypes of pairs in mannal_pairwise.txt 

```
parallel --verbose csvtk cut -t -f CHROM,POS,{} out_gt.tsv \> manual_pairwise_gt/{}.tsv < mannal_pairwise.txt
```

# rename the outputs

```
rename --subst , - manual_pairwise_gt/*.tsv
```


# run plot_pairwise.R to plot

```
for i in `cat mannal_pairwise.txt`
do
x=`echo ${i} | cut -d "," -f1`
y=`echo ${i} | cut -d "," -f2`
Rscript plot_pairwise.R ${x} ${y} manual_pairwise_gt/${x}-${y}.tsv
done
```

# move pngs to manual_pairwise_plots folder

```
mv *.png manual_pairwise_plots
```