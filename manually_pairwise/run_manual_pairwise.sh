#!/usr/bin/env bash

# Extract genotypes of pairs in mannal_pairwise.txt  
echo "extract genotypes"
gt_file=/projects/b1059/projects/Ye/repo/concordance-nf/concordance-20190602/variation/out_gt.tsv
parallel --verbose csvtk cut -t -f CHROM,POS,{} ${gt_file} \> manual_pairwise_gt/{}.tsv < mannal_pairwise.txt

# Rename the outputs
echo "rename ..."
rename --subst , - manual_pairwise_gt/*.tsv

# run plot_pairwise.R to plot
echo "run plot_pairwise.R ..."
for i in `cat mannal_pairwise.txt`
do
x=`echo ${i} | cut -d "," -f1`
y=`echo ${i} | cut -d "," -f2`
Rscript plot_pairwise.R ${x} ${y} manual_pairwise_gt/${x}-${y}.tsv
done

# move pngs to manual_pairwise_plots folder
echo "move pngs to manual_pairwise_plots folder"
mv *.png manual_pairwise_plots
