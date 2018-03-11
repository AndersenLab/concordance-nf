#!/bin/bash
# Run this on Quest

cat SM_sample_sheet.tsv | gawk -v OFS=$'\t' '
  function basename(file) {
    sub(".*/", "", file)
    return file
  }
NR % 100 == 1 { print $1, $2, $3, basename($4), basename($5), $6 }' > test_data/SM_sample_sheet.tsv

parallel --dryrun --verbose "zcat {} | head -n 100000 | gzip > test_data/{/}" ::: `cat SM_sample_sheet.tsv | awk 'NR % 100 == 1 { print $4; print $5 }'`