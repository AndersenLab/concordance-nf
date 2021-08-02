#!/usr/bin/env Rscript
library(gsheet)


args <- commandArgs(trailingOnly = TRUE)

# args[1] = species

if(args[1] == "c_elegans") {
    master <- gsheet::gsheet2tbl("https://docs.google.com/spreadsheets/d/10x-CcKNCl80F9hMcrGWC4fhP_cbekSzi5_IYBY2UqCc/edit#gid=538533765")
} else if(args[1] == "c_briggsae") {
    master <- gsheet::gsheet2tbl("https://docs.google.com/spreadsheets/d/1IJHMLwuaxS_sEO31TyK5NLxPX7_qSd0bHNKverAv8-0/edit")
} else if(args[1] == "c_tropicalis") {
    master <- gsheet::gsheet2tbl("https://docs.google.com/spreadsheets/d/1mqXOlUX7UeiPBe8jfAwFZnqlzhb7X-eKGK_TydT7Gx4/edit")
} else {
    print("ERROR: Please choose 'c_elegans', 'c_briggsae', or 'c_tropicalis' for '--params.species'")
}

readr::write_tsv(master, glue::glue("{args[1]}_species_sheet_{Sys.Date()}.tsv"))