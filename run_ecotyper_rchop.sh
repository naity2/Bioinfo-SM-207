#! /usr/bin/bash

cd /home/ubuntu/app/ecotyper
Rscript EcoTyper_recovery_bulk.R \
    -d Lymphoma \
    -m /home/ubuntu/data/bioinfo_sm207/schmitz_tpm.txt \
    -t 8 \
    -o /home/ubuntu/data/bioinfo_sm207/ecotyper_outs

Rscript EcoTyper_recovery_bulk.R \
    -d Lymphoma \
    -m /home/ubuntu/data/bioinfo_sm207/reddy_tpm.txt \
    -t 8 \
    -o /home/ubuntu/data/bioinfo_sm207/ecotyper_outs