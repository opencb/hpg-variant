#!/bin/bash

PLINK=plink
PLINKSEQ=pseq

$PLINK --file $1 --make-bed --out temp
$PLINKSEQ $1-proj new-project
$PLINKSEQ $1-proj load-plink --file temp --id loaded
$PLINKSEQ $1-proj write-vcf > $2
# rm temp.*
