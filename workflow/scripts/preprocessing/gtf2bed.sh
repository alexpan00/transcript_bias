#!/bin/bash
set -euo pipefail
# $1: gtf file
# $2: bed file
# $3: tama gtf2bed script
gffread -T $1 -o ${2}.gtf
python $3 ${2}.gtf ${2} 
rm ${2}.gtf 