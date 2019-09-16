#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail
set -o noclobber

cd $1
source /.mounts/labs/PDE/public/bam-qc-metrics/miniconda3/bin/activate base && find . -type f -exec python -mjson.tool --sort-keys {} +
