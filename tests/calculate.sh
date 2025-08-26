#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

module load jq
# extract the data for selected metrics
find . -xtype f -exec  jq '. | with_entries(select(.key | contains("mean","average","reads","total")))' {} \; 
