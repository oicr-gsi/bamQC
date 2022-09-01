#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

module load jq
# load python only if it's not already loaded
module is-loaded python || module load python/3.10.4
# remove the Picard header because it includes temporary paths
find . -xtype f -exec jq 'del(.picard | .header)' {} \; | python3 -mjson.tool --sort-keys
