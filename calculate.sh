#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

module load jq
# attempt to load different module versions in case of conflict
module load python/3.9 || module load python/3.7 || module load python/3.6
# remove the Picard header because it includes temporary paths
find . -xtype f -exec jq 'del(.picard | .header)' {} \; | python3 -mjson.tool --sort-keys
