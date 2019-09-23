#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail
set -o noclobber

cd $1
module load python/3.6
find . -type f -exec python3 -mjson.tool --sort-keys {} +
