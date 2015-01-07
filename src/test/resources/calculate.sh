#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail
set -o noclobber

cd $1
find . -type f -exec python -mjson.tool {} +
find . -type f -exec jsonToGenericRunReport.pl {} + | grep -v "Generic run report generated on"
