#!/bin/bash
# The analysis scripts carry absolute paths from the cluster they were run on. This
# rewrites them to a new location in one pass.
#
#   ./retarget.sh /path/to/analysis_root
#
# The scripts expect this repository at <analysis_root>/revision_repo, and the input
# data files (see code/prepare_inputs.R) directly under <analysis_root>.
set -euo pipefail
OLD="/work/users/x/u/xuejun1/Integration_paper_Sim"
NEW="${1:?usage: ./retarget.sh /path/to/analysis_root}"
NEW="${NEW%/}"
n=$(grep -rl "$OLD" code 2>/dev/null | wc -l)
grep -rl "$OLD" code 2>/dev/null | xargs -r sed -i "s#$OLD#$NEW#g"
echo "rewrote $n file(s): $OLD -> $NEW"
left=$(grep -rl "$OLD" code 2>/dev/null | wc -l)
echo "remaining references to the old root: $left"
