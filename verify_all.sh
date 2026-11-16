#!/bin/sh
# Rerun every check listed in CHECKS.md.  Plain-Python steps use python3; steps marked (+) need numpy/igraph/pynauty/ortools,
# set PY to an interpreter that has them (default: python3).  Total runtime: about one to two hours.
set -e
PY=${PY:-python3}
cd "$(dirname "$0")"
for f in connsets/*.txt; do python3 build_from_logs.py "$f"; done           # every connection set, all vertex pairs
python3 paper/check_constructions.py                                        # Constructions 4.1-4.5 as printed
(cd family && python3 fibre_data.py > fibre_data.out && python3 family_theory.py > family_theory.out \
   && python3 shared_structure.py > shared_structure.out && python3 spread_inside.py > spread_inside.out)
(cd family && $PY family_structure.py > family_structure.out && $PY nauty_check.py > nauty_check.out \
   && $PY spread_link.py > spread_link.out && $PY oa_test.py > oa_test.out)                                   # (+)
(cd paper/checks && $PY clique_numbers.py 8 > clique_numbers.out && $PY spread_lattice.py spread_lattice.json > spread_lattice.out \
   && $PY known_constructions.py > known_constructions.out)                                                    # (+)
(cd family && for m in "5 0,1,4,5 0,1,2" "6 0,1,6 0,1" "7 0,4 0"; do set -- $m; $PY low_range/exhaustive_low.py $1 $2 $3 > low_range/exhaustive_m$1.out; done)   # (+)
echo "all checks ran; compare the .out files with the committed ones (git diff)"
