#!/bin/bash
cd /home/claude/work/drake
for a in "allsat 0" "allsat 1" "twosat 0" "twosat 1"; do set -- $a; python dk6_run.py $1 $2 > out/dk6/log_$1_$2.txt 2>&1 & done; wait
python dk6_wind.py > out/dk6/log_wind.txt 2>&1 &
for S in 0 10 25; do python dk6_detect.py $S 30 $S > out/dk6/log_detect_$S.txt 2>&1 & done; wait
echo ALLDONE
