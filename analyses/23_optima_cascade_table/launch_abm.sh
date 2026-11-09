#!/bin/bash
# Seeds FIRST..LAST at N agents, detached: bash launch_abm.sh 1 16 100000
cd "$(dirname "$0")"
for s in $(seq $1 $2); do
  setsid nohup conda run -n hp_rc3 python run.py abm $s ${3:-100000} > logs/abm_${3:-100000}_$s.log 2>&1 < /dev/null & disown
done
