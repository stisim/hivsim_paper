#!/bin/bash
# Launch seeds $1..$2 as agents, each detached (N agents: $3, default 100000)
IFS=$' \t\n'
cd "$(dirname "$0")"
for s in $(seq $1 $2); do
  setsid nohup /home/robyn/miniconda/envs/hp_rc3/bin/python -u run.py abm $s ${3:-100000} > logs/abm_${3:-100000}_$s.log 2>&1 < /dev/null &
done
