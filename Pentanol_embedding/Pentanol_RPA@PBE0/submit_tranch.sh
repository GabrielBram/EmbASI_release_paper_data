#!/bin/bash

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <nnodes> <nnodesperjob> <numsubjobs> <basis> <nncut>"
    exit 1
fi

let nnodes=$(( $1-1 ))
let nnodespersubjob=$2
let subjobperjob=$(( $1 / $2 ))

hl="PBE0"
ll="PBE0"
toten_method="RPA"
basis=$4
let nncut=$5

for ((tranche_start=0 ; tranche_start<$3 ; tranche_start+=$subjobperjob)); do
    let tranche_end=$(( $tranche_start + $subjobperjob - 1 ))
    sbatch submission_dimer_tranche_base.script $nnodespersubjob $tranche_start $tranche_end $hl $ll $toten_method $basis $nncut
done

sbatch submission_monomer_base.script $nnodespersubjob $hl $ll $toten_method $basis $nncut
