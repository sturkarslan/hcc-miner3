#!/bin/bash
# Submit the full leave-one-cohort-out chain from the project root: bash scripts/slurm/08b_loco_submit.sh
# per held-out cohort: prepare -> MINER mechinf -> technical-regulon filter -> MINER subtypes_filtered -> causal;
# then one comparison job after all three.
set -euo pipefail
S=scripts/slurm
last=()
for H in TCGA CLCA LICA_FR; do
  M=loco_no$H
  a=$(sbatch --parsable -J loco_prep_$H $S/08b_loco_hccprep.sbatch prepare --held $H)
  b=$(sbatch --parsable -J loco_mech_$H --dependency=afterok:$a $S/04_run_miner.sbatch --matrix $M --steps mechinf)
  c=$(sbatch --parsable -J loco_filt_$H --dependency=afterok:$b $S/08b_loco_hccprep.sbatch filter --held $H)
  d=$(sbatch --parsable -J loco_subt_$H --dependency=afterok:$c $S/04_run_miner.sbatch --matrix $M --steps subtypes_filtered)
  e=$(sbatch --parsable -J loco_caus_$H --dependency=afterok:$d $S/05_causal_inference.sbatch --matrix $M)
  echo "$H: prepare $a -> mechinf $b -> filter $c -> subtypes $d -> causal $e"
  last+=($e)
done
f=$(sbatch --parsable -J loco_compare --dependency=afterok:$(IFS=:; echo "${last[*]}") $S/08b_loco_compare.sbatch)
echo "compare: $f"
