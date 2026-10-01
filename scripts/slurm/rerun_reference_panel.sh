#!/bin/bash
# Re-run every step that depends on the published-signature panel, in order, after the custom signatures
# (config/subtype_signatures_custom.tsv) or the reference panel (config/reference_panel.yaml) change.
# From the project root: bash scripts/slurm/rerun_reference_panel.sh
#   07a signature library -> 07b NTP calls (incl. Montironi / Sia / WNT) + state/program mapping
#   -> 07c program/state annotation, prognosis beyond known signatures, head-to-head, proposed labels
#   -> 07d causal-flow labels and 08 external head-to-head (parallel) -> 09 publication figures
# Steps 04-06, 07e, 07f and 08b do not use the signatures and are not re-run.
set -euo pipefail
S=scripts/slurm
mkdir -p logs
a=$(sbatch --parsable -J panel_07a $S/07a_subtype_signatures.sbatch)
b=$(sbatch --parsable -J panel_07b --dependency=afterok:$a $S/07b_subtype_mapping.sbatch)
c=$(sbatch --parsable -J panel_07c --dependency=afterok:$b $S/07c_figures.sbatch)
d=$(sbatch --parsable -J panel_07d --dependency=afterok:$c $S/07d_causal_figures.sbatch)
e=$(sbatch --parsable -J panel_08 --dependency=afterok:$c $S/08_external_validation.sbatch)
f=$(sbatch --parsable -J panel_09 --dependency=afterok:$d:$e $S/09_publication_figures.sbatch)
echo "07a $a -> 07b $b -> 07c $c -> 07d $d + 08 $e -> 09 $f"
echo "After 07c: review results/07_post/figures/program_labels_proposed.tsv and update config/program_labels.tsv,"
echo "then re-run 09 alone (sbatch $S/09_publication_figures.sbatch) so the figures use the curated labels."
