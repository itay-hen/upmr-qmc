This is a new source snapshot for the 24 November 2026 manuscript revision. Original project sources were not replaced. Bins retain raw sums and closed counts. Estimates and delete-one-bin jackknifes combine those raw sums; --rawbins 1 prints the inputs needed to reconstruct every estimate. Worm move selection retains weight two.

## revision_checks/

Outputs behind manuscript v18: corrected Table h3 worm and "all" rows (raw-sum jackknife), the new seed-41 and beta=1 seed-45 runs, the C_worm replay, and the fixed-N connectivity checks. `I_finish_checks.py` recomputes the summary JSON from these files and the production outputs in `runs/`.
Note: `tables.tex` has not been regenerated; the corrected h3 rows are in `I_h3_worm_corrected.txt` and `I_h3_all_corrected.txt`.
