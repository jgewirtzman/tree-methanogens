# pmoA/mmoX split — robustness of the central gene–flux results (R2 #1c)

Script `R_pmoa_mmox_robustness.R` (reuses the exact per-species objects from
`02_scale_dependent_gene_patterns.R`; n = 10 species). Full table:
`revision/outputs/pmoa_mmox_robustness.txt`.

## Result

| Model | Methanotroph term | R² | p |
|---|---|---|---|
| Abundance → flux | mcrA | 0.394 | 0.052 |
| Abundance → flux | pmoA + mmoX (combined) | 0.230 | 0.161 |
| Abundance → flux | **pmoA only** | **0.254** | 0.138 |
| Abundance → flux | **mmoX only** | **0.039** | 0.584 |
| Ratio → flux | mcrA : (pmoA + mmoX) | 0.465 | 0.030 |
| Ratio → flux | **mcrA : pmoA** | **0.442** | **0.036** |
| Ratio → flux | **mcrA : mmoX** | **0.197** | 0.198 |

## Interpretation (for the R2 #1c response)

Reporting pmoA and mmoX separately does **not** change any major conclusion:

1. The combined pmoA+mmoX abundance→flux relationship is essentially the **pmoA**
   relationship (pmoA-only R² = 0.25 ≈ combined 0.23); mmoX alone is negligible
   (R² = 0.04). So the combined methanotroph term is carried by pmoA.
2. The central result — the methanogen:methanotroph **ratio predicts stem flux** —
   is **robust**: mcrA:(pmoA+mmoX) R² = 0.47 (p = 0.030) and mcrA:pmoA R² = 0.44
   (p = 0.036) are effectively identical and both significant. Using mmoX alone in
   the denominator weakens it (R² = 0.20, p = 0.198), consistent with pmoA being the
   numerically dominant, more informative methanotroph gene here.

**Takeaway:** the referee's concern (many methanotrophs encode both genes) is
addressed — the flux relationships are driven by pmoA, and the ratio conclusion
holds whether the methanotroph pool is defined as pmoA+mmoX or pmoA alone. This is
a strengthening result, reported in SI (table + `SI_fig_pmoa_mmox_separate.png`).

(These R²/p are all ratio- or log-based and thus invariant to the pending ×10 and to
the dry/wet-weight basis.)
