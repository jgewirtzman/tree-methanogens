# Figure panel fixes (Referee 2 #6)

New scripts only; nothing edits the pipeline figure code. Each item below gives the
diagnosis, the fix, the drop-in output, and the caption/response text.

---

## Fig 2A — x-axis "non-intuitive and inconsistent"; outlier distortion; undefined dotted line
**R2:** "neither linear nor typical log scale... inconsistent across panels (e.g.
*Q. rubra* linear, *C. ovata* not)... distributions distorted by outliers (*A. rubrum*,
*Q. rubra*)... the dotted line should be explicitly defined."

**Diagnosis** (`code/01_flux_processing/static/04_height_effect_analysis.R:114-127, 205-241`):
each species panel was drawn on its **own free *linear* axis** whose three ticks were
literally that species' `(min, midpoint, max)`. So every panel had a different scale, a
single outlier maximum stretched the whole panel (squashing the rest near zero), and the
ticks looked neither standard-linear nor log. The dashed line is the **zero-net-flux**
reference (below = net uptake, above = net emission) — it was never defined.

**Fix** (`revision/analysis/R_fig2a_consistent_axis.R` → `fig2a_consistent_axis.png`):
one **shared signed pseudo-log (symlog) x-axis** for all panels — the same transform
family already adopted for the skewed regressions in this revision (`scales::pseudo_log_trans`,
σ = 0.005). Flux spans ~3 orders of magnitude (−0.028 to 6.2 nmol m⁻² s⁻¹) and includes
small negatives, so a shared *linear* axis would crush the low-flux species and a plain
log can't hold zero/negatives. Identical breaks on every panel (−0.01, 0, 0.01, 0.1, 1)
⇒ species directly comparable, outliers no longer distort, and the dashed zero line is
drawn and defined.

**Caption addition:** "CH₄ flux is shown on a shared signed pseudo-log (symlog) axis;
identical breaks on all panels allow direct comparison across species. The dashed line
marks zero net flux — points to its left indicate net uptake, to its right net emission."

---

## Fig 4 — "panel a and c missing legends for colors in stacked bar charts"
**R2** reviewed a version whose Fig 4 had **stacked bar charts in panels a and c**. That
figure was **restructured** in revision into the current two-panel form
(`code/02_ddpcr/util_combined_plot.R` + `04_species_barplots.R`): panel (a) = per-species
mcrA bars across four compartments; panel (b) = mcrA vs pmoA+mmoX scatter with marginal
densities. There are no stacked bars and no panel c any more.

**Both color mappings now carry titled legends** (verified in code and rendered figure):
- Panel (a): `scale_fill_viridis_c(name = "Phylogenetic Distance", …)`
  (`04_species_barplots.R:472`) — a viridis colorbar; endpoints are the nearest and
  farthest species (*Quercus*, *Tsuga*).
- Panel (b): "Sample Type" legend (Sapwood / Heartwood / Mineral / Organic).

So R2's "missing legends" is already resolved by the restructure. The one residual
ambiguity is that the panel-(a) colorbar shows only the endpoint genus names; the caption
should say what they mean:

**Caption addition (panel a):** "Bar fill encodes cophenetic phylogenetic distance among
species (viridis scale); the colorbar endpoints are labeled by the closest and most
distant taxa (*Quercus*, *Tsuga*)."

No new figure script is needed for Fig 4 (the fix is already in the current figure); this
is a caption clarification only.

---

## Fig 7c — missing CH₄ flux unit
**R2:** Fig 7c "missing the CH₄ flux unit."

**Diagnosis** (`code/06_figures/09_felled_oak_profiles.R:221`):
`ylab(expression(CH[4]~"Flux"))` — no unit. The flux is the goFlux `CH4_best.flux` (same
variable and pipeline as Fig 2; range here −0.005 to 0.40, median 0.039), i.e.
**nmol m⁻² s⁻¹**.

**Fix** (`revision/analysis/R_fig7c_flux_unit.R` → `fig7c_flux_unit_fixed.png`):
add the unit. **One-line patch** to apply in `09_felled_oak_profiles.R:221`:
```r
ylab(expression(CH[4]~"flux (nmol m"^-2~"s"^-1*")"))
```
The new script regenerates panel c standalone with the corrected label (and the same
internal-CH₄ point coloring) so the fix can be seen before merging.

---

## Response to Referee 2 (#6) — draft
> **Fig. 2A.** We have replaced the per-panel free axes with a single shared axis on all
> species panels. Because stem CH₄ flux spans roughly three orders of magnitude and
> includes small negative (uptake) values, we use a signed pseudo-log (symlog) scale with
> identical breaks on every panel, which makes species directly comparable and prevents
> individual outliers from distorting a panel. The dashed line, now defined in the
> caption, marks zero net flux (uptake to its left, emission to its right).
>
> **Fig. 4.** The figure has been restructured; the stacked bar panels have been replaced,
> and both colour mappings now carry explicit legends — a "Phylogenetic distance" colour
> bar for the per-species abundance panel (endpoints labelled by the closest and most
> distant taxa) and a "Sample type" legend for the abundance cross-plot.
>
> **Fig. 7c.** We have added the missing unit; the stem CH₄ flux axis now reads
> "CH₄ flux (nmol m⁻² s⁻¹)", consistent with Fig. 2.
