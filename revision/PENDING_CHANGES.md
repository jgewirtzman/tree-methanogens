# Pending changes still to apply (as of this session)

## A. DATA corrections — gate multiple figures + text numbers

### A1. ×10 ddPCR dilution correction — WAITING ON WYATT
- Template→reaction dilution (~2.5 µL in ~25 µL) may have been dropped, so absolute
  copies/g may be ~10× underreported.
- **Affects (absolute copies/g only):** Fig 4 (a,b), Fig 7b, pmoA/mmoX SI figure/tables,
  and every copies/g number in Methods/Results/captions.
- **NOT affected:** anything relative/ratio/log — Fig 5 (%), Fig 8 (ratios & correlations),
  pmoA:mmoX ratios, gene–flux R².
- **Status:** `DILUTION_10X` toggle already wired in `R_fig7_final.R` and
  `R_pmoa_mmox_separate.R` (set 1→10 when confirmed). Fig 4 would need the same toggle.

### A2. Dry vs wet weight harmonization (R2 #1 units) — ON HOLD (Jon)
- Wood copies/g is DRY-weight (cores freeze-dried); soil copies/g is WET/fresh-weight.
  Referee: use the same basis.
- **Fix:** apply moisture (VWC/gravimetric) correction to put soil on a dry-weight basis
  (or state both + a labeled conversion), throughout figures + text.
- **Affects:** Fig 4 (soil bars/scatter), pmoA/mmoX SI (soil points), Methods/Results
  copies/g numbers. Wood already dry. **Relative/among-wood and ratio results unaffected.**
- **Status:** HELD OFF per Jon. Needs the soil moisture data + a decided convention.

## B. FIGURES

### B1. Fig 4 (methanogen/methanotroph abundance) — pending A1/A2 only
- Otherwise fine (legend fix was caption-only). When A1/A2 resolve, rebuild as a new
  `R_fig4_final.R` with the ×10 toggle and the harmonized basis.

### B2. Fig 6 restructure — DECIDED, not executed
- Hydrogenotrophy synthesis becomes **Fig 6**; old PICRUSt *mcrA* heatmap moves to **SI**.
- Also (STATUS): promote **Methods S1–S3 to main**.
- Then **renumber** the whole main set and update all in-text "Fig N" references.

### B3. Fig 8 disks → SI — OPTIONAL (Jon's call)
- Move the radial "disk" panels (a–c) to SI, keep the scatter/R² panels (d–g) in main,
  if a leaner molecular section is wanted.

### B4. Fig 6 loose-vs-strict mcrA reconciliation — OPEN (Jon)
- Discussion cites strict families; Fig 6/S-heatmaps used the loose gene. Decide/annotate.

## C. TEXT / RESPONSE ASSEMBLY (drafts exist; fold in at text stage)
- Primer citation fix L117 (Luesken/McDonald → **Bourne 2001 + Fuse 1998**) — pending
  Wyatt confirm of primers/thermocycling (`primer_sequences.md`).
- Gene-name italics (*mcrA*, *pmoA*, *mmoX*) throughout.
- Define ddPCR at first use; per-gram-basis clause in each copies/g caption (R2 #1 units).
- Ranniku 2023 cite; Gauci discussion; mechanistic sentence fixes (`microbiology_corrections.md`).
- Plant-trait SI figure(s) + discussion paragraph (`plant_traits_SI_and_discussion.md`) —
  pick composite vs separate SI figs.
- FAPROTAX/PICRUSt caveats (R2 #3), L342 exponential defense — drafted.
- New main-figure captions (budget, hydrogenotrophy) — `new_figure_captions.md`.
- Assemble response-to-referees (`response_to_referees_skeleton.md`).

## Done this session (figures)
Fig 2 final (numeric axis, mean±SE + linear trends), Fig 5 (revised putative),
Fig 7 final (copies/g + flux unit + 10× toggle), Fig 9 budget (graticule, palette,
per-bar values), hydrogenotrophy synthesis (Fig 6 content), pmoA/mmoX SI (copies/g +
separate-gene robustness: conclusions hold, carried by pmoA).
