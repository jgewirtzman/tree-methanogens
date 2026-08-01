
# Random Forest (RF) Workflow for Plot‑Scale, Monthly CH₄ Flux (Prediction‑First)

> **Your constraints honored**
> - **No hardcoded filenames** → all sources/outputs are referred to as **PLACEHOLDER_…** (you will supply actual paths later).
> - **No sine/cosine month encoding** → use an **empirical smoothed seasonal index** (LOESS/spline over monthly residuals).
> - **Keep it simple** → no shading/microsite corrections; focus on robust, predictive scaling.
> - **Harmonize chambers** → correct monthly-chamber underestimation vs July chamber **before** training; do **not** include `chamber_type` in final RF.

---

## 0) Conventions & Units
- Model fluxes in **μmol m⁻² s⁻¹** (soil and stem). Convert to **mg CH₄ m⁻² d⁻¹** only at the end.
- Response transform: **asinh(y)** for training; back‑transform predictions with **sinh(ŷ)** (handles sign changes).
- Coordinates are **meters** in a consistent projected CRS for spatial blocking and raster sampling.

---

## 1) Inputs (placeholders; you will provide the actual sources)
- **PLACEHOLDER_TREE_JULY** — July stem‑flux campaign (~150 trees, ~15 spp).  
  Columns: `tree_id, species, dbh_m, datetime, month, stem_flux_umol_m2_s, air_temp_C, soil_temp_C, soil_moisture_abs, chamber_type, x, y`.
- **PLACEHOLDER_TREE_YEAR** — Year‑long monthly stem‑flux subset (~45 trees × ~12).  
  Same columns; **must include `chamber_type`**.
- **PLACEHOLDER_SOIL_YEAR** — Year‑long monthly soil chambers (36 sites).  
  Columns: `site_id, datetime, month, soil_flux_umol_m2_s, air_temp_C, soil_temp_C, soil_moisture_abs, x, y`.
- **PLACEHOLDER_INVENTORY** — Plot inventory (~7000 trees).  
  Columns: `tree_id, species, dbh_m, x, y`.
- **PLACEHOLDER_MOISTURE_DEC_RASTER** — December soil‑moisture raster (absolute) or table `{x, y, moisture_dec_abs}`.
- **PLACEHOLDER_DRIVERS** — Monthly plot‑level temps. Columns: `{month, air_temp_C_mean, soil_temp_C_mean}`.
- **PLACEHOLDER_TAXONOMY** — Taxon map: `species → genus → family → order → class`.
- **PLACEHOLDER_PLOT_AREA** — Scalar area of plot (m²).

---

## 2) Monthly Moisture Map (absolute) via Affine Calibration
**Goal:** Build an absolute moisture field for each month *t* that preserves December spatial ranks.

1. Read **December grid** `M_dec(x,y)` from **PLACEHOLDER_MOISTURE_DEC_RASTER**.
2. For each month *t=1..12*, gather point moisture observations (trees + soils measured in month *t*) with coordinates `(x_p, y_p)` and fit:
   \[ \theta_{p,t} = \alpha_t + \beta_t \cdot M_{dec}(x_p, y_p) + \varepsilon \]
   Use OLS or robust regression. Store `(alpha_t, beta_t)` and `R²_t`.
3. Define monthly absolute moisture surface:
   \[ \widehat{M}_t(x,y) = \alpha_t + \beta_t \cdot M_{dec}(x,y) \]
4. Clip \(\widehat{M}_t\) to plausible site bounds (e.g., 0–0.6 m³/m³ or site‑specific).

**Output placeholder:** `PLACEHOLDER_MOISTURE_AFFINE_TABLE` with `month, alpha_t, beta_t, R2_t`.

---

## 2.5) Chamber Calibration (Tree Flux Harmonization)
**Problem:** July and monthly stem chambers differ; monthly chamber likely **underestimates** flux.

**Approach (on asinh scale):**
1. Stack July + Year tree measurements; define:
   - Response: `z = asinh(stem_flux_umol_m2_s)`
   - Predictors (controls): `air_temp_C`, (`soil_temp_C` optional), `soil_moisture_abs`, species one‑hot
   - Indicator: `I_monthly = 1` for monthly chamber, `0` for July chamber (reference)
2. Fit OLS:
   \[ z = \beta_0 + \beta_1 I_{monthly} + \beta_2 airT + \beta_3 soilT + \beta_4 moisture + \sum_s \gamma_s species_s + \varepsilon \]
   `β1` estimates the systematic bias of the monthly chamber relative to July.
3. Correct **only** monthly‑chamber rows:
   - `z_corr = z − β1` if `I_monthly=1`; July rows unchanged.
   - Train the **final TreeRF on `z_corr`** (asinh) as target. **Do not include `chamber_type` in final RF.**
4. Save `β1`, SE, 95% CI, R². Check post‑correction residuals vs env. variables (should not be chamber‑structured).

**Bootstrap:** Re‑estimate `β1` in every replicate (see §8).

---

## 3) Empirical Smoothed Seasonal Index (no sine/cosine)
Create **SI_tree[month]** and **SI_soil[month]** from data.

**Steps (repeat separately for tree + soil):**
A. Fit a provisional RF **without any month feature** (use sections §4.1/§4.2 features minus SI).  
B. Obtain **out‑of‑fold predictions**; compute residuals `r = y_asinh − yhat_asinh`.  
C. Aggregate to monthly means: `r̄[m] = mean(r | month=m)`.  
D. Smooth `{(m, r̄[m]) : m=1..12}` with **LOESS or cubic smoothing spline**; enforce cyclic continuity (e.g., mirror months).  
E. The smoothed curve is **SI_tree[m]** / **SI_soil[m]** (asinh‑scale). Join it back as a single numeric feature.

---

## 4) Features & Targets

### 4.1 Tree model features (per observation)
- Taxonomy one‑hots: **species**, plus **genus** and **family** (collapse rare to `GENUS_OTHER`, `FAMILY_OTHER`).
- `dbh_m` (continuous).
- `air_temp_C_mean` (from **PLACEHOLDER_DRIVERS** by month). *(Optional `soil_temp_C_mean`.)*
- `soil_moisture_abs_at_tree`: observed in training; for prediction use \(\widehat{M}_t(x_i,y_i)\).
- `SI_tree[month]` (from §3, asinh‑scale residual index).
- **Engineered interactions (biologically motivated):**
  - `moisture_x_airT = soil_moisture_abs_at_tree * air_temp_C_mean`
  - *(optional)* `moisture_x_soilT = soil_moisture_abs_at_tree * soil_temp_C_mean`
  - `moisture_x_speciesOneHot_*` (moisture × each species one‑hot; ~15 columns)
  - *(optional)* `moisture_x_genusOneHot_*`
- **Taxon prior (numeric):** `taxon_prior_asinh` — median of **asinh‑scale residuals** at the **lowest available** taxonomic level (species→genus→family→order→class; see §4.3).

**Target:** `y_tree = asinh(stem_flux_umol_m2_s)` (after chamber correction where applicable).

### 4.2 Soil model features (per observation)
- `soil_temp_C_mean` (from **PLACEHOLDER_DRIVERS**). *(Optional `air_temp_C_mean`.)*
- `soil_moisture_abs_at_site`: observed in training; for prediction use \(\widehat{M}_t(x,y)\).
- `SI_soil[month]`.
- **Engineered interactions:**
  - `moisture_x_soilT = soil_moisture_abs_at_site * soil_temp_C_mean`
  - *(optional)* `moisture_x_airT = soil_moisture_abs_at_site * air_temp_C_mean`

**Target:** `y_soil = asinh(soil_flux_umol_m2_s)`.

**Extrapolation guard (both models):** Clip each numeric driver to training **[p1, p99]** before prediction.

### 4.3 Taxonomy‑aware handling (unmeasured species; medians, not means)
**Input:** **PLACEHOLDER_TAXONOMY**.

**Purpose:** For inventory species unseen in training, avoid a generic UNKNOWN by leveraging higher‑level taxa.

**Construction of `taxon_prior_asinh`:**
1) Fit a provisional **environment‑only** RF (no taxon one‑hots) on tree data; get out‑of‑fold `yhat_asinh`.  
2) Compute residuals `r = y_asinh − yhat_asinh`.  
3) Group by taxon level and compute **medians** (require Nmin≥5):  
   `med_species[s]`, `med_genus[g]`, `med_family[f]`, `med_order[o]`, `med_class[c]`.  
4) For each row (training + prediction), set `taxon_prior_asinh` to the **lowest available** median in priority order: species→genus→family→order→class→0.  
5) Train final TreeRF with **species/genus/family one‑hots + `taxon_prior_asinh`**. For unseen species, genus/family one‑hots and the prior still inform predictions.

---

## 5) Random Forest Models (simple, robust)
Two regressors: **TreeRF** and **SoilRF** (scikit‑learn `RandomForestRegressor` or equivalent).

**Baseline hyperparameters:**
- `n_estimators = 800`
- `max_depth = None`
- `min_samples_leaf = 5`
- `max_features = "sqrt"` (or `0.6`)
- `bootstrap = True`
- `oob_score = True`
- `random_state = 42`

**Training data:**
- **TreeRF:** stack **PLACEHOLDER_TREE_JULY** + **PLACEHOLDER_TREE_YEAR**, after **chamber correction** (§2.5). Target = `y_tree`.
- **SoilRF:** **PLACEHOLDER_SOIL_YEAR**. Target = `y_soil`.

Fit on **asinh** targets; back‑transform predictions with **sinh()** for reporting.

---

## 6) Cross‑Validation (minimal)
- Report **OOB R²** and RMSE on **back‑transformed units** (μmol m⁻² s⁻¹).
- Optionally **temporal block CV** (leave out months) to confirm the SI term captures residual seasonality.
- **Species hold‑out**: leave out a July‑only species to assess generalization.

---

## 7) Monthly Prediction & Area Scaling

### 7.1 Geometry
- Stem lateral area to **2 m** height:
  \[ S_i = \pi \cdot dbh_i \cdot 2 \]  (dbh in meters)
- Basal area:
  \[ BA_i = \pi \left(\frac{dbh_i}{2}\right)^2 \]
- Soil ground area:
  \[ A_{soil} = A_{plot} - \sum_i BA_i \]

### 7.2 Tree predictions (per month *t*)
For each inventory tree *i* (from **PLACEHOLDER_INVENTORY**):
1) Build features: taxon one‑hots, `dbh_m`, `air_temp_C_mean[t]` (and `soil_temp_C_mean[t]` if used), `SI_tree[t]`, and moisture \(\widehat{M}_t(x_i,y_i)\).  
2) Predict `ŷ_asinh = TreeRF(X_i,t)` → back‑transform `\(\widehat{F}^{tree}_{i,t} = sinh(ŷ_asinh)\)` (μmol m⁻² s⁻¹ per stem area).  
3) Tree‑level rate (μmol s⁻¹): \( R^{tree}_{i,t} = \widehat{F}^{tree}_{i,t} \cdot S_i \).  
4) Plot‑basis tree flux:
   \[ \Phi^{tree}_t = \frac{\sum_i R^{tree}_{i,t}}{A_{plot}} \]

### 7.3 Soil predictions (per month *t*)
**Option A — Raster/pixel integration (preferred):**
- For each pixel *p* (area `a_p`, center `(x_p,y_p)`): features = `soil_temp_C_mean[t]`, `SI_soil[t]`, `\widehat{M}_t(x_p,y_p)` (and optional air temp).  
- Predict \(\widehat{F}^{soil}_{p,t}\) (μmol m⁻² s⁻¹).  
- Area‑weighted mean: \(\bar{F}^{soil}_t = \sum_p \widehat{F}^{soil}_{p,t} a_p / \sum_p a_p\).  
- Scale by soil fraction: \(\Phi^{soil}_t = \bar{F}^{soil}_t \cdot (A_{soil} / A_{plot})\).

**Option B — Point grid:** predict on a regular grid; take mean and apply soil fraction.

### 7.4 Total monthly plot flux
\[\Phi^{plot}_t = \Phi^{tree}_t + \Phi^{soil}_t \]  (μmol m⁻² s⁻¹)

### 7.5 Unit conversion
To **mg CH₄ m⁻² d⁻¹**:  
\[\text{mg m}^{-2}\text{d}^{-1} = \Phi \; (\mu mol\, m^{-2}\, s^{-1}) \times 86400 \times 16 \times 10^{-3}\]

---

## 8) Uncertainty (spatially blocked bootstrap)
Use **spatial block resampling** to respect autocorrelation; fall back to simple bootstrap if coords missing.

**A) TreeRF bootstrap**
1) **Spatial blocks:** cluster tree `(x,y)` into `B_tree` blocks (k‑means with k≈8–15 or fixed grid).  
2) **Resample:** sample blocks with replacement; include all rows from selected blocks.  
3) **Per replicate:**  
   - Refit chamber bias `β1` (§2.5) and correct monthly‑chamber rows.  
   - Recompute `SI_tree` (§3).  
   - Refit TreeRF.

**B) SoilRF bootstrap**
1) Block soil sites `(x,y)` similarly; resample blocks.  
2) Recompute `SI_soil`; refit SoilRF.

**C) Moisture calibration uncertainty**
- For each month, bootstrap the calibration pairs `(θ_{p,t}, M_dec(x_p,y_p))` by spatial blocks if possible; refit `(alpha_t, beta_t)`.

**D) Prediction aggregation**
- In each replicate, compute `Φ_tree_t, Φ_soil_t, Φ_plot_t`.  
- Summarize **median & 95% CI** per month and annual.

**Pseudocode (sketch):**
```python
def spatial_bootstrap(df, coords=('x','y'), n_blocks=10, B=200):
    # 1) assign block_id via k-means or grid
    # 2) for each b in 1..B:
    #       choose blocks with replacement
    #       df_b = rows where block_id in chosen
    #       refit chamber bias, SI, RFs, moisture affine
    #       predict & scale → Φ_tree, Φ_soil, Φ_plot
    # 3) return monthly medians & CIs
```

---

## 9) Minimal Function Spec (guide for implementation)
- `fit_moisture_affine(month, points_df, M_dec_lookup) -> (alpha, beta, R2)`
- `make_Mhat(alpha_beta_tbl, M_dec_lookup) -> f(month, x, y) -> M_hat`
- `train_baseline_RF(df, features_wo_month) -> RF0`  *(for SI residuals)*
- `compute_empirical_SI(df, RF0) -> SI_table[{group='tree'|'soil'}][month, SI]`
- `build_features_tree(df, drivers, M_hat, SI_table, taxonomy, options) -> X_tree, y_asinh`
- `build_features_soil(df, drivers, M_hat, SI_table, options) -> X_soil, y_asinh`
- `compute_taxon_prior(df, taxonomy, env_only_oof) -> taxon_prior_asinh`
- `train_rf(X, y_asinh, params) -> RF`
- `predict_trees_month(inventory, month, RF, drivers, M_hat, SI_table, taxonomy) -> per_tree_flux`
- `predict_soil_month(grid_or_raster, month, RF, drivers, M_hat, SI_table) -> per_pixel_flux`
- `scale_tree_month(per_tree_flux, inventory, A_plot) -> Phi_tree_t`
- `scale_soil_month(per_pixel_flux, A_soil, A_plot) -> Phi_soil_t`
- `convert_units(muMol_m2_s) -> mg_m2_d`
- `spatial_bootstrap_pipeline(B, ...) -> monthly_summary, annual_summary`

---

## 10) Guardrails & Logging
- **Chamber correction:** predictions are on the harmonized scale via §2.5; **no** chamber indicator used at prediction time.
- **Clipping**: clip numeric drivers to training [p1, p99] before prediction.
- **Species coverage**: taxonomy encoding (species/genus/family) + `taxon_prior_asinh`; log count of inventory trees using genus/family fallback.
- **DBH sanity**: drop/cap impossible values; log affected rows and % basal area affected.
- **Negative flux** allowed; do not floor at zero.

---

## 11) Quality Control Metrics & Plots (outputs)
Produce these diagnostics (output names are placeholders):
1) **Chamber correction stability** — plot `β1` (±CI) vs `air_temp_C`, `soil_temp_C`, `soil_moisture_abs` (should be roughly constant).  
   Output: `PLACEHOLDER_QC_CHAMBER_STABILITY`.
2) **Prediction vs observation, by chamber type** (post‑correction).  
   Output: `PLACEHOLDER_QC_PRED_VS_OBS_BY_CHAMBER`.
3) **Species‑specific OOB R² (TreeRF)**; highlight July‑only species.  
   Output: `PLACEHOLDER_QC_SPECIES_R2`.
4) **Moisture calibration R² by month** with `(alpha_t, beta_t)`.  
   Output: `PLACEHOLDER_QC_MOISTURE_CALIB`.
5) **Partial dependence / ALE** for moisture and temperature in both models.  
   Output: `PLACEHOLDER_QC_PDP_ALE`.
6) *(Optional)* **Residual maps** (spatial patterns).  
   Output: `PLACEHOLDER_QC_RESIDUAL_MAPS`.

---

## 12) Outputs to Produce (placeholders)
- `PLACEHOLDER_MONTHLY_FLUXES` — table with `month, Phi_tree_umol_m2_s, Phi_soil_umol_m2_s, Phi_plot_umol_m2_s` and converted `mg_CH4_m2_d`.
- `PLACEHOLDER_DIAGNOSTICS` — OOB R², RMSE, feature importances for TreeRF/SoilRF.
- `PLACEHOLDER_MOISTURE_AFFINE_TABLE` — `(alpha_t, beta_t, R2_t)` per month.
- `PLACEHOLDER_SI_TABLES` — `SI_tree[month]`, `SI_soil[month]`.
- `PLACEHOLDER_TAXONOMY_PRIORS` — median residuals per species/genus/family/order/class.
- `PLACEHOLDER_COVERAGE_REPORT` — counts, taxonomy fallback %, soil area fraction.
- QC plots listed in §11.

---

### Notes on Simplicity and Extensibility
- You can swap RFs for **Quantile Regression Forests** with no design changes if you want direct prediction intervals.
- If raster processing is heavy, use a regular grid of points; document the grid density.

---

**End of Spec — Ready to Implement**
