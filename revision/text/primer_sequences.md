# pmoA / mmoX / mcrA primer sequences for Methods (R2 #1e / L248-249)

Source: Jon's SI Table 2 (from the methods paper, Arnold et al. 2024, *Methods Ecol. Evol.*) +
the methods-paper text. **ALL TO BE CONFIRMED WITH WYATT** (which primers + which thermocycling).

## The assays used in THIS manuscript
| Gene | Primers (5'→3') | Chemistry | Source |
|---|---|---|---|
| **mcrA** | F `ACGACYTRCAGGAYCAGTGY` · Probe `WGGWCCWAACTAYCCBAACTACG` · R `TGGTGWCCBACGTTCATYG` (118 bp) | probe (FAM) | **novel, this study** (Arnold et al. 2024) — replaces the EvaGreen ML assay |
| **pmoA** | `189f` GGNGACTGGGACTTCTGG · `mb661r` CCGGMGCAACGTCYTTACC | EvaGreen | **Bourne et al. 2001** |
| **mmoX** | `536f` CGCTGTGGAAGGGCATGAAGCG · `898r` GCTCGACCTTGAACTTGGAGCC | EvaGreen | **Fuse et al. 1998** |

(The mcrA *ML* primers — MLf `GGTGGTGTMGGATTCACACARTAYGCWACAGC` / MLr `TTCATTGCRTAGTTWGGRTAGTT`,
Luton et al. 2002 — were the initial EvaGreen mcrA approach, superseded here by the probe assay; cite
only if the text describes the switch.)

## Thermocycling (from SI Table 2 — CONFIRM w/ Wyatt; he notes it may follow Arnold et al. 2023, *Microbiol. Spectr.* 11(5), spectrum.02714-23)
- **pmoA:** 95 °C 5 min; ×40 [95 °C 30 s, 52 °C 35 s, 72 °C 50 s]; 4 °C 5 min; 90 °C 5 min. Ramp 1 °C/s. Std curve R²=0.98.
- **mmoX:** 95 °C 5 min; ×40 [95 °C 30 s, 52 °C 30 s, 72 °C 45 s]; 4 °C 5 min; 90 °C 5 min. Ramp 1.5 °C/s. R²=0.99.
- **mcrA (probe):** 95 °C 10 min; ×40 [94 °C 30 s, 48 °C 1 min 20 s]; 98 °C 10 min (already in L117).

## ⚠️ Citation correction needed (L117)
The manuscript currently cites **Luesken et al. 2011; McDonald & Murrell 1997; McDonald et al. 1995**
for pmoA/mmoX. Per the methods paper these should be **Bourne et al. 2001 (pmoA)** and **Fuse et al.
1998 (mmoX)**. Replace the citation and add the sequences. (Confirm w/ Wyatt before finalizing.)

## Drop-in Methods text (replaces the pmoA/mmoX sentence at L117)
> For pmoA and mmoX we used previously published EvaGreen ddPCR assays: pmoA primers 189f
> (5′-GGNGACTGGGACTTCTGG-3′) and mb661r (5′-CCGGMGCAACGTCYTTACC-3′) (Bourne et al. 2001), and mmoX
> primers 536f (5′-CGCTGTGGAAGGGCATGAAGCG-3′) and 898r (5′-GCTCGACCTTGAACTTGGAGCC-3′) (Fuse et al.
> 1998). Primer sequences and thermocycling conditions for all assays are given in Table S2.

## Suggested SI table (drop-in; mirrors SI Table 2)
| Target | Primers (5'→3') | PCR schedule | Ramp | Std curve |
|---|---|---|---|---|
| mcrA (probe) | F ACGACYTRCAGGAYCAGTGY / Probe WGGWCCWAACTAYCCBAACTACG / R TGGTGWCCBACGTTCATYG | 95 °C 10 min; ×40 [94 °C 30 s, 48 °C 80 s]; 98 °C 10 min | — | — |
| pmoA | 189f GGNGACTGGGACTTCTGG / mb661r CCGGMGCAACGTCYTTACC | 95 °C 5 min; ×40 [95 °C 30 s, 52 °C 35 s, 72 °C 50 s]; 4 °C 5 min; 90 °C 5 min | 1 °C/s | R²=0.98 |
| mmoX | 536f CGCTGTGGAAGGGCATGAAGCG / 898r GCTCGACCTTGAACTTGGAGCC | 95 °C 5 min; ×40 [95 °C 30 s, 52 °C 30 s, 72 °C 45 s]; 4 °C 5 min; 90 °C 5 min | 1.5 °C/s | R²=0.99 |

## Response to Referee 2 (#1e)
> The primer sequences and thermocycling conditions for all ddPCR assays (mcrA, pmoA, mmoX) are now
> provided in full in Table S2 and in the Methods, so readers need not consult the source references.
> pmoA was amplified with primers 189f/mb661r (Bourne et al. 2001) and mmoX with 536f/898r (Fuse et
> al. 1998); mcrA used our probe-based assay. We have corrected the primer citations accordingly.
