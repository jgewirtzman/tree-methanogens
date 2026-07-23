# Email to Wyatt Arnold — two ddPCR checks (final)

**Subject:** Two ddPCR methods checks for the tree-methane revision

Hi Wyatt,

Working through the New Phytologist revision and I've got two ddPCR methods details to confirm with you.

**1) Absolute copies/g — are we missing a ~10× template→reaction dilution factor?**

We may have missed a dilution correction step. The ddPCR raw files report concentrations per µL and per 20 µL reaction. I'd been assuming "per µL" meant per µL of eluent, and converting to copies per gram as:

copies/g ≈ (Conc, copies/µL) × 75 µL elution ÷ sample mass (g)

But now I'm thinking the QX200 reports copies/µL of the ~25 µL reaction well, which contained 2.5 µL of template in 22.5 µL reaction mix (i.e. a ~10× dilution in what the machine actually measures) — so the eluate concentration should be the well value × ~10 before we scale by the 75 µL elution volume. If that's right, our absolute mcrA/pmoA/mmoX copies/g are ~10× underreported. Just want to confirm that checks out — unless you'd set the machine to report per µL of eluent, or there's a correction already applied in the software?

**2) pmoA/mmoX primers + thermocycling — confirming for a reviewer who asked.**

The methods paper cites Bourne et al. 2001 and Fuse et al. 1998 but doesn't actually restate the primer sequences or reaction conditions — though you've got the details in the SI of your Microbiology Spectrum paper. Just want to confirm we used the identical chemistry for the trees as you'd used there — or I can't remember if there was any new thermocycler optimization. If there was, any chance you still have the notes/details? The Spectrum paper reports:

- pmoA: 189f (GGNGACTGGGACTTCTGG) / mb661r (CCGGMGCAACGTCYTTACC) — Bourne et al. 2001
- mmoX: 536f (CGCTGTGGAAGGGCATGAAGCG) / 898r (GCTCGACCTTGAACTTGGAGCC) — Fuse et al. 1998

and conditions:
- pmoA: 95 °C 5 min → 40× (95 °C 30 s / 52 °C 35 s / 72 °C 50 s) → 4 °C 5 min → 90 °C 5 min; ramp 1 °C/s
- mmoX: 95 °C 5 min → 40× (95 °C 30 s / 52 °C 30 s / 72 °C 45 s) → 4 °C 5 min → 90 °C 5 min; ramp 1.5 °C/s

Do those match what we ran on the tree samples? (mcrA is our probe assay, so that one's separate.) Also flagging: our current draft mistakenly cites Luesken 2011 / McDonald for pmoA/mmoX, which I'll swap to Bourne/Fuse unless you tell me otherwise.

Thanks — these are the last two things gating the methods.

Jon
