# ==============================================================================
# chamber_constants.R -- closed-loop volumes shared by every goFlux auxfile
# ------------------------------------------------------------------------------
# Vtot = chamber + tubing + analyzer. The chamber and tubing terms differ by
# chamber design and stay in the scripts that build each auxfile; the analyzer
# term is the same instrument everywhere, so it is defined once, here.
#
# ANALYZER_VOLUME_CM3 is the internal volume of the ABB/LGR GLA131-GGA
# microportable ("LGR3", files UGGA_micro_*). The lab convention, adopted in
# ch4-data-filtering and shared across projects, is 28 cm3 (0.028 L).
#
# Until 2026-09-30 every pipeline here used 70 cm3. That is the value in goFlux's
# example auxfile, and it describes the larger benchtop LGR UGGA, not this
# instrument. Flux is proportional to Vtot, so the correction lowers every flux
# by 42 cm3 / Vtot: about 0.5-0.8 % for the 5-8 L semi-rigid and soil chambers
# and 1.8-6.8 % for the 0.5-2.3 L rigid chambers.
#
# Sourced by: 02_flux/semirigid/03_prep_{tree,soil}_auxfile.R, 04_goflux_soils.R,
#   02_flux/static/01_prep_auxfile{,_2023}.R, 02_flux/apply_auxfile_vtot.R,
#   03_merge/01_fix_soil_flux.R
# ==============================================================================
ANALYZER_VOLUME_CM3 <- 28
