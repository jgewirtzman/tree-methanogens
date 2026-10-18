# ==============================================================================
# chamber_constants.R -- chamber, tubing and analyzer volumes shared by the
#                        goFlux auxfiles (Vtot = chamber + tubing + analyzer)
# ------------------------------------------------------------------------------
# ANALYZER_VOLUME_CM3: internal volume of the ABB/LGR GLA131-GGA microportable
# ("LGR3", files UGGA_micro_*), 28 cm3 (lab convention, shared across projects).
# (goFlux's example auxfile uses 70 cm3, the larger benchtop UGGA.)
#
# TUBING_CM3: 12 ft of 1/8" internal-diameter tubing, pi x (1/16 in)^2 x 144 in,
# in cm3 (x 16.387): about 29.0 cm3. Used by the semi-rigid stem and soil chambers.
#
# Soil collars: one geometry for every collar (25.43 cm interior diameter, 7.53 L
# headspace). Fits from saved windows use SOIL_VTOT_L, not the Vtot stored in the
# window file: the December 2020 windows were picked while the collar volume was
# the pre-correction 17.75 L (SOIL_CHAMBER_VOLUME_PRECORRECTION_L, cap 14.42 L + collar
# headspace 3.33 L, still written to the soil auxfile).
# ==============================================================================
ANALYZER_VOLUME_CM3 <- 28
TUBING_CM3 <- pi * (1/16)^2 * 12 * 12 * 16.387

SOIL_COLLAR_AREA_CM2 <- 507.7
SOIL_COLLAR_VOLUME_L <- 7.53
SOIL_CHAMBER_VOLUME_PRECORRECTION_L <- 17.75
SOIL_VTOT_L <- (SOIL_COLLAR_VOLUME_L * 1000 + TUBING_CM3 + ANALYZER_VOLUME_CM3) / 1000
