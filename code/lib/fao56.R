# fao56.R -- FAO-56 daily reference evapotranspiration (one definition)
# ------------------------------------------------------------------------------
# et0_fao56(date, Tmax, Tmin, Tmean, Rs, ea): Penman-Monteith reference ET (mm day-1)
# for a grass reference at the plot, and the Hargreaves-Samani value as a check.
#   Rs  incoming shortwave, MJ m-2 day-1;  ea  actual vapour pressure, kPa.
# Wind is not available at the plot (the tower's wind column is unusable): u2 = 2 m s-1,
# the FAO-56 default for missing wind. Used by code/04_drivers/01_wb_reference_et.R
# (tower record).
FAO56_LAT_DEG <- 41.9899      # plot centre
FAO56_ELEV_M  <- 215          # Yale-Myers, approximate
FAO56_ALBEDO  <- 0.23         # FAO-56 reference grass
FAO56_SIGMA   <- 4.903e-9     # Stefan-Boltzmann, MJ K-4 m-2 day-1
FAO56_U2      <- 2.0          # m s-1, FAO-56 default for missing wind

et0_fao56 <- function(date, Tmax, Tmin, Tmean, Rs, ea) {
  doy <- as.integer(format(as.Date(date), "%j"))
  phi <- FAO56_LAT_DEG * pi/180
  dr  <- 1 + 0.033*cos(2*pi*doy/365)
  dec <- 0.409*sin(2*pi*doy/365 - 1.39)
  ws  <- acos(pmin(pmax(-tan(phi)*tan(dec), -1), 1))
  Ra  <- (24*60/pi) * 0.0820 * dr * (ws*sin(phi)*sin(dec) + cos(phi)*cos(dec)*sin(ws))
  Rso <- (0.75 + 2e-5*FAO56_ELEV_M) * Ra
  esat <- function(t) 0.6108*exp(17.27*t/(t + 237.3))
  es <- (esat(Tmax) + esat(Tmin))/2
  delta <- 4098*esat(Tmean)/(Tmean + 237.3)^2
  P <- 101.3*((293 - 0.0065*FAO56_ELEV_M)/293)^5.26
  gamma <- 0.665e-3 * P
  Rns <- (1 - FAO56_ALBEDO) * Rs
  frac <- pmin(pmax(Rs/Rso, 0.25), 1.0)
  Rnl <- FAO56_SIGMA * (((Tmax + 273.16)^4 + (Tmin + 273.16)^4)/2) *
         (0.34 - 0.14*sqrt(pmax(ea, 0))) * (1.35*frac - 0.35)
  Rn <- Rns - Rnl
  u2 <- FAO56_U2
  data.frame(
    ET0_pm = pmax((0.408*delta*Rn + gamma*900/(Tmean + 273)*u2*(es - ea)) / (delta + gamma*(1 + 0.34*u2)), 0),
    ET0_hs = pmax(0.0023 * (Ra/2.45) * (Tmean + 17.8) * sqrt(pmax(Tmax - Tmin, 0)), 0))
}
