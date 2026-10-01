# soil_collars.R -- soil collars outside the modelled stand.
# One definition, used by 05_model/01_load_and_prep_data.R (excludes them from the
# soil model) and zenodo/01_compile_datasets.R (records why they are not in training).
# Both collars sit in the wetland margin outside the censused stand, so the soil model,
# which is applied only inside the stand, is not trained on them.
OUT_OF_STAND_COLLARS <- c("WS_1-1", "WS_1-2")
