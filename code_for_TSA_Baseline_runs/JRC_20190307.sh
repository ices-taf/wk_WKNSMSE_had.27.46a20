
# Scenario A - Alt OM

/opt/R/bin/R CMD BATCH --vanilla --quiet '--args iters=1000 years=20 nblocks=200 par_env=2 n_workers=200 HCRoption=1' Run_MSE_AltOM_A.R output/reports/OMBaseline_MPbase_full_HCRA_AltOM_1000.Rout

echo "[A]"

# Scenario C - Alt OM

/opt/R/bin/R CMD BATCH --vanilla --quiet '--args iters=1000 years=20 nblocks=200 par_env=2 n_workers=200 HCRoption=3' Run_MSE_AltOM_C.R output/reports/OMBaseline_MPbase_full_HCRC_AltOM_1000.Rout

echo "[C]"

# Scenario AD - Alt OM

/opt/R/bin/R CMD BATCH --vanilla --quiet '--args iters=1000 years=20 nblocks=200 par_env=2 n_workers=200 HCRoption=4' Run_MSE_AltOM_AD.R output/reports/OMBaseline_MPbase_full_HCRAD_AltOM_1000.Rout

echo "[AD]"
