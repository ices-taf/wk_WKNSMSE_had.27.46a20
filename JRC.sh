
# Scenario A
R CMD BATCH --vanilla --quiet '--args iters=999 years=20 nblocks=200 par_env=2 n_workers=200 HCRoption=1' Run_MSE_gridsearch_A_extra.R output/reports/OMBaseline_MPbase_full_HCRA_gridsearch_999_extra.Rout

echo "[A]"

# Scenario C
R CMD BATCH --vanilla --quiet '--args iters=999 years=20 nblocks=200 par_env=2 n_workers=200 HCRoption=3' Run_MSE_gridsearch_C_extra.R output/reports/OMBaseline_MPbase_full_HCRC_gridsearch_999_extra.Rout

echo "[C]"

# Scenario AD
R CMD BATCH --vanilla --quiet '--args iters=999 years=20 nblocks=200 par_env=2 n_workers=200 HCRoption=4' Run_MSE_gridsearch_AD_extra.R output/reports/OMBaseline_MPbase_full_HCRAD_gridsearch_999_extra.Rout

echo "[AD]"

# Scenario A - Alt OM
R CMD BATCH --vanilla --quiet '--args iters=1000 years=20 nblocks=200 par_env=2 n_workers=200 HCRoption=1' Run_MSE_AltOM_A.R output/reports/OMBaseline_MPbase_full_HCRA_AltOM_1000.Rout

echo "[A - Alt OM]"

# Scenario C - Alt OM
R CMD BATCH --vanilla --quiet '--args iters=1000 years=20 nblocks=200 par_env=2 n_workers=200 HCRoption=3' Run_MSE_AltOM_C.R output/reports/OMBaseline_MPbase_full_HCRC_AltOM_1000.Rout

echo "[C - Alt OM]"

# Scenario AD - Alt OM
R CMD BATCH --vanilla --quiet '--args iters=1000 years=20 nblocks=200 par_env=2 n_workers=200 HCRoption=4' Run_MSE_AltOM_AD.R output/reports/OMBaseline_MPbase_full_HCRAD_AltOM_1000.Rout

echo "[AD - Alt OM]"
