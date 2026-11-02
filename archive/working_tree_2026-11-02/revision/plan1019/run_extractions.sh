#!/usr/bin/env bash
# Extractions needed for the remaining plan items. They need your logins and cannot run from Claude's workspace
# (Copernicus Marine is not reachable from there). Run from the repository root, in the environment you used for
# extract_monthly.py (copernicusmarine installed and logged in). Every step is restartable.
set -e
cd "$(dirname "$0")/../.."      # repository root

# 1) 1993-2000 for the full Martinez-Moreno / Hogg period (65-35 S), both products  (~0.8 GB)
python revision/extract_monthly.py --product allsat --start-year 1993 --end-year 2000
python revision/extract_monthly.py --product twosat --start-year 1993 --end-year 2000

# 2) Antarctic shelf and slope for the Hogg/Hughes sea-level index (80-58 S)
python revision/extract_monthly.py --product allsat --start-year 2001 --end-year 2025 --lat-min -80 --lat-max -58 --out-dir revision/data/monthly_allsat_south
python revision/extract_monthly.py --product twosat --start-year 1993 --end-year 2025 --lat-min -80 --lat-max -58 --out-dir revision/data/monthly_twosat_south
#    bathymetry to 80 S (needs a global ETOPO 2022 / GEBCO file you download yourself):
# python revision/plan1019/h0_prepare_bathymetry_80S.py /path/to/ETOPO_2022_bedrock.nc revision/data/bathymetry_80s.nc

# 3) Earlier product version (two-sat vDT2018, C3S; needs ~/.cdsapirc and `pip install cdsapi`)
# python revision/plan1019/extract_c3s_twosat.py --version vdt2018 --start-year 2001 --end-year 2020
# python revision/plan1019/extract_c3s_twosat.py --version vdt2024 --start-year 2001 --end-year 2020   # same-source control

# Afterwards (or ask Claude to run them):
# M1_Y0=1993 python revision/plan1019/m1_comparison.py revision revision/plan1019/outputs/m1_1993
# python revision/plan1019/h0_antarctic_slope_index.py revision/data/monthly_twosat_south revision/data/bathymetry_80s.nc revision/plan1019/outputs/h0/h0_twosat.csv 1993 2025
# python revision/plan1019/h0_antarctic_slope_index.py revision/data/monthly_allsat_south revision/data/bathymetry_80s.nc revision/plan1019/outputs/h0/h0_allsat.csv 2001 2025
