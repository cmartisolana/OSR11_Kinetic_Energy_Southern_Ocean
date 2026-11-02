# Half-power point EKE (from Chambers et al., 2018) for front detection
# Cristina Martí-Solana, October 2025
#
# EKE definition:
#   EKE = ½ · mean_t( (u - ū_lt)² + (v - v̄_lt)² )
#   where ū_lt, v̄_lt are the long-term means over the FULL study period.
#
# Ice masking:
#   A permanent ice mask (max flag_ice over the full study period) is computed
#   once and applied consistently to the long-term mean and every monthly field.

import argparse
import os
import xarray as xr
import numpy  as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Internal imports
from utils.chambers_fronts import detect_fronts_from_dataset
from utils.regions_and_masks import preprocessing_region, compute_max_ice_mask

# ── Configuration ──
# Input: daily L4 gridded altimetry (CMEMS product SEALEVEL_GLO_PHY_L4_MY_008_047,
# dataset cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.125deg_P1D), organised
# locally as <input_dir>/<YYYY>/<MM>/*.nc. Download via the copernicusmarine
# toolbox (see README) and point --input-dir (or ALTIMETRY_DIR) at the folder.
_parser = argparse.ArgumentParser(description="Half-power point EKE preprocessing.")
_parser.add_argument(
    "--input-dir",
    default=os.environ.get("ALTIMETRY_DIR"),
    help="Directory with daily altimetry NetCDF files organised as YYYY/MM/*.nc "
         "(default: $ALTIMETRY_DIR).",
)
_args = _parser.parse_args()
if not _args.input_dir:
    raise SystemExit("No input directory given. Use --input-dir or set $ALTIMETRY_DIR.")
input_dir = _args.input_dir
LONMIN, LONMAX = -180, 180
LATMIN, LATMAX = -65, -35

year_folders = sorted([f for f in os.listdir(input_dir) if f not in ['.DS_Store']])

# ── Pass 1: Collect ALL files for ice mask + long-term mean ──
all_files = []
for yr in year_folders:
    yr_path = os.path.join(input_dir, yr)
    for mo in sorted(os.listdir(yr_path)):
        if mo in ['.DS_Store']: continue
        mo_path = os.path.join(yr_path, mo)
        for f in sorted(os.listdir(mo_path)):
            if f.endswith('.nc'):
                all_files.append(os.path.join(mo_path, f))

print(f"Total files in study period: {len(all_files)}")

def _keep_uv_ice(ds):
    return ds[['ugos', 'vgos', 'flag_ice']]

# Permanent ice mask — True wherever flag_ice==1 on ANY day in the full period
ds_ice = xr.open_mfdataset(all_files, combine='by_coords',
                            preprocess=lambda ds: ds[['flag_ice']])
ice_mask = compute_max_ice_mask(ds_ice, LONMIN, LONMAX, LATMIN, LATMAX)
ds_ice.close()
print(f"Permanent ice mask: {int(ice_mask.sum())} cells masked "
      f"({100*float(ice_mask.mean()):.1f}%)")

# Long-term mean velocity (with permanent ice mask applied)
ds_all = xr.open_mfdataset(all_files, combine='by_coords', preprocess=_keep_uv_ice)
ds_all = preprocessing_region(ds_all, LONMIN, LONMAX, LATMIN, LATMAX, ice_mask=ice_mask)
print("Computing long-term mean velocities...")
u_ltm = ds_all['ugos'].mean(dim='time', skipna=True).compute()
v_ltm = ds_all['vgos'].mean(dim='time', skipna=True).compute()
ds_all.close()
print("Long-term mean computed.")

# ── Pass 2: Monthly processing ──
for year in year_folders[-5:]:   # process only the last 5 years for testing
    year_path = os.path.join(input_dir, year)
    month_paths = sorted([os.path.join(year_path, m)
                          for m in os.listdir(year_path) if m not in ['.DS_Store']])

    for month in month_paths:
        files = sorted([f for f in os.listdir(month) if f.endswith('.nc')])
        files_path = [os.path.join(month, f) for f in files]

        ds = xr.open_mfdataset(files_path, combine='by_coords', preprocess=_keep_uv_ice)

        # Apply permanent ice mask
        ds = preprocessing_region(ds, LONMIN, LONMAX, LATMIN, LATMAX, ice_mask=ice_mask)

        # ── EKE from full-period anomalies ──
        # u'(t) = u(t) - ū_lt;  EKE = ½·mean_t(u'² + v'²)
        u_anom = ds['ugos'] - u_ltm
        v_anom = ds['vgos'] - v_ltm
        eke = 0.5 * (u_anom**2 + v_anom**2).mean(dim='time', skipna=True)

        # Monthly mean and KE
        ds_mean = ds.mean(dim='time')
        ds_mean['ke']  = 0.5 * (ds_mean['ugos']**2 + ds_mean['vgos']**2)
        ds_mean['eke'] = eke
        ds_mean = ds_mean.load()

        print()
        print(f"Processing year {year}, month {os.path.basename(month)}")

        results = detect_fronts_from_dataset(ds_mean, use_numba=True, sample_every=1)

        # Filter to valid latitude range
        filtered_results = {}
        for lon, envs in results.items():
            valid = [e for e in envs if LATMIN <= e['lat_mid'] <= LATMAX]
            if valid:
                filtered_results[lon] = valid

        n_fronts = sum(len(envs) for envs in filtered_results.values())

        if n_fronts == 0:
            print(f"No fronts detected for {year}/{os.path.basename(month)}")
            continue

        # Build output arrays
        front_lons = np.zeros(n_fronts)
        front_lats = np.zeros(n_fronts)
        envelope_lat_south = np.zeros(n_fronts)
        envelope_lat_north = np.zeros(n_fronts)
        envelope_lon = np.zeros(n_fronts)
        peak_lat = np.zeros(n_fronts)
        peak_val = np.zeros(n_fronts)
        width_km = np.zeros(n_fronts)

        lat_vals = ds_mean.latitude.values
        idx = 0
        for lon, envs in filtered_results.items():
            for e in envs:
                front_lons[idx] = e['lon_mid']
                front_lats[idx] = e['lat_mid']
                envelope_lat_south[idx] = lat_vals[e['isouth']]
                envelope_lat_north[idx] = lat_vals[e['inorth']]
                envelope_lon[idx] = lon
                peak_lat[idx] = e['peak_lat']
                peak_val[idx] = e['peak_val']
                width_km[idx] = e['width_km']
                idx += 1

        ds_output = xr.Dataset(
            {
                'ke':  (('latitude', 'longitude'), ds_mean['ke'].values),
                'eke': (('latitude', 'longitude'), ds_mean['eke'].values),
                'front_lon':          (['front'], front_lons),
                'front_lat':          (['front'], front_lats),
                'envelope_lon':       (['front'], envelope_lon),
                'envelope_lat_south': (['front'], envelope_lat_south),
                'envelope_lat_north': (['front'], envelope_lat_north),
                'peak_lat':           (['front'], peak_lat),
                'peak_val':           (['front'], peak_val),
                'width_km':           (['front'], width_km),
            },
            coords={
                'latitude':  ds_mean.latitude.values,
                'longitude': ds_mean.longitude.values,
                'front':     np.arange(n_fronts),
            },
            attrs={
                'title':             'Southern Ocean Fronts - Chambers (2018) Method',
                'description':       'Half-power points and 25% envelope boundaries',
                'method':            'Chambers, D. P. (2018)',
                'eke_reference':     'long-term mean over full study period',
                'ice_masking':       'permanent (max flag_ice over full study period)',
                'threshold':         '0.02 m²/s²',
                'min_width':         '100 km',
                'peak_fraction':     '0.25',
                'year':              year,
                'month':             os.path.basename(month),
                'n_fronts_detected': n_fronts,
            }
        )
        ds_output['ke'].attrs  = {
            'long_name': 'Mean Kinetic Energy', 'units': 'm²/s²',
            'description': 'KE = ½(ū² + v̄²); monthly-mean absolute geostrophic velocities',
        }
        ds_output['eke'].attrs = {
            'long_name': 'Eddy Kinetic Energy', 'units': 'm²/s²',
            'description': "EKE = ½·mean_t((u-ū_lt)²+(v-v̄_lt)²); anomalies from full-period mean",
        }
        ds_output['front_lon'].attrs         = {'long_name': 'Half-power point longitude', 'units': 'degrees_east'}
        ds_output['front_lat'].attrs         = {'long_name': 'Half-power point latitude',  'units': 'degrees_north'}
        ds_output['envelope_lon'].attrs      = {'long_name': 'Envelope longitude', 'units': 'degrees_east'}
        ds_output['envelope_lat_south'].attrs = {'long_name': 'Southern boundary (25% of peak)', 'units': 'degrees_north'}
        ds_output['envelope_lat_north'].attrs = {'long_name': 'Northern boundary (25% of peak)', 'units': 'degrees_north'}
        ds_output['peak_lat'].attrs          = {'long_name': 'Peak EKE latitude', 'units': 'degrees_north'}
        ds_output['peak_val'].attrs          = {'long_name': 'Peak EKE value',    'units': 'm²/s²'}
        ds_output['width_km'].attrs          = {'long_name': 'Envelope width',    'units': 'km'}

        output_dir = f"./outputs/monthly_half-power_points/{year}"
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, f"{os.path.basename(month)}.nc")
        ds_output.to_netcdf(output_file)
        print(f"Saved {n_fronts} fronts → {output_file}")