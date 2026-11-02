# Definition of regions and masks for altimetry data processing
# date:     September 2025
# author:   Cristina Martí-Solana (IMEDEA CSIC-UIB)

import numpy   as np
import xarray  as xr


def compute_max_ice_mask(dataset, lonmin, lonmax, latmin, latmax):
    """Compute a permanent ice mask from the maximum flag_ice over all time.

    A grid cell is considered permanently ice-affected if ``flag_ice == 1``
    on **any** timestep within the provided dataset (full study period).
    The resulting 2-D Boolean DataArray can be reused across all months so
    that a consistent, study-period-wide mask is applied everywhere.

    Parameters
    ----------
    dataset : xarray.Dataset
        Full (possibly lazy/dask) dataset containing ``flag_ice``.
    lonmin, lonmax, latmin, latmax : float
        Region boundaries matching the rest of the pipeline.

    Returns
    -------
    ice_mask : xarray.DataArray (bool, 2-D: latitude × longitude)
        True where the cell is ice-affected at least once in the period.
    """
    region = dataset.sel(
        longitude=slice(lonmin, lonmax),
        latitude=slice(latmin, latmax),
    )
    if 'flag_ice' not in region:
        raise ValueError("'flag_ice' variable not found in dataset.")

    if 'time' in region['flag_ice'].dims:
        ice_mask = (region['flag_ice'].max(dim='time') > 0).compute()
    else:
        ice_mask = (region['flag_ice'] > 0).compute()

    return ice_mask


def preprocessing_region(dataset, lonmin, lonmax, latmin, latmax, ice_mask=None):
    """
    Preprocess the dataset to select a specific region and apply ice mask.

    Ice masking strategy — maximum ice mask:
        If a pre-computed 2-D ``ice_mask`` (from ``compute_max_ice_mask``)
        is supplied it is applied directly — this is the recommended approach
        for the full pipeline so that a single, study-period-wide mask is
        used consistently for long-term means, monthly EKE, and KE.

        If ``ice_mask`` is not supplied, a per-dataset maximum mask is
        derived on-the-fly: a grid cell is masked for the *entire* period
        if ``flag_ice == 1`` on **any** timestep within that dataset.

    Parameters
    ----------
    dataset : xarray.Dataset
        The input dataset containing altimetry data.
    lonmin, lonmax, latmin, latmax : float
        Region boundaries.
    ice_mask : xarray.DataArray or None
        Optional pre-computed 2-D boolean mask (True = ice-contaminated).
        Produced by ``compute_max_ice_mask``.

    Returns
    -------
    xarray.Dataset
        Region-selected dataset with ice cells set to NaN.
    """
    # Select the region based on provided boundaries
    region_dataset = dataset.sel(
        longitude=slice(lonmin, lonmax),
        latitude=slice(latmin, latmax)
    )

    if ice_mask is not None:
        # Use the supplied study-period-wide mask (recommended)
        region_dataset = region_dataset.where(~ice_mask)
    elif 'flag_ice' in region_dataset:
        # Fallback: derive max mask from the dataset itself
        if 'time' in region_dataset['flag_ice'].dims:
            ice_ever = region_dataset['flag_ice'].max(dim='time') > 0
        else:
            ice_ever = region_dataset['flag_ice'] > 0
        region_dataset = region_dataset.where(~ice_ever)

    return region_dataset


def compute_acc_mask_from_ssh(mean_adt, lat, lon, ssh_south=-0.6, ssh_north=0.2):
    """Create a 2-D boolean mask for the ACC band defined by SSH contours.

    The ACC is approximated as the region where the time-mean absolute dynamic
    topography (ADT) lies between ``ssh_south`` and ``ssh_north``.

    Default contour values follow the Sokolov & Rintoul (2009, JGR-Oceans)
    SSH-contour approach for the AVISO/CMEMS DUACS all-sat product:

        ssh_south ≈ −0.6 m  (near the Southern ACC Front / SACCF)
        ssh_north ≈  0.2 m  (north of the Sub-Antarctic Front / SAF)

    These values enclose the dynamically active part of the ACC and are
    consistent with the range reported in Park et al. (2009) and
    Gille (2014) for the same altimetry product.  They can be tuned via
    the keyword arguments — e.g. use ``ssh_south=-0.8`` to restrict to
    the core ACC jets.

    Parameters
    ----------
    mean_adt : 2-D array-like or xr.DataArray  (latitude × longitude)
        Time-mean absolute dynamic topography [m].
    lat : 1-D array-like
        Latitude coordinate [degrees_north], must match the first axis of
        ``mean_adt``.
    lon : 1-D array-like
        Longitude coordinate [degrees_east], must match the second axis.
    ssh_south : float, optional
        Lower SSH contour defining the southern ACC boundary [m].
        Default: −0.6 m.
    ssh_north : float, optional
        Upper SSH contour defining the northern ACC boundary [m].
        Default: 0.2 m.

    Returns
    -------
    acc_mask : xr.DataArray (bool, latitude × longitude)
        ``True`` where the cell lies **inside** the ACC band.
        Apply as ``data.where(acc_mask)`` to set outside-ACC cells to NaN,
        or pass to ``preprocessing_region`` as an additional spatial mask.

    References
    ----------
    Sokolov, S., & Rintoul, S. R. (2009). Circumpolar structure and
        distribution of the Antarctic Circumpolar Current fronts: 2.
        Variability and relationship to sea surface height.
        *Journal of Geophysical Research: Oceans*, 114, C11019.
        https://doi.org/10.1029/2008JC005248
    """
    lat = np.asarray(lat)
    lon = np.asarray(lon)

    if hasattr(mean_adt, 'values'):
        adt_np = mean_adt.values
    else:
        adt_np = np.asarray(mean_adt, dtype=np.float64)

    acc_mask = (
        np.isfinite(adt_np)
        & (adt_np >= ssh_south)
        & (adt_np <= ssh_north)
    )

    coords = {'latitude': lat, 'longitude': lon}
    return xr.DataArray(
        acc_mask,
        dims=['latitude', 'longitude'],
        coords=coords,
        name='acc_mask',
        attrs={
            'long_name': 'ACC SSH-contour mask',
            'description': (
                f'True inside ACC band: {ssh_south} m ≤ ADT ≤ {ssh_north} m. '
                'Sokolov & Rintoul (2009) SSH-contour approach.'
            ),
            'ssh_south_m': ssh_south,
            'ssh_north_m': ssh_north,
        },
    )


def filter_fronts_by_acc_mask(df, acc_mask):
    """Keep only frontal detections that fall inside the ACC SSH mask.

    Performs a nearest-neighbour lookup on the ``acc_mask`` grid for every
    row in ``df`` identified by ``front_lat`` / ``front_lon``.

    Parameters
    ----------
    df : pd.DataFrame
        Front-detection DataFrame with columns ``front_lat`` and ``front_lon``.
    acc_mask : xr.DataArray (bool, latitude × longitude)
        ACC mask produced by :func:`compute_acc_mask_from_ssh`.

    Returns
    -------
    pd.DataFrame
        Filtered copy of ``df`` containing only inside-ACC fronts.
    """
    lats   = acc_mask.latitude.values
    lons   = acc_mask.longitude.values
    mask_np = acc_mask.values  # bool (lat, lon)

    def _in_acc(row):
        ilat = int(np.argmin(np.abs(lats - row['front_lat'])))
        ilon = int(np.argmin(np.abs(lons % 360 - row['front_lon'] % 360)))
        return bool(mask_np[ilat, ilon])

    inside = df.apply(_in_acc, axis=1)
    return df[inside].copy()