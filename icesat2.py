



# ICESat-2 L4 Monthly Gridded Sea Ice Thickness, Version 3
# Data set id: IS2SITMOGR4
# DOI: 10.5067/ZCSU8Y5U1BQW

from datetime import datetime, timedelta
import xarray as xr
import numpy as np
import cartopy
import cartopy.crs as ccrs
from metpy.units import units
import os



def construct_filename(date, VersionID = '003'):

    # construct filename
    file = f'IS2SITMOGR4_01_{date.year}{date.strftime("%m")}_006_{VersionID}.nc'

    return file


def extract_variables(ds):
    """Extract variables from dataset."""
    data = {}
    data['ds'] = ds
    data['proj'] = grab_projection(ds)

    for var in list(ds.variables.keys()):
        if var not in ['x', 'y', 'crs', 'ice_type', 'sea_ice_conc', 'region_mask']:
            data[var] = ds[var].values * units(ds[var].units)
            
    xx, yy = np.meshgrid(ds['x'].values, ds['y'].values)
    data['xx'], data['yy'] = xx * units(ds['x'].units), yy * units(ds['x'].units)
    
    data['ice_type'] = ds.ice_type.values
    data['sea_ice_conc'] = 100*ds.sea_ice_conc.values * units('%')
    data['region_mask'] = ds.region_mask     
    
    data['region_mask_flags'] = [int(flag) for flag in ds.region_mask.flag_values.split(', ')]
    
    meanings = ds.region_mask.flag_meanings.split(' ')
    while '' in meanings:
        meanings.remove('')
    
    data['region_mask_flag_meanings'] = meanings
    
    return data


def open_local_file(date,
                    main_path = '/Volumes/Seagate_Jewell/KenzieStuff/NSIDC-IS2SITMOGR4/', 
                    filenametype = 'icemotion_daily_nh_25km_{}0101_{}1231_v4.1.nc', 
                    include_units = False, quiet = True):

    """Import NSIDC Polar Pathfinder (sea ice drift NSIDC-0116, doi:10.5067/INAWUWO7QH7B) lats, lons, u, v cropped to within given lat/lon range.

INPUT:
- date: single datetime object 
- main_path: directory where PPD files are locally stored.
- filenametype: naming convention for PPD files (default: 'icemotion_daily_nh_25km_{}0101_{}1231_v4.1.nc' where {} will be replaced with year of dt_obj)
- include_units: bool, whether or not to return data with units
- quiet: bool, whether or not to supress print statements (default: True)

OUTPUT:
Dictionary "data" containing:
- proj: cartopy projection from PP drift data projection info
- ds: xarray data frame containing data from year of date
- all other variables included in file

Latest recorded update:
10-31-2024
    """

    # open file corresponding to date's year and month
    filename = construct_filename(date, VersionID = '003')
    ds = xr.open_dataset(os.path.join(main_path, filename))

    if not quiet:
        print(f'>>> opening {filename}')

    # Extract variables from dataset
    data = extract_variables(ds)

    # remove units if desired
    if not include_units:
        for key in data.keys():
            if key not in ['ds','proj','crs', 'ice_type', 'region_mask', 'region_mask_flags', 'region_mask_flag_meanings']:
                data[key] = data[key].magnitude


    return data



def grab_projection(ds, quiet = True):

    """Grab projection info from NSIDC L4 ICESat-2 monthly gridded data 
    (data set ID: IS2SITMOGR4, DOI: 10.5067/ZCSU8Y5U1BQW)

INPUT: 
- ds: data opened with xarray
- quiet: bool, whether or not to supress print statements (default: True)

OUTPUT:
- projection: cartopy projection from data projection info

Latest recorded update:
10-31-2024
    """
    
    spatial = ds.crs
    
    # grab parameters from crs spatial attributes
    semimajor = spatial.attrs['semimajor_radius']
    semiminor = spatial.attrs['semiminor_radius']
    central_longitude = spatial.attrs['straight_vertical_longitude_from_pole']
    standard_parallel = spatial.attrs['standard_parallel']

    if not quiet:
        print(f'>>> data provided in {spatial.grid_mapping_name} projection from the {spatial.long_name}')
        print(f'  - semi_major_axis: {semimajor}')
        print(f'  - semi_minor_axis: {semiminor}')
        print(f'  - central_longitude: {central_longitude}')
        print(f'  - central_latitude: {standard_parallel}')
        print(f'  - proj4text: {spatial.proj4text}')

    # create ice projection from info
    projection = ccrs.NorthPolarStereo(central_longitude=central_longitude, 
                                       true_scale_latitude = standard_parallel,
                                       globe=ccrs.Globe(semimajor_axis = semimajor, semiminor_axis = semiminor))
    
    return projection