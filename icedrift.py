
# Functions for importing and processing NSIDC Polar Pathfinder sea ice drift data (NSIDC-0116, doi: 10.5067/INAWUWO7QH7B)

# DEPENDENCIES:
import xarray as xr
import numpy as np
import numpy.ma as ma
import cartopy
import cartopy.crs as ccrs
from datetime import datetime, timedelta
from metpy.units import units
import pandas as pd
import os


# FUNCTIONS:
#---------------------------------------------------------------------
def convert_vectors(lon = [], u_EASE = [], v_EASE = []):
    
    """Convert u/v (polar EASE grid) velocity components from NSIDC Polar Pathfinder data to 
    eastward and northward vector components.
    
INPUT:
- u_EASE: M x N grid of along-x velocity component on EASE grid
- v_EASE: M x N grid of along-y velocity component on EASE grid
    (u_EASE: toward the right on the grid)
    (v_EASE: upward (toward the top) on the grid)
- lon: M x N grid of longitude grid (0 to 360) associated with u, v

OUTPUT:
- u:  M x N grid of east component of velocity
- v:  M x N grid of north component of velocity

Latest recorded update:
10-25-2024
    """
    # convert EASE grid vector components to northward, eastward
    # ref: https://nsidc.org/support/how/how-convert-horizontal-and-vertical-components-east-and-north
    #------------------------------------------------------------------
    u = u_EASE * np.cos(lon/180*np.pi)  +  v_EASE * np.sin(lon/180*np.pi)
    v = -u_EASE * np.sin(lon/180*np.pi)  +  v_EASE * np.cos(lon/180*np.pi)
                
    return u, v


def grab_projection(ds, quiet = True):
    
    """Grab projection info from NSIDC PP sea ice drift (NSIDC-0116) data (doi: 10.5067/MPYG15WAA4WX)

INPUT: 
- ds: sea ice drift data opened with xarray
- quiet: bool, whether or not to supress print statements (default: True)

OUTPUT:
- ice_projection: cartopy projection from data projection info

Latest recorded update:
10-25-2024
    """
    
    spatial = ds.crs
    
    # grab parameters from crs spatial attributes
    semimajor = float(spatial.proj4text[spatial.proj4text.find('+a=')+3:].split(' ')[0])
    semiminor = float(spatial.proj4text[spatial.proj4text.find('+b=')+3:].split(' ')[0])
    central_longitude = float(spatial.proj4text[spatial.proj4text.find('+lon_0=')+7:].split(' ')[0])
    central_latitude = float(spatial.proj4text[spatial.proj4text.find('+lat_0=')+7:].split(' ')[0])

    if not quiet:
        print(f'>>> data provided in {spatial.grid_mapping_name} projection from the {spatial.long_name}')
        print(f'  - semi_major_axis: {semimajor}')
        print(f'  - semi_minor_axis: {semiminor}')
        print(f'  - central_longitude: {central_longitude}')
        print(f'  - central_latitude: {central_latitude}')
        print(f'  - proj4text: {spatial.proj4text}')

    # create ice projection from info
    projection = ccrs.LambertAzimuthalEqualArea(central_longitude=central_longitude, 
                                                central_latitude=central_latitude,
                                                globe=ccrs.Globe(semimajor_axis = semimajor, semiminor_axis = semiminor))
    
    return projection


def generate_filename(year, main_path, filenametype):
        """Generate filename for NSIDC Polar Pathfinder data based on year of data requested.
        
        INPUT:
        - year: year of data to import
        - main_path: directory where PPD files are locally stored
        - filenametype: naming convention for PPD files (as in 'icemotion_daily_nh_25km_{}0101_{}1231_v4.1.nc' where {} will be replaced with year of date)
        
        OUTPUT:
        - file: string with path to file corresponding to year of data requested
        """
        
        # incomplete year of data in 1978
        if year == 1978: 
            file = main_path + 'icemotion_daily_nh_25km_19781101_19781231_v4.1.nc'
        else:
            file = main_path + filenametype.format(year,year)
        
        return file

def identify_existing_files_dates(main_path):

    """Identify which files exist in the main directory and which dates they include.
    """
    
    # determine which files exist in the main directory
    existing_files = [file for file in os.listdir(main_path) if file.endswith('.nc')]
    if len(existing_files) == 0:
        print('No files found in the specified directory. Please check the path.')
        
    # determine which dates these files include
    existing_dates = np.array([], dtype=np.datetime64)

    for file in existing_files:
        datei = datetime.strptime(file.split('_')[4], '%Y%m%d') 
        datef = datetime.strptime(file.split('_')[5], '%Y%m%d') 
        existing_dates = np.append(existing_dates, pd.date_range(datei, datef))
        
    return existing_files, pd.to_datetime(existing_dates)
        

def open_local_file(dates,
                    main_path = '/Volumes/Jewell_EasyStore/NSIDC-0116_PPdrift/', 
                    filenametype = 'icemotion_daily_nh_25km_{}0101_{}1231_v4.1.nc', 
                    crop = [0,None,0,None],
                    include_units = False):

    """Import NSIDC Polar Pathfinder (sea ice drift NSIDC-0116, doi:10.5067/INAWUWO7QH7B) lats, lons, u, v cropped to within given lat/lon range.
    Note that any desired dates not found in local directory will be ignored, but listed in "missing_dates" key.

INPUT:
- dates: single datetime object, or list, array of desired datetimes
- main_path: directory where PPD files are locally stored.
- filenametype: naming convention for PPD files (default: 'icemotion_daily_nh_25km_{}0101_{}1231_v4.1.nc' where {} will be replaced with year of dt_obj)
- crop: indices along dim1 ("i"), dim2 ("j") to crop to area of interest [ai, bi, aj, bj]
- include_units: bool, whether or not to return data with units

OUTPUT:
Dictionary "data" containing:
- lon: M x N longitude grid (0 to 360) associated with u, v
- lat: M x N latitude grid associated with u, v
- e: M x N grid of eastward vector components of ice drift
- n: M x N grid of northward vector components of ice drift
- xx: M x N grid of x values of EASEgrid, corresponding to u_EASE, v_EASE
- yy: M x N grid of y values of EASEgrid, corresponding to u_EASE, v_EASE
- u: M x N grid of along-x component of ice drift
- v: M x N grid of along-y component of ice drift
- error: M x N grid of estimated error variance (ice motion error measure)
- proj: cartopy projection from PP drift data projection info
- ds: xarray data frame containing data from year of date
- missing_dates: list of dates that were not found in the local directory

Latest recorded update:
03-19-2025
    """

    # determine which files exist in the main directory
    existing_files, existing_dates = identify_existing_files_dates(main_path)
    if len(existing_files) == 0:
        print('No files found in the specified directory. Please check the path.')

    # strip time info from dates
    if isinstance(dates, datetime):
        dates = [dates]
    dates_no_hours = np.array([date-timedelta(hours=date.hour) for date in dates])

    # if we have just one date, add it to a list
    if isinstance(dates_no_hours, datetime):
        dates_no_hours = [dates_no_hours]
    
    # determine if any dates are missing from local files
    include_dates = []
    missing_dates = []
    for date in dates_no_hours:
        if date not in existing_dates:
            missing_dates.append(date)
        else:
            include_dates.append(date)

    include_dates = list(set(include_dates))
    missing_dates = list(set(missing_dates))

    # print(include_dates)
    # determine years of data to import
    #----------------------------------
    years = list(set([date.year for date in include_dates]))

    # open files
    #-----------
    # just one:
    if len(years) == 1:
        # open file corresponding to date's year
        file = generate_filename(years[0], main_path, filenametype)
        ds = xr.open_dataset(file)
        ds.close()

    # or several:
    else:
        # open file corresponding to date's year
        files = []
        for year in years:
            files.append(generate_filename(year, main_path, filenametype))
        ds = xr.open_mfdataset(files)
        ds.close()

    # convert from CFtime to standard calendar
    # and grab dates provided
    ds = ds.convert_calendar('standard').sel(time=include_dates)

    # save all to dict
    data = {}
    data['ds'] = ds
    data['proj'] = grab_projection(ds) # cartopy projection

    # projected drift components and error variance
    # remove time info if only one date
    if len(include_dates) == 1:
        data['u'] = ds.u.values[0,:,:]
        data['v'] = ds.v.values[0,:,:]
        data['error'] = ds.icemotion_error_estimate.values[0,:,:]
    else:
        data['u'] = ds.u.values
        data['v'] = ds.v.values
        data['error'] = ds.icemotion_error_estimate.values
    
    # projected coords
    data['xx'], data['yy'] = np.meshgrid(ds.x.values, ds.y.values)

    # lat/lon coords(strip time dimension)
    data['lat'] = ds.latitude.sel(time = ds.time[0]).values
    data['lon'] = ds.longitude.sel(time = ds.time[0]).values
    data['lon'][data['lon']<0]+=360
    
    # convert projected coordinates to east, north velocities
    e, n = convert_vectors(lon = data['lon'], 
                           u_EASE = data['u'], 
                           v_EASE = data['v'])
    data['e'] = e
    data['n'] = n
        
    # crop data values here
    ai, bi, aj, bj = crop[0], crop[1], crop[2], crop[3]

    data['xx'] = data['xx'][ai:bi, aj:bj]
    data['yy'] = data['yy'][ai:bi, aj:bj]
    data['lat'] = data['lat'][ai:bi, aj:bj]
    data['lon'] = data['lon'][ai:bi, aj:bj]
    data['missing_dates'] = missing_dates

    if len(data['u'].shape) == 3:
        data['u'] = data['u'][:,ai:bi, aj:bj]
        data['v'] = data['v'][:,ai:bi, aj:bj]
        data['e'] = data['e'][:,ai:bi, aj:bj]
        data['n'] = data['n'][:,ai:bi, aj:bj]
        data['error'] = data['error'][:,ai:bi, aj:bj]
    else:
        data['u'] = data['u'][ai:bi, aj:bj]
        data['v'] = data['v'][ai:bi, aj:bj]
        data['e'] = data['e'][ai:bi, aj:bj]
        data['n'] = data['n'][ai:bi, aj:bj]
        data['error'] = data['error'][ai:bi, aj:bj]

    if include_units:
        data['u'] = data['u'] * units(ds['u'].units)
        data['v'] = data['v'] * units(ds['v'].units)
        data['e'] = data['e'] * units(ds['u'].units)
        data['n'] = data['n'] * units(ds['v'].units)
        data['error'] = data['error'] * units('cm/s') # according to documentation
        data['xx'] = data['xx'] * units(ds['x'].units)
        data['yy'] = data['yy'] * units(ds['y'].units)
        data['lat'] = data['lat'] * units(ds['latitude'].units)
        data['lon'] = data['lon'] * units(ds['longitude'].units)


    return data
