#!/usr/bin/env python
# coding: utf-8

# If inputting urls as a json then we must have filename metadata included as:
#{"1588269600000": "http://tds.renci.org:8080/thredds//dodsC/2020/nam/2020043018/hsofs/hatteras.renci.org/ncfs-dev-hsofs-nam-master/nowcast/maxele.63.nc"}
# If simply inputting a raw url (not expected to be a common approach) then the utime value will be set to 0000

import os
import sys
import time
#import re
import datetime as dt
import numpy.ma as ma
import pandas as pd
import xarray as xr

from pylab import *
import matplotlib.tri as Tri
from matplotlib import pyplot
import netCDF4

from utilities.utilities import utilities as utilities

import rasterio as rio
from rasterio.transform import from_origin
from rasterio.plot import show
import geopandas as gpd

utilities.log.info('Begin TIF generation')
utilities.log.info("netCDF4 Version = {}".format(netCDF4.__version__))
utilities.log.info("Pandas Version = {}".format(pd.__version__))
utilities.log.info("Matplotlib Version = {}".format(matplotlib.__version__))
utilities.log.info("rasterio Version = {}".format(rio.__version__))
utilities.log.info("Geopandas Version = {}".format(gpd.__version__))

# define url functionality 
# http://tds.renci.org:8080/thredds/dodsC/2020/nam/2020012706/hsofs/hatteras.renci.org/ncfs-dev-hsofs-nam-master/namforecast/maxele.63.nc

# TODO: Jeff, the original intent of the regex that checkEnumuation replaced
# was to catch the "max" in the variable name.  If max is in the name,
# then the ADCIRC solution component being "geotiffed" in this run
# does NOT contain time, and thus the time index does not exist.
# We do not want to do that check in this "def" way.
# We may want validate the variable name parameter, but rather like this:
# possible_vars = ['zeta_max',
#                  'vel_max',
#                  'inun_max']
# if var not in possible_vars:
#     error...


def checkEnumuation(v):
    if v.lower() in ('zeta_max'):
        return True
    if v.lower() in ('vel_max'):
        return True
    if v.lower() in ('inun_max'):
        return True
    return False

# Define some basic grid functionality and defaults
def validate_url(url):  # TODO
    """
    """
    pass

def get_adcirc_grid(nc):
    """
    """
    agdict = {}
    agdict['lon'] = nc.variables['x'][:]
    agdict['lat'] = nc.variables['y'][:]
    # nv is the triagle list
    agdict['ele'] = nc.variables['element'][:,:] - 1
    agdict['latmin'] = np.mean(nc.variables['y'][:])  # needed for scaling lon/lat plots
    agdict['depth'] = nc.variables['depth'][:]
    return agdict

def get_adcirc_grid_from_ds(ds):
    """
    """
    agdict = {}
    agdict['lon'] = ds['x'][:]
    agdict['lat'] = ds['y'][:]
    # nv is the triagle list
    agdict['ele'] = ds['element'][:,:] - 1
    agdict['latmin'] = np.mean(ds['y'][:])  # needed for scaling lon/lat plots
    agdict['depth'] = ds['depth'][:]
    return agdict

def get_adcirc_time(nc):
    """
    """
    atdict = {}
    temp = nc.variables['time']
    dates = num2date(temp[:], temp.units, only_use_python_datetimes=False)
    #dates = [date.strftime('%Y%m%d%H%M') for date in dates]
    dates = [pd.to_datetime(date.strftime('%Y%m%d%H%M')) for date in dates]
    atdict['time'] = dates
    return atdict

def get_adcirc_time_from_ds(ds):
    """
    """
    return {'time': ds['time']}

def get_adcirc_slice(nc,v,it=None):
    """
    """
    advardict = {}
    var = nc.variables[v]
    if re.search('max', v):
        var_d = var[:] # the actual data
    else:
        var_d = var[it,:] # the actual data
    var_d[var_d.mask] = np.nan
    advardict['var'] = var_d
    return advardict

def get_adcirc_slice_from_ds(ds,v,it=0):
    """
    """
    advardict = {}
    var = ds.variables[v]
    if re.search('max', v):
        var_d = var[:]  # the actual data
    else:
        var_d = var[it, :]  # the actual data
    #var_d[var_d.mask] = np.nan
    advardict['data'] = var_d
    return advardict


def default_inter_grid():
    """
    A simple grid for testing. See also read_inter_grid
    to fetch grid data from the yaml
    Returns:
        targetgrid
        target crs 
    """
    upperleft_lo = -77.
    upperleft_la = 35.7
    res = 100  # resolution in meters
    nx = 1000
    ny = 1000
    crs = 'epsg:6346'  
    targetgrid = {'Latitude': [upperleft_la],
                  'Longitude': [upperleft_lo],
                  'res': res,
                  'nx': nx,
                  'ny': ny}
    return targetgrid, crs


def read_inter_grid_yaml():
    """
    fetch geotiff grid parameters from the yaml
    Returns:
        targetgrid dict
        target crs 
    """
    try:
        config = utilities.load_config()['TARGETGRID']
    except:
        utilities.log.error("TARGETGRID: load_config yaml failed")
    upperleft_lo = config['upperleft_lo']
    upperleft_la = config['upperleft_la']
    res = config['res']
    nx = config['nx']
    ny = config['ny']
    crs = config['coordrefsys']
    targetgrid = {'Latitude': [upperleft_la],
                  'Longitude': [upperleft_lo],
                  'res': res,
                  'nx': nx,
                  'ny': ny,
                  'crs': crs}

    return targetgrid

# Define geopandas processors
# project grid coords, before making Triangulation object
def adcircgrid_to_geopandas(agdict, targetepsg):

    utilities.log.info('Making DataFrame of ADCIRC Grid Coords')
    df_Adcirc = pd.DataFrame({'Latitude': agdict['lat'],
                              'Longitude': agdict['lon']})

    utilities.log.info('Converting to Geopandas DF')
    t0 = time.time()
    gdf = gpd.GeoDataFrame(
        df_Adcirc, geometry=gpd.points_from_xy(agdict['lon'], agdict['lat']))

    # init crs is LonLat, WGS84
    utilities.log.info('Adding WGS84 crs')
    gdf.crs = {'init': 'epsg:4326'}
    utilities.log.info('Converting to {}'.format(targetepsg))
    gdf = gdf.to_crs({'init': targetepsg})

    utilities.log.info('Time to create geopandas was {}'.format(time.time()-t0))
    utilities.log.debug('GDF data set {}'.format(gdf))
    return gdf


# project interpolation grid to target crs
def compute_target_grid(targetgrid, targetepsg):
    """
    Results:
        meshdict. Values for upperleft_x, upperleft_y, x,y,xx,yy,xxm,yym
    """
    df_target = pd.DataFrame(data=targetgrid)
    gdf_target = gpd.GeoDataFrame(
        df_target, geometry=gpd.points_from_xy(df_target.Longitude,
                                               df_target.Latitude))
    # print(gdf_target.crs)

    # init projection is LonLat, WGS84
    gdf_target.crs = {'init': 'epsg:4326'}

    # convert to "targetepsg"
    gdf_target = gdf_target.to_crs({'init': targetepsg})
    # print(gdf_target.crs)

    # compute spatial grid for raster (for viz purposes below)
    upperleft_x = gdf_target['geometry'][0].x
    upperleft_y = gdf_target['geometry'][0].y
    x = np.arange(upperleft_x, upperleft_x+targetgrid['nx']*targetgrid['res'], targetgrid['res'])
    y = np.arange(upperleft_y, upperleft_y-targetgrid['ny']*targetgrid['res'], -targetgrid['res'])
    xx, yy = np.meshgrid(x, y)

    # get centroid coords
    xm = (x[1:] + x[:-1]) / 2
    ym = (y[1:] + y[:-1]) / 2
    xxm, yym = np.meshgrid(xm, ym)
    utilities.log.debug('compute_mesh: lon {}. lat {}'.format(upperleft_x, upperleft_y))
    
    meshdict = {'uplx': upperleft_x,
                'uply': upperleft_y,
                'x': x,
                'y': y,
                'xx': xx,
                'yy': yy,
                'xxm': xxm,
                'yym': yym}

    targetgrid.update(meshdict)

    return targetgrid

# def merge_dicts(self, *dict_args):
#     """
#     Given any number of dicts, shallow copy and merge into a new dict,
#     precedence goes to key value pairs in latter dicts.
#     """
#     result = {}
#     for dictionary in dict_args:
#         result.update(dictionary)
#     return result


def reg_grid_params(self):
    # load list of reg grid params
    return self.config["TARGETGRID"]["RECT"]

# get ADCIRC grid parts;  this need only be done once, as it can be time-consuming over the network
def get_triangulation(agdict, gdf):
    """
    Extract ADCIRC grid parts from THREDDS dataset
    Results:
        tri:
    """
    # get converted ADCIRC node coordinates
    xtemp = gdf['geometry'].x
    ytemp = gdf['geometry'].y
    tri = Tri.Triangulation(xtemp, ytemp, triangles=agdict['ele'])
    return tri



#################################################################
## Start doing some work

# Do we need varname as in input argument ?
# Yes we do make it an enumeration of three possible things.
# Want to be able to simply input a URL
# urlinput='http://tds.renci.org:8080/thredds//dodsC/2020/nam/2020042912/hsofs/hatteras.renci.org/ncfs-dev-hsofs-nam-master/namforecast/maxele.63.nc'
#
def get_target_grid():
    # set up target interpolation grid
    upperleft_lo = -93.95
    upperleft_la = 30.15
    lower_right_lo = -89.5
    lower_right_la = 28.95
    res = 1000  # resolution in meters
    nx = 500
    ny = 200
    targetgrid = {'Latitude': [upperleft_la], 'Longitude': [upperleft_lo], 'res': res, 'nx': nx, 'ny': ny,
                  'crs': 'epsg:6346'}
    return targetgrid


def write_tif(targetgriddict, zi_lin, targetepsg, filename='test.tif'):
    """
    Construct the new TIF file and store it to disk in filename
    """
    xm0, ym0 = targetgriddict['uplx'], targetgriddict['uply']
    transform = from_origin(xm0 - targetgriddict['res'] / 2,
                            ym0 + targetgriddict['res'] / 2,
                            targetgriddict['res'],
                            targetgriddict['res'])
    utilities.log.info('TIF transform {}'.format(transform))

    nx = targetgrid['nx']
    ny = targetgrid['ny']
    crs = targetepsg
    md = {'crs': crs,
          'driver': 'GTiff',
          'height': ny,
          'width': nx,
          'count': 1,
          'dtype': zi_lin.dtype,
          'nodata': -99999,
          'transform': transform}
    # output a geo-referenced tiff
    dst = rio.open(filename, 'w', **md)
    try:
        dst.write(zi_lin, 1)
        utilities.log.info('Wrote TIF file to {}'.format(filename))
    except:
        utilities.log.error('Failed to write TIF file to {}'.format(filename))
    dst.close()

def main(args):
    """

    """
    # utilities.log.info(args)
    # filename = args.filename
    # png_filename = args.png_filename
    # varname = args.varname
    # showInterpolatedPlot = args.showInterpolatedPlot
    # showRasterizedPlot = args.showRasterizedPlot
    # showPNGPlot = args.showPNGPlot
    #
    # showGDALPlot = False # Tells GDAL to load the tif and display it
    # # Add in option to simply upload a url
    #
    # utilities.log.info('Start ADCIRC2Geotiff')
    # config = utilities.load_config()

    outputfilename='test'
    url = args.url
    targetepsg = 'epsg:6346'
    if url is None:
        exit(0)

    t0 = time.time()
    ds = xr.open_dataset(url, drop_variables=['neta', 'nvel'])
    ds.attrs['agrid']

    agdict = get_adcirc_grid_from_ds(ds)
    pklfile='{}.gdf.pkl'.format(ds.attrs['agrid'])

    if os.path.exists(pklfile):
        # use it
        gdf = pd.read_pickle(pklfile)
    else:  # generate/store
       gdf = adcircgrid_to_geopandas(agdict, targetepsg)
       gdf.to_pickle('{}.gdf.pkl'.format(ds.attrs['agrid']))

    tri = get_triangulation(agdict, gdf)

    targetgrid = get_target_grid()
    targetgriddict = compute_target_grid(targetgrid, targetepsg)

    utilities.log.info('compute_geotiff_grid took {} secs'.format(time.time() - t0))

    xxm, yym = targetgriddict['xxm'], targetgriddict['yym']

    advardict = get_adcirc_slice_from_ds(ds, 'zeta_max')
    vmin = np.nanmin(advardict['data'])
    vmax = np.nanmax(advardict['data'])
    utilities.log.info('Min/Max in ADCIRC Slice: {}/{}'.format(vmin, vmax))

    # construct the interpolator
    t0 = time.time()
    interp_lin = Tri.LinearTriInterpolator(tri, advardict['data'])
    utilities.log.info('Finished linearInterpolator in {} secs'.format(time.time()-t0))

    # zi_lin is the interpolated data that form the rasterized data down below
    zi_lin = interp_lin(xxm, yym)
    utilities.log.debug('zi_lin {}'.format(zi_lin))

    t0 = time.time()
    write_tif(targetgriddict, zi_lin, targetepsg, outputfilename)
    utilities.log.info('compute_mesh took {} secs'.format(time.time()-t0))

    t0 = time.time()
    write_png(filename, png_filename)
    utilities.log.info('write_png took {} secs'.format(time.time()-t0))

    if showInterpolatedPlot:
        plot_triangular_vs_interpolated(meshdict, varname, tri, zi_lin, advardict)

    if showRasterizedPlot:
        plot_tif(filename)

    if showPNGPlot:
        plot_png(png_filename)

    if showGDALPlot:
        # Can we reread the file using GDAL ?
        from osgeo import gdal
        ds = gdal.Open(filename).ReadAsArray()
        im = plt.imshow(ds)
        plt.show()

    utilities.log.info('Finished') 


if __name__ == '__main__':
    from argparse import ArgumentParser
    import sys
    parser = ArgumentParser()
    # parser.add_argument('--experiment_name', action='store', dest='experiment_name', default=None,
    #                     help='Highlevel Experiment-tag value')
    # parser.add_argument('--tif_filename', action='store', dest='filename', default='test.tif',
    #                     help='String: tif output file name will be prepended by new path. Must include extension')
    # parser.add_argument('--png_filename', action='store', dest='png_filename', default='test.png',
    #                     help='String: png output file name will be prepended by new path. Must include extension')
    # parser.add_argument('--showInterpolatedPlot', type=str2bool, action='store', dest='showInterpolatedPlot', default=False,
    #                     help='Boolean: Display the comparison of Trangular and interpolated plots')
    # parser.add_argument('--showRasterizedPlot', type=str2bool, action='store', dest='showRasterizedPlot', default=False,
    #                     help='Boolean: Display the generated and saved tif plot')
    # parser.add_argument('--showPNGPlot', type=str2bool, action='store', dest='showPNGPlot', default=False,
    #                     help='Boolean: Display the generated and saved png plot')
    # parser.add_argument('--varname', action='store', dest='varname', default='zeta_max',
    #                     help='String: zeta_max, vel_max, or inun_max')
    parser.add_argument('--url', action='store', dest='url', default=None,
                        help='String: simply input a URL for processing')
    # parser.add_argument('--urljson', action='store', dest='urljson', default=None,
    #                     help='String: Filename with a json of urls to loop over.')
    args = parser.parse_args()
    sys.exit(main(args))


