#!/usr/bin/env python
# coding: utf-8

# If inputing url as a json then we must have filename metadata included as:
#{"1588269600000": "http://tds.renci.org:8080/thredds//dodsC/2020/nam/2020043018/hsofs/hatteras.renci.org/ncfs-dev-hsofs-nam-master/nowcast/maxele.63.nc"}
# If simply inputting a raw url (not expected to be a common approach) then the utime value will be set to 0000

import os
import sys
import re
import time
#import logging
import datetime as dt
#import numpy.ma as ma
import pandas as pd
import numpy as np
#from pylab import *
import matplotlib.tri as Tri
#import matplotlib.pyplot as plt
import netCDF4

import rasterio as rio
from rasterio.transform import from_origin
from rasterio.plot import show
import geopandas as gpd

from utilities.utilities import utilities as utilities

utilities.log.info('Begin TIF generation')
utilities.log.info('netCDF4 Version = {}'.format(netCDF4.__version__))
utilities.log.info('Pandas Version = {}'.format(pd.__version__))
utilities.log.info('rasterio Version = {}'.format(rio.__version__))
utilities.log.info('Geopandas Version = {}'.format(gpd.__version__))

# define url functionality 
# http://tds.renci.org:8080/thredds/dodsC/2020/nam/2020012706/hsofs/hatteras.renci.org/ncfs-dev-hsofs-nam-master/namforecast/maxele.63.nc

def checkInputVar(v):
    allowable_vars = ['zeta_max', 'vel_max', 'inun_max']
    if v.lower() in allowable_vars: return True
    return False

def get_adcirc_grid(nc):
    agdict = {}
    agdict['lon'] = nc.variables['x'][:]
    agdict['lat'] = nc.variables['y'][:]
    # nv is the triangle list
    agdict['ele'] = nc.variables['element'][:, :] - 1
    agdict['latmin'] = np.mean(nc.variables['y'][:])  # needed for scaling lon/lat plots
    agdict['depth'] = nc.variables['depth'][:]
    return agdict

def get_adcirc_slice(nc, v, it=None):
    advardict = {}
    var = nc.variables[v]
    if re.search('max', v):
        var_d = var[:]  # the actual data
    else:
        var_d = var[it, :]  # the actual data
    var_d[var_d.mask] = np.nan
    advardict['data'] = var_d
    return advardict

def get_interpolation_target(gridname=None, yamlfile=os.path.join(os.path.dirname(__file__), '..', 'config', 'main.yml')):
    """
    fetch geotiff grid parameters from the yaml
    Returns:
        targetgrid dict
        target crs 
    """
    try:
        config = utilities.load_config(yamlfile)['REGRID']
    except:
        utilities.log.error("REGRID: load_config yaml failed")

    if not gridname:
        gridname = 'DEFAULT'
    if gridname not in config.keys():
        gridname = 'DEFAULT'

    targetgrid = {'Latitude': [config[gridname]['upperleft_la']],
                  'Longitude': [config[gridname]['upperleft_lo']],
                  'res': config[gridname]['res'],
                  'nx': config[gridname]['nx'],
                  'ny': config[gridname]['ny']}
    return targetgrid, config[gridname]['adcirc_crs'], config[gridname]['target_crs']

def fetch_XY_fromGeopandas(gdf):
    """
    Simply grab the current lat,lon values from the gdp object
    """
    xtemp = gdf['geometry'].x
    ytemp = gdf['geometry'].y
    return xtemp, ytemp


def computeInundation(advardict, agdict):
    """
    Compute the inundation depth (water amount above ADCIRC's ground level).
    FEMA needs this in ft.
    :param advardict:
    :param agdict:
    :return: masked inundation array
    """
    v = advardict['data']
    d = agdict['depth']
    v = np.where(d > 1, 0, v)
    d = np.where(d > 0, 0, d)
    # conv m to ft
    inun = (v + d) * 3.2808
    # return inundation depth masked out over "open water" (>1)
    return np.where(agdict['depth'] > 1, np.nan, inun)

# Define geopandas processors
# project grid coords, before making Triangulation object
def construct_geopandas(agdict, targetepsg):
    utilities.log.info('Computing GeoPandas DF from ADCIRC grid')
    df_Adcirc = pd.DataFrame(
        {'Latitude': agdict['lat'],
        'Longitude': agdict['lon']})

    t0 = time.time()
    gdf = gpd.GeoDataFrame(
        df_Adcirc, geometry=gpd.points_from_xy(agdict['lon'], agdict['lat']))

    # init crs is LonLat, WGS84
    adcircepsg = agdict['crs']
    utilities.log.info('Adding {} crs to initial GDF'.format(adcircepsg))
    gdf.crs = {'init': adcircepsg}
    utilities.log.info('Converting GDF from {} to {}'.format(adcircepsg,targetepsg))
    gdf = gdf.to_crs({'init': targetepsg})
    utilities.log.info('Time to create GDF was {}'.format(time.time()-t0))
    # utilities.log.debug('GDF data set {}'.format(gdf))
    return gdf

# project interpolation grid to target crs
def compute_geotiff_grid(targetgrid, adcircepsg, targetepsg):
    """
    Results:
        meshdict. Values for upperleft_x, upperleft_y, x,y,xx,yy,xxm,yym
    """
    df_target = pd.DataFrame(data=targetgrid)
    gdf_target = gpd.GeoDataFrame(
        df_target, geometry=gpd.points_from_xy(df_target.Longitude, df_target.Latitude))

    # init projection is LonLat, WGS84
    gdf_target.crs = {'init': adcircepsg}

    # convert to "targetepsg"
    utilities.log.info('Converting GDF from {} to {}'.format(adcircepsg,targetepsg))
    gdf_target = gdf_target.to_crs({'init': targetepsg})

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
    return {'uplx': upperleft_x,
            'uply': upperleft_y,
            'x': x,
            'y': y,
            'xx': xx,
            'yy': yy,
            'xxm': xxm,
            'yym': yym}

# get ADCIRC grid parts;  this need only be done once, as it can be time-consuming over the network
def extract_url_grid(url):
    """
    Extract ADCIRC grid parts from THREDDS dataset
    Results:
    """
    nc = netCDF4.Dataset(url)
    # print(nc.variables.keys())
    agdict = get_adcirc_grid(nc)
    return nc, agdict

def write_tif(meshdict, zi_lin, targetgrid, targetepsg, filename='test.tif'):
    """
    Construct the new TIF file and store it to disk in filename
    """
    xm0, ym0 = meshdict['uplx'], meshdict['uply']
    transform = from_origin(xm0 - targetgrid['res'] / 2, ym0 + targetgrid['res'] / 2, targetgrid['res'], targetgrid['res'])
    utilities.log.debug('TIF transform {}'.format(transform))
    nx = targetgrid['nx']
    ny = targetgrid['ny']
    md = {'crs': targetepsg,
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

def fetchGridName(nc):
    """
    Return a (hopefully) unique grid name, based on the content of
    the netCDF file global attribute agrid, stripping out all spec chars
    except _
    """
    return re.sub(r'[^A-Za-z0-9_]+', '', getattr(nc,'agrid'))

#################################################################

def main(args):
    """
    Entry point for computing geotiff files from an ADCIRC mnetCDF file
    """
    utilities.log.info(args)
    varname = args.varname
    writePNG = False
    filename = '.'.join([varname,'tiff']) if args.filename is None else args.filename
    png_filename = '.'.join([varname,'png']) if args.png_filename is None else args.png_filename

    # Add in option to simply upload a url

    if not checkInputVar(varname):
        utilities.log.error('Variable {} not yet supported.'.format(varname))

    utilities.log.info('Start ADCIRC2Geotiff')
    main_yaml_file =  os.path.join(os.path.dirname(__file__), '..', 'config', 'main.yml') 
    main_config = utilities.load_config(yaml_file=main_yaml_file)

    raster_yaml_file = args.rconfigfile
    raster_config = utilities.load_config(yaml_file=raster_yaml_file)

    # if s3path not passed in, disable SEND2AWS
    if not args.s3path:
        main_config['S3']['SEND2AWS'] = False

    if main_config['S3']['SEND2AWS']:

        from utilities.s3_utilities import utilities as s3_utilities
        s3_resource = s3_utilities.s3
        utilities.log.debug(s3_resource)
        utilities.log.debug(s3_utilities.config)

        thisBucket = s3_utilities.config['S3_UPLOAD_Main_Bucket']
        thisRegion = s3_utilities.config['region_name']

        utilities.log.info('Bucket={}'.format(thisBucket))
        utilities.log.info('Region={}'.format(thisRegion))
        if not s3_utilities.bucket_exists(thisBucket):
            res = s3_utilities.create_bucket(thisRegion, thisBucket)
            utilities.log.info('Bucket {} created.'.format(thisBucket))
        else:
            utilities.log.info('Bucket {} already exists.'.format(thisBucket))

    url = args.url

    rootdir = utilities.fetchBasedir(main_config['DEFAULT']['RDIR'])

    targetgrid, adcircepsg, targetepsg = get_interpolation_target(gridname=args.gridname, yamlfile=raster_yaml_file)

    t0 = time.time()
    nc, agdict = extract_url_grid(url)
    agdict['crs'] = adcircepsg
    #utilities.log.info(f'Reading URL and Building ADCIRC grid took {time.time()-t0} secs')

    # Fetch grid name for building gdf filename
    gridname = args.gridname
    if not args.gridname:
        gridname = fetchGridName(nc)
    utilities.log.info('ADCIRC grid name is {}'.format(gridname))

    # Construct a geopandas object on the input URL grid
    t0 = time.time()
    gdf_pklfile = '{}.gdf.pkl'.format(gridname)

    # Need to check in specified dirs
    f = os.path.join(args.pkldir, gdf_pklfile)
    if not os.path.exists(f):
        gdf = construct_geopandas(agdict, targetepsg)
        if not os.path.exists(args.pkldir): os.makedirs(args.pkldir)
        utilities.writePickle(gdf, filename=f)
        utilities.log.info('Wrote Geopandas file to {}'.format(f))
    else:
        utilities.log.info('{} exists.  Using it...'.format(f))
        gdf = pd.read_pickle(f)

    # Extract the lat,lon values of the current gdf object
    xtemp, ytemp = fetch_XY_fromGeopandas(gdf)
    utilities.log.info('Extracted X and Y from the current (input) GDF')

    # Build Triangulate object for interpolating the input geopandas object
    t0 = time.time()
    tri = Tri.Triangulation(xtemp, ytemp, triangles=agdict['ele'])
    utilities.log.info('Tri interpolation object took {} secs to compute'.format(time.time()-t0))

    # Set up grid for geotiff data
    t0 = time.time()
    meshdict = compute_geotiff_grid(targetgrid, adcircepsg, targetepsg)
    utilities.log.info('compute_geotiff took {} secs'.format(time.time() - t0))
    xxm, yym = meshdict['xxm'], meshdict['yym']

    orig_filename = filename
    orig_png_filename = png_filename

    # filename = os.path.join(rootdir, orig_filename)
    filename = orig_filename
    png_filename = os.path.join(rootdir, orig_png_filename)

    nc = netCDF4.Dataset(url)
    if varname == 'inun_max':
        advardict = get_adcirc_slice(nc, 'zeta_max')
        # compute inundation and replace advardict['data']
        advardict['data'] = computeInundation(advardict, agdict)
    else:
        advardict = get_adcirc_slice(nc, varname)

    vmin = np.nanmin(advardict['data'])
    vmax = np.nanmax(advardict['data'])
    utilities.log.info('Min/Max in ADCIRC Slice: {}, {}'.format(vmin,vmax))

    # construct the interpolator
    t0 = time.time()
    interp_lin = Tri.LinearTriInterpolator(tri, advardict['data'])
    utilities.log.info('Finished linearInterpolator in {} secs'.format(time.time()-t0))

    # zi_lin is the interpolated data that form the rasterized data down below
    zi_lin = interp_lin(xxm, yym)
    utilities.log.debug('zi_lin {}'.format(zi_lin))

    utilities.log.info('Outputting tiff file {}'.format(filename))
    write_tif(meshdict, zi_lin, targetgrid, targetepsg, filename)

    if main_config['S3']['SEND2AWS']:
        resp = s3_utilities.upload(thisBucket, args.s3path, filename)
        if not resp:
            utilities.log.info('Upload to s3://{}:/{}/{} failed.'.format(thisBucket,args.s3path,args.filename))
        else:
            utilities.log.info('Upload to s3://{}:/{}/{} succeeded.'.format(thisBucket,args.s3path,args.filename))
        pass

    #if png_filename is not None:
    #    utilities.log.info('Outputting png file {}'.format(png_filename))
    #    write_png(filename, png_filename)

    utilities.log.info('Finished')

if __name__ == '__main__':

    from argparse import ArgumentParser
    import sys
    parser = ArgumentParser()

    parser.add_argument('--pkldir', action='store', dest='pkldir', default='pklfiles',
                        help='String: Directory of pre-computed grid pkl files.')
    parser.add_argument('--png_filename', action='store', dest='png_filename', default=None,
                        help='String: png output file name.')
    parser.add_argument('--tif_filename', action='store', dest='filename', default=None,
                        help='String: tiff output file name.')
    parser.add_argument('--varname', action='store', dest='varname', default='zeta_max',
                        help='String: zeta_max, vel_max, or inun_max')
    parser.add_argument('--urljson', action='store', dest='urljson', default=None,
                        help='String: Filename with a json of urls to loop over.')
    parser.add_argument('--url', action='store', dest='url', default=None,
                        help='String: url.')
    parser.add_argument('--rasterconfigfile', action='store', dest='rconfigfile', default=None,
                        help='config yml for raster parameters')
    # add a user-supplied grid name to override the computed grid name from the agrid global attr in
    # the ADCIRC netCDF file.  ASGS will know the correct gridname at runtime, so pass it in via --gridname
    # otherwise, it is up to the user to ensure that whatever the gridname determined automatically, there needs to
    # be a set of parameters in the config file that matches the gridname
    parser.add_argument('--gridname', action='store', dest='gridname', default=None,
                        help='String: ADCIRC gridname to use for caching pkl file and getting raster parameters')

    parser.add_argument('--s3path', action='store', dest='s3path', type=str,
                        help='String: object path in s3 bucket.')

    args = parser.parse_args()
    sys.exit(main(args))
