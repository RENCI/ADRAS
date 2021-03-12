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
import matplotlib.pyplot as plt
import netCDF4

import rasterio as rio
from rasterio.transform import from_origin
from rasterio.plot import show
import geopandas as gpd

from utilities.utilities import utilities as utilities

utilities.log.info(f'Begin TIF generation')
utilities.log.info(f'netCDF4 Version = {netCDF4.__version__}')
utilities.log.info(f'Pandas Version = {pd.__version__}')
utilities.log.info(f'rasterio Version = {rio.__version__}')
utilities.log.info(f'Geopandas Version = {gpd.__version__}')

# define url functionality 
# http://tds.renci.org:8080/thredds/dodsC/2020/nam/2020012706/hsofs/hatteras.renci.org/ncfs-dev-hsofs-nam-master/namforecast/maxele.63.nc
ensname = 'namforecast'


def checkInputVar(v):
    allowable_vars = ['zeta_max', 'vel_max', 'inun_max']
    if v.lower() in allowable_vars: return True
    return False


# Deprecated
def get_url(dict_cfg, varname, datestr, hourstr, yearstr, enstag):
    """
    Build a URL for finding ADCIRC information.
    Parameters:
        datestr: "YYYYmmdd"
        hourstr: "hh"
        yearstr: "YYYY"
        input dict.
    Returns:
        url as a full path string
    """
    # dstr = dt.datetime.strftime(d, "%Y%m%d%H")
    url = dict_cfg['baseurl']+dict_cfg['dodsCpart'] % \
        (yearstr,
         datestr+hourstr,
         dict_cfg['AdcircGrid'],
         dict_cfg['Machine'],
         dict_cfg['Instance'],
         enstag,
         varname)
    return url


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
    utilities.log.info(f'Adding {adcircepsg} crs to initial GDF')
    gdf.crs = {'init': adcircepsg}
    utilities.log.info(f'Converting GDF from {adcircepsg} to {targetepsg}')
    gdf = gdf.to_crs({'init': targetepsg})
    utilities.log.info(f'Time to create GDF was {time.time()-t0}')
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
    utilities.log.info(f'Converting GDF from {adcircepsg} to {targetepsg}')
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


#  Not used anymore
def construct_url(varname, geo_yamlfile=os.path.join(os.path.dirname(__file__), '..', 'config', 'main.yml')):
    """
    Assembles several method into an aggregate method to grab parameters
    from the yaml and construct a url. This is skipped if the user
    specified a URL on input.
    """
    geo_yml= utilities.load_config(geo_yamlfile)
    varnamedict = geo_yml['VARFILEMAP']
    varfile = varnamedict[varname]
    utilities.log.info('map dict {}'.format(varnamedict))
    # Specify time parameters of interest
    timedict = geo_yml['TIME']
    doffset = timedict['doffset']
    hoffset = timedict['hoffset']

    # Set time for url
    thisdate = dt.datetime.utcnow() + dt.timedelta(days=doffset) + dt.timedelta(hours=hoffset)
    cyc = "%02d" % (6 * int(thisdate.hour / 6))     # Hour/cycle specification
    dstr = dt.datetime.strftime(thisdate, "%Y%m%d")
    ystr = dt.datetime.strftime(thisdate, "%Y")  # NOTE slight API change to get_url

    # Fetch url
    urldict = geo_yml['ADCIRC']
    utilities.log.info('url dict {}'.format(urldict))

    # TODO: we will also need to elevate "namforecast" to be a default/input parameter, since this
    # can take several different values that depend in the ASGS configuration.
    # However, namforecast is a reasonable default
    url = get_url(urldict, varfile, dstr, str(cyc), ystr, ensname)
    utilities.log.info('Validated url {}'.format(url))
    utilities.log.info('Datetime {}'.format(dstr))
    utilities.log.debug('Constructed URL {}'.format(url))
    return url, dstr, cyc


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


# Assemble some optional plot methods
def plot_triangular_vs_interpolated(meshdict, varname, tri, zi_lin, advardict):
    """
    A rough plotting routine generally used for validation studies.
    The real data results are the tif fle that gets generated later
    """
    xm0, ym0 = meshdict['uplx'], meshdict['uply']
    x, y = meshdict['x'], meshdict['y']
    xx, yy = meshdict['xx'], meshdict['yy']
    xxm, yym = meshdict['xxm'], meshdict['yym']
    nlev = 11
    vmin = np.floor(np.nanmin(zi_lin))
    vmax = np.ceil(np.nanmax(zi_lin))
    levels = np.linspace(0., vmax, nlev+1)
    utilities.log.info('Levels are {}, vmin {}, vmax {}'.format(levels, vmin, vmax))
    #
    v = advardict['data']
    utilities.log.debug('nanmin {}, nammix {}'.format(np.nanmin(v),np.nanmax(v)))
    #
    cmap = plt.cm.get_cmap('jet', 8)
    # Start the plots
    fig, ax = plt.subplots(1, 2, figsize=(20, 20), sharex=True, sharey=True)
    # tcf = ax[0].tricontourf(tri, v,cmap=plt.cm.jet,levels=levels)
    if True:
        # fig, ax = plt.subplots(1, 2, figsize=(20, 20), sharex=True, sharey=True)
        tcf = ax[0].tripcolor(tri, v, cmap=cmap, vmin=vmin, vmax=vmax, shading='flat')
        ax[0].set_aspect('equal')
        ax[0].plot(xx[0, ], yy[0, ], color='k', linewidth=.25)
        ax[0].plot(xx[-1, ], yy[-1, ], color='k', linewidth=.25)
        ax[0].plot(xx[:, 0], yy[:, 0], color='k', linewidth=.25)
        ax[0].plot(xx[:, -1], yy[:, -1], color='k', linewidth=.25)
        fig.colorbar(tcf, ax=ax[0], orientation='horizontal')
        ax[0].set_title('ADCIRC {}'.format(varname), fontsize=14)
    pcm = ax[1].pcolormesh(xxm, yym, zi_lin, cmap=cmap,  shading='faceted', vmin=vmin, vmax=vmax)
    ax[1].set_aspect('equal')
    ax[1].plot(xx[0, ], yy[0, ], color='k')
    ax[1].plot(xx[-1, ], yy[-1, ], color='k')
    ax[1].plot(xx[:, 0], yy[:, 0], color='k')
    ax[1].plot(xx[:, -1], yy[:, -1], color='k')
    ax[1].set_xlim([min(x), max(x)])
    ax[1].set_ylim([min(y), max(y)])
    fig.colorbar(pcm, ax=ax[1], orientation='horizontal')
    ax[1].set_title('Interpolated ADCIRC {}'.format(varname), fontsize=14)
    plt.show()


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


def write_png(filenametif, filenamepng):
    """
    Construct the new PNG file  from the tif file
    Translate will automatically create and save the new file,
    """
    from osgeo import gdal
    scale = '-scale min_val max_val'
    options_list = [
        '-ot Byte',
        '-of PNG',
        scale
    ]
    options_string = " ".join(options_list)
    gdal.Translate(filenamepng,
           filenametif,
           options=options_string)
    utilities.log.info('Storing a PNG with the name {}'.format(filenamepng))


def plot_tif(filename):
    """
    Read TIF file that has been previously generated
    """
    dataset = rio.open(filename)
    band1 = dataset.read(1, masked=True)
    show(band1, cmap='jet')
    msk = dataset.read_masks(1)


def plot_png(filename):
    """
    Read TIF file that has been previously generated
    """
    import matplotlib.image as mpimg

    img = mpimg.imread(filename)
    plt.imshow(img)
    plt.show()


def fetchGridName(nc):
    """
    Return a (hopefully) unique grid name, based on the content of
    the netCDF file global attribute agrid, stripping out all spec chars
    except _
    """
    return re.sub(r'[^A-Za-z0-9_]+', '', getattr(nc,'agrid'))

#################################################################
# urlinput='http://tds.renci.org:8080/thredds//dodsC/2020/nam/2020042912/hsofs/hatteras.renci.org/ncfs-dev-hsofs-nam-master/namforecast/maxele.63.nc'


def main(args):
    """
    Entry point for computing geotiff files from an ADCIRC mnetCDF file
    """

    utilities.log.info(args)
    varname = args.varname
    showInterpolatedPlot = args.showInterpolatedPlot
    showRasterizedPlot = args.showRasterizedPlot
    showPNGPlot = args.showPNGPlot
    writePNG = False
    showGDALPlot = False # Tells GDAL to load the tif and display it
    filename = '.'.join([varname,'tiff']) if args.filename is None else args.filename
    png_filename = '.'.join([varname,'png']) if args.png_filename is None else args.png_filename

    # Add in option to simply upload a url

    if not checkInputVar(varname):
        utilities.log.error(f'Variable {varname} not yet supported.')

    utilities.log.info('Start ADCIRC2Geotiff')
    yaml_file = os.path.join(os.path.dirname(__file__), '..', 'config', 'main.yml')
    main_config = utilities.load_config(yaml_file=yaml_file)

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

        utilities.log.info(f'Bucket={thisBucket}')
        utilities.log.info(f'Region={thisRegion}')

        if not s3_utilities.bucket_exists(thisBucket):
            res = s3_utilities.create_bucket(thisRegion, thisBucket)
            utilities.log.info(f'Bucket {thisBucket} created.')
        else:
            utilities.log.info(f'Bucket {thisBucket} already exists.')

    if args.urljson is not None:
        setGetURL = False
        if not os.path.exists(args.urljson):
            utilities.log.error('urljson file not found.')
            sys.exit(1)
        utilities.log.info('Read data from supplied urljson')
        urls = utilities.read_json_file(args.urljson)
    elif args.url is not None:
        # If here we still need to build a dict for ADCIRC
        url = args.url
        dte='manual'  # The times will be determined from the real data
        urls={dte:url}
        utilities.log.info('Explicit URL provided {}'.format(urls))
    else:
        utilities.log.error('No Proper URL specified')

    rootdir = utilities.fetchBasedir(main_config['DEFAULT']['RDIR'])

    ################################################################
    #
    # We pick the very first url simply to set up the initial/final grids and build an interpolator
    # Assumes all other urls have the same grid
    ################################################################

    targetgrid, adcircepsg, targetepsg = get_interpolation_target(gridname=args.gridname)
    # use the first url to get the ADCIRC grid parts, assumes all urls point to the same
    # ADCIRC grid...
    url = list(urls.values())[0]

    t0 = time.time()
    nc, agdict = extract_url_grid(url)
    agdict['crs'] = adcircepsg
    utilities.log.info(f'Reading URL and Building ADCIRC grid took {time.time()-t0} secs')

    # Fetch grid name for building gdf filename
    gridname = args.gridname
    if not args.gridname:
        gridname = fetchGridName(nc)
    utilities.log.info(f'ADCIRC grid name is {gridname}')

    # Construct a geopandas object on the input URL grid
    t0 = time.time()
    gdf_pklfile = f'{gridname}.gdf.pkl'

    # Need to check in specified dirs
    f = os.path.join(rootdir, 'pklfiles', gdf_pklfile)
    if not os.path.exists(f):
        gdf = construct_geopandas(agdict, targetepsg)
        # gdf.to_pickle(f)
        if not os.path.exists('pklfiles'): os.makedirs('pklfiles')
        # if filename is passed in, all other arguments are ignored
        utilities.writePickle(gdf, filename=f, rootdir=rootdir, fileroot=gdf_pklfile, subdir='', iometadata='')
        utilities.log.info(f'Wrote Geopandas file to {f}')
        utilities.log.info('Construct geopandas object took {} secs'.format(time.time() - t0))
    else:
        utilities.log.info(f'{f} exists.  Using it...')
        gdf = pd.read_pickle(f)

    # Extract the lat,lon values of the current gdf object
    xtemp, ytemp = fetch_XY_fromGeopandas(gdf)
    utilities.log.info('Extracted X and Y from the current (input) GDF')

    # Build Triangulate object for interpolating the input geopandas object
    t0 = time.time()
    tri = Tri.Triangulation(xtemp, ytemp, triangles=agdict['ele'])
    utilities.log.info('Tri interpolation object took {} secs to compute'.format(time.time()-t0))

    # Setup desired grid for geotiff data

    t0 = time.time()
    meshdict = compute_geotiff_grid(targetgrid, adcircepsg, targetepsg)
    utilities.log.info('compute_geotiff took {} secs'.format(time.time() - t0))
    xxm, yym = meshdict['xxm'], meshdict['yym']

    orig_filename = filename
    orig_png_filename = png_filename

    # For now assume only one url was provided

    for utime, url in urls.items():
        iometadata = utime # Grabs the key and uses it to differentiate urls
        # filename = "{}/{}_{}".format(rootdir,iometadata,orig_filename)
        # filename = os.path.join(rootdir, orig_filename)
        filename = orig_filename
        # png_filename = "{}/{}_{}".format(rootdir,iometadata,orig_png_filename)
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
        utilities.log.info(f'Min/Max in ADCIRC Slice: {vmin}/{vmax}')

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
                utilities.log.info(f'Upload to s3://{thisBucket}:/{args.s3path}/{args.filename} failed.')
            else:
                utilities.log.info(f'Upload to s3://{thisBucket}:/{args.s3path}/{args.filename} succeeded.')

            pass

        if writePNG:
            utilities.log.info('Outputting png file {}'.format(png_filename))
            write_png(filename, png_filename)

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

    parser.add_argument('--tif_filename', action='store', dest='filename', default=None,
                        help='String: tiff output file name will be prepended by new path. Must include extension')

    parser.add_argument('--png_filename', action='store', dest='png_filename', default=None,
                        help='String: png output file name will be prepended by new path. Must include extension')

    parser.add_argument('--showInterpolatedPlot', action='store_true',
                        help='Boolean: Display the comparison of Triangular and interpolated plots')

    parser.add_argument('--showRasterizedPlot', action='store_true',
                        help='Boolean: Display the generated and saved tif plot')

    parser.add_argument('--showPNGPlot', action='store_true',
                        help='Boolean: Display the generated and saved png plot')

    parser.add_argument('--varname', action='store', dest='varname', default='zeta_max',
                        help='String: zeta_max, vel_max, or inun_max')

    parser.add_argument('--urljson', action='store', dest='urljson', default=None,
                        help='String: Filename with a json of urls to loop over.')

    parser.add_argument('--url', action='store', dest='url', default=None,
                        help='String: url.')

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
