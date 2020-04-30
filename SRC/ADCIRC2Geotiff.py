#!/usr/bin/env python
# coding: utf-8

# This notebook demonstrates an approach for computing geotiffs from ADCIRC FEM output in netCDF. 

import os, sys, time
#import re
import datetime as dt
import numpy.ma as ma
import pandas as pd

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


def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def checkEnumuation(v):
    if v.lower() in ('zeta_max'):
        return True
    if v.lower() in ('vel_max'):
        return True
    if v.lower() in ('inun_max'):
        return True
    return False

def build_dict_from_yaml():
    """
    Read the config yaml file and populate the dict.
    used downstream to build URLs.
    """
    try:
        config = utilities.load_config()['ADCIRC']
    except:
        utilities.log.error("load_config yaml failed")
    datagrid = config['AdcircGrid']
    machine = config['Machine']
    instance = config['Instance']
    baseurl =  config['baseurl']
    catpart = config['catPart']
    #dodscpart = '/dodsC/%s/nam/%s/%s/%s/%s/%s/%s'
    dodscpart = config['dodsCpart'] 
    fortnumber = config['fortNumber']
    # NOTE the use of dodsCpart has changed wrt ADDA work 
    urldict = {
        'AdcircGrid': datagrid,
        'Machine': machine, 
        'Instance': instance, 
        'baseurl': baseurl,
        'catPart': catpart,
        'dodsCpart': dodscpart,
        'fortNumber': fortnumber}
    return urldict

def get_url(dict_cfg,varname,datestr,hourstr,yearstr,enstag):
    """
    Build a URL for finding ADCIRC information.
    Parameters:
        datastr: "YYYYmmdd"
        hourstr: "hh"
        yearstr: "YYYY"
        input dict from build_dict_from_yaml.
    Returns:
        url as a full path string
    """
    # dstr = dt.datetime.strftime(d, "%Y%m%d%H")
    url = dict_cfg['baseurl']+dict_cfg['dodsCpart'] % (yearstr,
    datestr+hourstr,
    dict_cfg['AdcircGrid'],
    dict_cfg['Machine'],
    dict_cfg['Instance'],
    enstag,
    varname)
    return url

def validate_url(url):  # TODO
    status = True # Innocent until proven guilty
    if not status:
        utilities.log.error("Invalid URL:".format(url))
    return True 

# Define some basic grid functionality and defaults

def get_adcirc_grid(nc):
    agdict = {}
    agdict['lon'] = nc.variables['x'][:]
    agdict['lat'] = nc.variables['y'][:]
    # nv is the triagle list
    agdict['ele'] = nc.variables['element'][:,:] - 1
    agdict['latmin'] = np.mean(nc.variables['y'][:])  # needed for scaling lon/lat plots
    return agdict

def get_adcirc_slice(nc,v,it=None):
    advardict = {}
    var = nc.variables[v]
    if re.search('max', v):
        var_d = var[:] # the actual data
    else:
        var_d = var[it,:] # the actual data
    var_d[var_d.mask] = np.nan
    advardict['data'] = var_d
    return advardict

def default_inter_grid():
    """
    A simple grid for testing. See also read_inter_grid
    to fetch grid data from the yaml
    Returns:
        targetgrid
        targetgridepsg (reference)
    """
    upperleft_lo = -77.
    upperleft_la = 35.7
    res = 100  # resolution in meters
    nx = 1000
    ny = 1000
    coorref = 'epsg:6346'  
    targetgrid={'Latitude': [upperleft_la], 'Longitude': [upperleft_lo],'res': res,'nx': nx,'ny': ny}
    return targetgrid, coorref

def read_inter_grid_yaml():
    """
    fetch grid data from the yaml
    Returns:
        targetgrid
        targetgridepsg (reference)
    """
    try:
        config = utilities.load_config()['REGRID']
    except:
        utilities.log.error("REGRID: load_config yaml failed")
    upperleft_lo = config['upperleft_lo']
    upperleft_la = config['upperleft_la']
    res = config['res']
    nx = config['nx']
    ny = config['ny']
    coorref = config['reference']
    targetgrid={'Latitude': [upperleft_la], 'Longitude': [upperleft_lo],'res': res,'nx': nx,'ny': ny}
    return targetgrid, coorref

# Define geopandas processors
# project grid coords, before making Triangulation object
def construct_geopandas(agdict, targetepsg):
    print('Making DataFrame of ADCIRC Grid Coords')
    df_Adcirc = pd.DataFrame(
        {'Latitude': agdict['lat'],
        'Longitude': agdict['lon']})
    utilities.log.info('Converting to Geopandas DF')
    t0 = time.time()
    gdf = gpd.GeoDataFrame(
        df_Adcirc, geometry=gpd.points_from_xy(agdict['lon'], agdict['lat']))
    print(gdf.crs)
    # init projection is LonLat, WGS84
    utilities.log.info('Adding WGS84 projection information')
    gdf.crs={'init' :'epsg:4326'}
    utilities.log.info('Converting to {}'.format(targetepsg))
    gdf = gdf.to_crs({'init': targetepsg})
    print(gdf.crs)
    xtemp=gdf['geometry'].x
    ytemp=gdf['geometry'].y
    utilities.log.info('Time to create geopandas was {}'.format(time.time()-t0))
    utilities.log.debug('GDF data set {}'.format(gdf))
    return xtemp, ytemp

# project interpolation grid to target crs

def compute_mesh(targetgrid, targetepsg):
    """
    Results:
        meshdict. Values for upperleft_x, upperleft_y, x,y,xx,yy,xxm,yym
    """
    df_target = pd.DataFrame(data=targetgrid)
    gdf_target = gpd.GeoDataFrame(
        df_target, geometry=gpd.points_from_xy(df_target.Longitude, df_target.Latitude))
    # print(gdf_target.crs)
    #init projection is LonLat, WGS84
    gdf_target.crs={'init' : 'epsg:4326' }
    # convert to "targetepsg"
    gdf_target = gdf_target.to_crs({'init': targetepsg})
    # print(gdf_target.crs)
    # compute spatial grid for raster (for viz purposes below)
    upperleft_x=gdf_target['geometry'][0].x
    upperleft_y=gdf_target['geometry'][0].y
    x=np.arange(upperleft_x,upperleft_x+targetgrid['nx']*targetgrid['res'],targetgrid['res'])
    y=np.arange(upperleft_y,upperleft_y-targetgrid['ny']*targetgrid['res'],-targetgrid['res'])
    xx,yy = np.meshgrid(x,y)
    # get centroid coords
    xm=(x[1:] + x[:-1]) / 2
    ym=(y[1:] + y[:-1]) / 2
    xxm,yym = np.meshgrid(xm,ym)
    #print(gdf_target['Longitude'][0], upperleft_x)
    #print(gdf_target['Latitude'][0], upperleft_y)
    utilities.log.debug('compute_mesh: lon {}. lat {}'.format(upperleft_x,upperleft_y)) 
    # also need geometry
    meshdict ={'uplx':upperleft_x, 'uply':upperleft_y,
        'x':x, 'y':y, 'xx':xx, 'yy':yy, 'xxm':xxm, 'yym':yym} 
    return meshdict

# Aggregation of individual methods

def construct_url(varname):
    """
    Assembles several method into an aggreegate method to grab parameters
    from the yaml and construct a url. THis is skipped is the user
    specified a URL on input.
    """
    varnamedict = utilities.load_config()['VARFILEMAP'] 
    varfile = varnamedict[varname]
    utilities.log.info('map dict {}'.format(varnamedict))
    # Specify time parameters of interest
    timedict = utilities.load_config()['TIME'] 
    doffset = timedict['doffset']
    hoffset = timedict['hoffset']
    period = timedict['period']
    utilities.log.info('time dict {}'.format(timedict))
    # Set time for ADCVIC data fetch
    thisdate = dt.datetime.utcnow() + dt.timedelta(days=doffset) + dt.timedelta(hours=hoffset)
    cyc = "%02d"%(period * int(thisdate.hour / period))     # Hour specification   
    dstr = dt.datetime.strftime(thisdate, "%Y%m%d") #
    ystr = dt.datetime.strftime(thisdate, "%Y") # NOTE slight API change to get_url
    # Fetch url
    urldict = build_dict_from_yaml()
    utilities.log.info('url dict {}'.format(urldict))
    url=get_url(urldict,varfile,dstr,str(cyc),ystr,'namforecast')
    utilities.log.info('Validated dict {}'.format(url))
    utilities.log.info('Datetime {}'.format(dstr))
    utilities.log.debug('Constructed URL {}'.format(url))
    return url, dstr, cyc

# get ADCIRC grid parts;  this need only be done once, as it can be time-consuming over the network

def construct_grid(url, varname):
    """
    Construct the grid parts for use by all subsequent plots
    Results:
        advardict:
        tri:
        xtemp,ytemp:
        targetgrid, targetepsg:
    """
    nc = netCDF4.Dataset(url)
    print(nc.variables.keys())
    agdict = get_adcirc_grid(nc)
    advardict = get_adcirc_slice(nc, varname)
    #
    targetgrid, targetepsg = read_inter_grid_yaml()
    # targetgrid, targetepsg = default_inter_grid() # A builtin option but same as yaml example
    xtemp, ytemp = construct_geopandas(agdict, targetepsg)
    t0 = time.time()
    tri = Tri.Triangulation(xtemp, ytemp, triangles=agdict['ele'])
    deltat = time.time()-t0
    vmin=np.nanmin(advardict['data'])
    vmax=np.nanmax(advardict['data'])
    utilities.log.info('Min/Max in ADCIRC Slice: {}/{}'.format(vmin,vmax))
    print("Min/Max in ADCIRC Slice: {}/{}".format(vmin,vmax))
    return xtemp, ytemp, tri, targetgrid, targetepsg, advardict

# Assemble some optionmal plot methods

def plot_triangular_vs_interpolated(meshdict, varname, tri, zi_lin, advardict):
    """
    A rough plotting routine generally used for validation studies.
    The real data results are the tif fle that gets generated later
    """
    xm0,ym0 = meshdict['uplx'],meshdict['uply']
    x,y = meshdict['x'],meshdict['y']
    xx,yy = meshdict['xx'],meshdict['yy']
    xxm,yym = meshdict['xxm'],meshdict['yym']
    nlev=11
    vmin=np.floor(np.nanmin(zi_lin))
    vmax=np.ceil(np.nanmax(zi_lin))
    levels = linspace(0.,vmax,nlev+1)
    utilities.log.info('Levels are {}, vmin {}, vmax {}'.format(levels,vmin, vmax))
    #
    v = advardict['data']
    utilities.log.debug('nanmin {}, nammix {}'.format(np.nanmin(v),np.nanmax(v))) 
    #
    cmap=plt.cm.get_cmap('jet', 8)
    # Start the plots
    fig, ax = plt.subplots(1,2,figsize=(20,20), sharex=True, sharey=True)
    #f, axs = plt.subplots(2,2)
    # tcf = ax[0].tricontourf(tri, v,cmap=plt.cm.jet,levels=levels)
    if True:
        fig, ax = plt.subplots(1,2,figsize=(20,20), sharex=True, sharey=True)
        tcf = ax[0].tripcolor(tri, v, cmap=cmap, vmin=vmin, vmax=vmax, shading='flat')
        ax[0].set_aspect('equal')
        ax[0].plot(xx[0,],yy[0,],color='k',linewidth=.25)
        ax[0].plot(xx[-1,],yy[-1,],color='k',linewidth=.25)
        ax[0].plot(xx[:,0],yy[:,0],color='k',linewidth=.25)
        ax[0].plot(xx[:,-1],yy[:,-1],color='k',linewidth=.25)
        fig.colorbar(tcf,ax=ax[0],orientation='horizontal')
        ax[0].set_title('ADCIRC {}'.format(varname), fontsize=14)
    pcm = ax[1].pcolormesh(xxm, yym, zi_lin, cmap=cmap,  shading='faceted', vmin=vmin, vmax=vmax)
    ax[1].set_aspect('equal')
    ax[1].plot(xx[0,],yy[0,],color='k')
    ax[1].plot(xx[-1,],yy[-1,],color='k')
    ax[1].plot(xx[:,0],yy[:,0],color='k')
    ax[1].plot(xx[:,-1],yy[:,-1],color='k')
    ax[1].set_xlim([min(x),max(x)])
    ax[1].set_ylim([min(y),max(y)])
    fig.colorbar(pcm,ax=ax[1],orientation='horizontal')
    ax[1].set_title('Interpolated ADCIRC {}'.format(varname), fontsize=14)
    plt.show()

def write_tif(meshdict, zi_lin, targetgrid, targetepsg, filename='test.tif'):
    """
    Construct the new TIF file and store it to disk in filename
    """
    xm0,ym0 = meshdict['uplx'],meshdict['uply']
    transform = from_origin(xm0 - targetgrid['res'] / 2, ym0 + targetgrid['res'] / 2, targetgrid['res'], targetgrid['res'])
    utilities.log.info('TIF transform {}'.format(transform))
    print(transform)
    nx = targetgrid['nx']
    ny = targetgrid['ny']
    crs = targetepsg
    md={'crs':crs, 'driver':'GTiff','height':ny, 'width':nx,'count':1,'dtype':zi_lin.dtype,'nodata':-99999,'transform':transform}
    # output a geo-referenced tiff
    dst = rio.open(filename, 'w', **md)
    try:
        dst.write(zi_lin, 1)
        utilities.log.info('Wrote TIF file to {}'.format(filename))
    except:
        utilities.log.error('Failed to write TIF file to {}'.format(filename))
    dst.close()

def plot_tif(filename='test.tif'):
    """
    Read TIF file that has been previously generated
    """
    dataset = rio.open(filename)
    band1 = dataset.read(1,masked=True)
    show(band1,cmap='jet')
    msk=dataset.read_masks(1)

#################################################################
## Start doing some work

# Do we need varname as in input argument ?
# Yes we we do make it an enumeration of three possible things.
# Want to be able to simply inpout a URL
# urlinput='http://tds.renci.org:8080/thredds//dodsC/2020/nam/2020042912/hsofs/hatteras.renci.org/ncfs-dev-hsofs-nam-master/namforecast/maxele.63.nc'
#

# varname = 'zeta_max'
# varname = 'maxvel.63.nc'
# varname  'maxinundepth.63.nc'

def main(args):
    """
    Prototype script to construct geotif files ferm the ADCIRC trangulation grid
    """

    utilities.log.info(args)
    experimentTag = args.experiment_name
    filename = args.filename
    varname = args.varname
    showInterpolatedPlot = args.showInterpolatedPlot
    showRasterizedPlot = args.showRasterizedPlot

    utilities.log.info('Start ADSVIZ')
    config = utilities.load_config()

    url, dstr, cyc = construct_url(varname)
    if not validate_url(url):
        utilities.log.info('URL is invalid {}'.format(url))

    # Now construct filename destination using the dstr,cyc data

    iometadata = '_'.join([dstr,cyc])
    if experimentTag==None:
        rootdir = utilities.fetchBasedir(config['DEFAULT']['RDIR'], basedirExtra='APSVIZ_'+iometadata)
    else:
        rootdir = utilities.fetchBasedir(config['DEFAULT']['RDIR'], basedirExtra='APSVIZ_'+experimentTag+'_'+iometadata)
    filename = '/'.join([rootdir,filename])
    utilities.log.info('Using outputfilename of {}'.format(filename))

    # Build final pieces for the subsequent plots
    t0 = time.time()
    xtemp, ytemp, tri, targetgrid, targetepsg, advardict = construct_grid(url, varname)
    utilities.log.info('Building ADCIRC grid took {} secs'.format(time.time()-t0))

    ## everything above this line is a once-per-grid cost.  Only need the grid parts, not the actual solution
    ## Start the raster processing and build some plots.

    t0 = time.time()
    interp_lin = Tri.LinearTriInterpolator(tri,advardict['data'])
    utilities.log.info('Finished linearInterpolator in {} secs'.format(time.time()-t0))

    t0 = time.time()
    meshdict = compute_mesh(targetgrid, targetepsg)
    utilities.log.info('compute_mesh took {} secs'.format(time.time()-t0))
    xxm,yym = meshdict['xxm'], meshdict['yym']

    # zi_lin is the thing that will get rasterized down below
    zi_lin = interp_lin(xxm, yym)
    print(zi_lin.shape)
    utilities.log.debug('zi_lin {}'.format(zi_lin))

    t0 = time.time()
    write_tif(meshdict, zi_lin, targetgrid, targetepsg, filename)
    utilities.log.info('compute_mesh took {} secs'.format(time.time()-t0))

    if (showInterpolatedPlot):
        plot_triangular_vs_interpolated(meshdict, varname, tri, zi_lin, advardict)

    if (showRasterizedPlot):
        plot_tif(filename)

    utilities.log.info('Finished') 

if __name__ == '__main__':
    from argparse import ArgumentParser
    import sys
    parser = ArgumentParser()
    parser.add_argument('--experiment_name', action='store', dest='experiment_name', default=None,
                        help='Highlevel Experiment-tag value')
    parser.add_argument('--tif_filename', action='store', dest='filename', default='test.tif',
                        help='String: tif output file name will be prepended by new path')
    parser.add_argument('--showInterpolatedPlot', type=str2bool, action='store', dest='showInterpolatedPlot', default=True,
                        help='Boolean: Display the comparison of Trangular and interpolated plots')
    parser.add_argument('--showRasterizedPlot', type=str2bool, action='store', dest='showRasterizedPlot', default=True,
                        help='Boolean: Display the generated and saved tif plot')
    parser.add_argument('--varname', action='store', dest='varname', default='zeta_max',
                        help='String: zeta_max, vel_max, or inun_max')
    args = parser.parse_args()
    sys.exit(main(args))
