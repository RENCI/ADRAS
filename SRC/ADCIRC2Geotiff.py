#!/usr/bin/env python

# SPDX-FileCopyrightText: 2022 Renaissance Computing Institute. All rights reserved.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-License-Identifier: LicenseRef-RENCI
# SPDX-License-Identifier: MIT

# coding: utf-8

# If inputting url as a json then we must have filename metadata included as:
#{"1588269600000": "http://tds.renci.org:8080/thredds//dodsC/2020/nam/2020043018/hsofs/hatteras.renci.org/ncfs-dev-hsofs-nam-master/nowcast/maxele.63.nc"}
# If simply inputting a raw url (not expected to be a common approach) then the utime value will be set to 0000

import os
import sys
import re
import time
import datetime as dt
#import numpy.ma as ma
import pandas as pd
import numpy as np
import matplotlib.tri as Tri
import matplotlib.pyplot as plt
import matplotlib as mpl

import netCDF4

import rasterio as rio
from rasterio.transform import from_origin
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio.warp import calculate_default_transform, reproject

import geopandas as gpd
from affine import Affine

from PIL import Image
from colour import Color

import cartopy.crs as ccrs
import contextily as cx

#from utilities.adcirc_utilities import Utilities as adcirc_utilities
import utilities.adcirc_utilities as adcirc_utilities
from utilities.utilities import utilities as utilities
import logging
from utilities.logging import LoggingUtil 

global logger

#logger.info('Begin TIF generation')
#logger.info('netCDF4 Version = {}'.format(netCDF4.__version__))
#logger.info('Pandas Version = {}'.format(pd.__version__))
#logger.info('rasterio Version = {}'.format(rio.__version__))
#logger.info('Geopandas Version = {}'.format(gpd.__version__))

# define url functionality 
# http://tds.renci.org:8080/thredds/dodsC/2020/nam/2020012706/hsofs/hatteras.renci.org/ncfs-dev-hsofs-nam-master/namforecast/maxele.63.nc

allowable_vars = ['depth', 'zeta_max', 'vel_max', 'inun_max', 'wind_max']

def checkInputVar(v):
    if v.lower() in allowable_vars: return True
    return False

def get_interpolation_target(gridname=None, yamlfile=os.path.join(os.path.dirname(__file__), '..', 'config', 'main.yml')):
    """
    get geotiff grid parameters from the yaml
    Returns:
        targetgrid dict
        target crs 
    """
    try:
        config = utilities.load_config(yamlfile)['REGRID']
    except:
        logger.error("REGRID: load_config yaml failed")

    if not gridname:
        gridname = 'DEFAULT'
    if gridname not in config.keys():
        gridname = 'DEFAULT'

    targetgrid = {'Center_Latitude':  [config[gridname]['center_la']],
                  'Center_Longitude': [config[gridname]['center_lo']],
                  'res':        config[gridname]['res'],
                  'nx':         config[gridname]['nx'],
                  'ny':         config[gridname]['ny'],
                  'theta':      config[gridname]['theta']}

    return targetgrid, config[gridname]['adcirc_crs'], config[gridname]['target_crs']

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
    #return inun
    # return inundation depth masked out over "open water" (>1)
    return np.where(agdict['depth'] > 1, np.nan, inun)

def construct_geopandas(agdict, targetepsg):
    """
    Define geopandas processors
    project grid coords, before making Triangulation object
    """
    logger.info('Computing GeoPandas DF from ADCIRC grid')
    df_Adcirc = pd.DataFrame(
        {'Latitude': agdict['lat'],
        'Longitude': agdict['lon']})

    t0 = time.time()
    gdf = gpd.GeoDataFrame(
        df_Adcirc, geometry=gpd.points_from_xy(agdict['lon'], agdict['lat']))

    # init crs is LonLat, WGS84
    adcircepsg = agdict['crs']
    logger.info('Adding {} crs to initial GDF'.format(adcircepsg))
    gdf.crs = {'init': adcircepsg}
    logger.info('Converting GDF from {} to {}'.format(adcircepsg,targetepsg))
    gdf = gdf.to_crs({'init': targetepsg})
    logger.info('Time to create GDF was {}'.format(time.time()-t0))
    return gdf

def compute_geotiff_grid(targetgrid, adcircepsg, targetepsg):
    """
    project raster grid to target crs

    Results:
        rasdict. Dict of raster grid parameters and coords
        Values for upperleft_x, upperleft_y, x,y,xx,yy,xxm,yym
    """
    df_target = pd.DataFrame(data=targetgrid)
    gdf_target = gpd.GeoDataFrame(
        df_target, geometry=gpd.points_from_xy(df_target.Center_Longitude, df_target.Center_Latitude))

    # init projection is LonLat, WGS84
    gdf_target.crs = {'init': adcircepsg}

    # convert to "targetepsg"
    logger.info('Converting GDF from {} to {}'.format(adcircepsg,targetepsg))
    gdf_target = gdf_target.to_crs({'init': targetepsg})

    # compute spatial grid for raster 
    center_x = gdf_target['geometry'][0].x
    center_y = gdf_target['geometry'][0].y

    # length in x,y
    lxo2=targetgrid['nx']*targetgrid['res']/2
    lyo2=targetgrid['ny']*targetgrid['res']/2
    logger.info(f'lxo2={lxo2}, lyo2={lyo2}')
    logger.info(f'cenx={center_x}, ceny={center_y}')

    # compute mesh of cell centers
    x = np.arange(-np.floor(targetgrid['nx'] / 2),  np.floor(targetgrid['nx'] / 2) + 1)     * targetgrid['res'] + center_x
    y = np.arange( np.floor(targetgrid['ny'] / 2), -np.floor(targetgrid['ny'] / 2) - 1, -1) * targetgrid['res'] + center_y
    xx, yy = np.meshgrid(x, y)

    # move to origin and columnate
    xx0=(xx-center_x).ravel()
    yy0=(yy-center_y).ravel()

    # apply rotation
    ang=targetgrid['theta']*np.pi/180
    logger.debug(f"rotation ang={ang} radians")
    r = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
    Z = np.array([xx0, yy0])
    Zr=r.dot(Z)

    # translate back to origin
    xxmr=Zr[0,:].reshape(xx.shape)+center_x
    yymr=Zr[1,:].reshape(yy.shape)+center_y

    return {'uplx': xxmr[0,0],
            'uply': yymr[0,0],
            'cenx': center_x,
            'ceny': center_y,
            'x':    x,
            'y':    y,
            'xx':   xx,
            'yy':   yy,
            'xxm':  xxmr,
            'yym':  yymr,
            'nx':   x.shape,
            'ny':   y.shape}


def write_tif(rasdict, zi_lin, targetgrid, targetepsg, filename='test.tif'):
    """
    Construct the new TIF file and store it to disk in filename
    """
    xm0, ym0 = rasdict['uplx'], rasdict['uply']

    # specify transform
    a=targetgrid['res']
    b=0
    c=xm0 - a/2
    d=0
    e=-targetgrid['res']
    f=ym0 + a/2
    aft=Affine(a,b,c,b,e,f)*Affine.rotation(-targetgrid['theta'])
    logger.debug(f"TIF transform {aft}")

    md = {'crs':       targetepsg,
          'driver':    'GTiff',
          'height':    zi_lin.shape[0],
          'width':     zi_lin.shape[1],
          'count':     1,
          'dtype':     zi_lin.dtype,
          'nodata':    -99999,
          'transform': aft}
          
    # output a geo-referenced tiff
    dst = rio.open(filename, 'w', **md)
    try:
        dst.write(zi_lin, 1)
        logger.info(f"Wrote TIF file to {filename}")
    except:
        logger.error(f"Failed to write TIF file to {filename}")
    dst.close()

def getGridName(nc):
    """
    Return a (hopefully) unique grid name, based on the content of
    the netCDF file global attribute agrid, stripping out all spec chars
    except _
    """
    return re.sub(r'[^A-Za-z0-9_]+', '', getattr(nc,'agrid'))

def utm2WGS84(infile):
    #Define new projection
    new_crs = CRS.from_epsg(4326)  # WGS84

    #Read UTM data and reproject it to new projection
    with rio.open(infile) as dataset:
        new_transform, width, height = calculate_default_transform(
            dataset.crs, new_crs, dataset.width, dataset.height, *dataset.bounds)
        profile = dataset.profile.copy()
        profile.update(
            crs=new_crs,
            transform=new_transform,
            nodata=np.nan,
            height=height,
            width=width
        )
        imgdata = np.array([dataset.read(i) for i in dataset.indexes])
        new_imgdata = np.zeros((1,height,width), np.float64)
    
        reproject(source=imgdata,
                  destination=new_imgdata,
                  src_transform=dataset.transform,
                  src_crs=dataset.crs,
                  dst_transform=new_transform,
                  dst_crs=new_crs,
                  dst_nodata=np.nan,
                  resampling=Resampling.nearest)
    
    #Write reprojected data to tiff file
    output_file = infile[:-4]+'llr.tiff'
    with rio.open(output_file, 'w', **profile) as dst:
        dst.write(new_imgdata)

def get_discrete_cmap(vmin, vmax, bar):
    bottomvalue = vmin
    topvalue =  vmax
    bottomcolor = Color('#0000ff')
    topcolor = Color('#ff0000')
    colorramp=list(bottomcolor.range_to(topcolor, 32))

    hexlist = []
    for color in colorramp:
        hexlist.append(color.hex_l)

    if bar == False:
        hexlist[0] = '#00000000'
    elif bar == True:
        hexlist = hexlist
    else:
        print('You need to specify if color map is for a color bar True, or not False.')

    cmp = mpl.colors.ListedColormap(hexlist)

    return cmp

def create_colorbar(cmap,minval,maxval,barPathFile):
    """Create tick marks for values in meters"""
    valrange = abs(minval - maxval)

    ticks = [minval, valrange/4, valrange/2, valrange/1.33, maxval]

    tick1m = '<'+str("{:.2f}".format(ticks[0]))
    tick2m = str("{:.2f}".format(ticks[1]))
    tick3m = str("{:.2f}".format(ticks[2]))
    tick4m = str("{:.2f}".format(ticks[3]))
    tick5m = str("{:.2f}".format(ticks[4]))+'>'

    ticks_labels = [tick1m,tick2m,tick3m,tick4m,tick5m]

    """Get color map and plot range"""
    cmap = plt.cm.get_cmap(cmap)
    norm = mpl.colors.Normalize(vmin=minval, vmax=maxval)

    """Plot color bar and axis"""
    unit = 'm'
    fig, ax = plt.subplots(figsize=(1, 8))
    cbar = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, ticks=ticks, orientation='vertical')
    ax.tick_params(direction='out', length=10, width=2, labelsize=17, colors='black', grid_color='black', grid_alpha=1.0)
    cbar.set_label(unit, rotation=0, fontsize=17)
    cbar.ax.set_yticklabels(ticks_labels, va="center")
    fig.patch.set_facecolor('xkcd:pale grey')
    
    """Save colorbar image and close plot"""
    fig.savefig(barPathFile, bbox_inches = 'tight', pad_inches = 0.25)
    plt.close()

def scale_bar(ax, length=None, location=(0.5, 0.05), linewidth=3):
    """
    ax is the axes to draw the scalebar on.
    length is the length of the scalebar in km.
    location is center of the scalebar in axis coordinates.
    (ie. 0.5 is the middle of the plot)
    linewidth is the thickness of the scalebar.
    """
    #Get the limits of the axis in lat long
    llx0, llx1, lly0, lly1 = ax.get_extent(ccrs.PlateCarree())

    #Make tmc horizontally centred on the middle of the map,
    #vertically at scale bar location
    sbllx = (llx1 + llx0) / 2
    sblly = lly0 + (lly1 - lly0) * location[1]
    tmc = ccrs.TransverseMercator(sbllx, sblly, approx=True)

    #Get the extent of the plotted area in coordinates in metres
    x0, x1, y0, y1 = ax.get_extent(tmc)

    #Turn the specified scalebar location into coordinates in metres
    sbx = x0 + (x1 - x0) * location[0]
    sby = y0 + (y1 - y0) * location[1]

    #Calculate a scale bar length if none has been given
    #(Theres probably a more pythonic way of rounding the number but this works)
    if not length:
        length = (x1 - x0) / 5000 #in km
        ndim = int(np.floor(np.log10(length))) #number of digits in number
        length = round(length, -ndim) #round to 1sf
        #Returns numbers starting with the list
        def scale_number(x):
            if str(x)[0] in ['1', '2', '5']: return int(x)
            else: return scale_number(x - 10 ** ndim)
        length = scale_number(length)

    #Generate the x coordinate for the ends of the scalebar
    bar_xs = [sbx - length * 500, sbx + length * 500]
    #Plot the scalebar
    ax.plot(bar_xs, [sby, sby], transform=tmc, color='k', linewidth=linewidth)
    #Plot the scalebar label
    ax.text(sbx, sby, str(length) + ' km', transform=tmc,
            horizontalalignment='center', verticalalignment='bottom')

def get_concat_h_cut(infile, im1, im2):
    newsize = (im2.width, im1.height)
    im2 = im2.resize(newsize)
    dst = Image.new('RGB', (im1.width + im2.width, min(im1.height, im2.height)))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (im1.width, 0))
    dst.save(infile)

def createColorMap(infile,vmin,vmax):
    #Convert data in UTM projection to EPSG:4326 projection
    utm2WGS84(infile)

    #Create color map and use to create color bar
    cmap = get_discrete_cmap(vmin, vmax, True)
    create_colorbar(cmap,vmin,vmax,infile[:-4]+'colorbar.png')

    #Read data
    ds = rio.open(infile[:-4]+'llr.tiff')

    #Define map extent, center longitude, and projection
    extent = [ds.bounds.left, ds.bounds.right, ds.bounds.bottom, ds.bounds.top]
    cm_lon = (ds.transform * (ds.width/2, ds.height/2))[0]
    tran = ccrs.PlateCarree(central_longitude = cm_lon)
    proj = tran

    #Create color map for plot figure
    cmap = get_discrete_cmap(vmin, vmax, False)

    #Define dimensions and projection of plot image
    ysize = 23.4
    xsize = ysize*(ds.height/ds.width)
    fig = plt.figure(figsize=(ysize,xsize))
    ax = plt.axes(projection=proj)

    try:
        #Get OSM image, and convert it to EPSG 4326
        w, s, e, n = (ds.bounds.left,ds.bounds.bottom,ds.bounds.right,ds.bounds.top)
        cximg, cxext = cx.bounds2img(w, s, e, n, ll=True, source=cx.providers.OpenStreetMap.Mapnik)
        wimg, wext = cx.warp_tiles(cximg, cxext, "EPSG:4326")

        #Create map image from OSM warped data
        ax.imshow(wimg, extent=wext)

        #Create map image, with color and transparent background, from raw tiff in epsg 4326 projection
        ax.imshow(ds.read(1), extent=extent, transform=tran, interpolation='nearest', cmap=cmap, vmin=vmin, vmax=vmax)

    except:
        #Create map image, with color and transparent background, from raw tiff in epsg 4326 projection
        ax.imshow(ds.read(1), extent=extent, transform=tran, interpolation='nearest', cmap=cmap, vmin=vmin, vmax=vmax)

    #Add gridlines and set alpha to 0
    ax.gridlines(crs=ccrs.PlateCarree(central_longitude = cm_lon), draw_labels=True, linewidth=2,
                 color='gray', alpha=0.5, linestyle='--')
    ax.patch.set_alpha(0)

    #Add scale bar and north arrow
    scale_bar(ax, length=100, location=(0.8, 0.05), linewidth=3)
    x, y, arrow_length = 0.8, 0.2, 0.1
    ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
                arrowprops=dict(facecolor='black', width=5, headwidth=15),
                ha='center', va='center', fontsize=20,
                xycoords=ax.transAxes)

    #Save figure to png file and close
    fig.savefig(infile[:-4]+'ll.png', bbox_inches='tight')
    plt.close(fig)

    #Add color bar to image and save and png file
    im1 = Image.open(infile[:-4]+'ll.png')
    im2 = Image.open(infile[:-4]+'colorbar.png')
    get_concat_h_cut(infile[:-4]+'png', im1, im2)

    #Remove temp files
    os.remove(infile[:-4]+'llr.tiff')
    os.remove(infile[:-4]+'colorbar.png')
    os.remove(infile[:-4]+'ll.png')

#################################################################

def main(args):
    """
    Entry point for computing geotiff files from an ADCIRC netCDF file
    """

    log_level: int = int(os.getenv('LOG_LEVEL', logging.DEBUG))
    log_path: str = os.getenv('LOG_PATH', os.path.join(os.path.dirname(__file__), 'logs'))
    if not os.path.exists(log_path): os.mkdir(log_path)

    global logger
    logger = LoggingUtil.init_logging("APSVIZ.adras.hazus", level=log_level, line_format='long', log_file_path=log_path)
    logger.info(f"Start ADCIRC2Geotiff")

    varname = args.varname
    logger.info(f"varname={varname}")

    writePNG = False
    filename = '.'.join([varname,'tiff']) if args.filename is None else args.filename
    png_filename = '.'.join([varname,'png']) if args.png_filename is None else args.png_filename
#    filename = args.datadir + '/' + filename
#    png_filename = args.datadir + '/' + png_filename

    if not checkInputVar(varname):
        logger.error(f"Variable {varname} not yet supported.")
        sys.exit(1)

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
        logger.debug(s3_resource)
        logger.debug(s3_utilities.config)

        thisBucket = s3_utilities.config['S3_UPLOAD_Main_Bucket']
        thisRegion = s3_utilities.config['region_name']

        logger.info(f"Bucket={thisBucket}")
        logger.info(f"Region={thisRegion}")
        
        if not s3_utilities.bucket_exists(thisBucket):
            res = s3_utilities.create_bucket(thisRegion, thisBucket)
            logger.info(f"Bucket {thisBucket} created.")
        else:
            logger.info(f"Bucket {thisBucket} already exists.")

    url = args.url

    rootdir = utilities.getBasedir(main_config['DEFAULT']['RDIR'])

    targetgrid, adcircepsg, targetepsg = get_interpolation_target(gridname=args.gridname, yamlfile=raster_yaml_file)

    logger.info(f"Extracting ADCIRC grid from nc object...")
    t0 = time.time()
    nc, agdict = adcirc_utilities.extract_url_grid(url)
    agdict['crs'] = adcircepsg

    # get grid name for building geopandas (gdf) filename
    gridname = args.gridname
    if not args.gridname:
        gridname = getGridName(nc)
    logger.info(f"ADCIRC grid name is {gridname}")

    # Construct or load existing geopandas object on the input URL grid
    t0 = time.time()
    gdf_pklfile = '{}.gdf.pkl'.format(gridname)

    # Need to check in specified dirs
    f = os.path.join(args.pkldir, gdf_pklfile)
    gdf = construct_geopandas(agdict, targetepsg)

    if not os.path.exists(f):
        gdf = construct_geopandas(agdict, targetepsg)
        if not os.path.exists(args.pkldir): os.makedirs(args.pkldir)
        utilities.writePickle(gdf, filename=f)
        logger.info(f"Wrote Geopandas file to {f}")
    else:
        logger.info(f"{f} exists. Using it...")
        gdf = pd.read_pickle(f)

    # Extract the lat,lon values of the current gdf object
    logger.info(f"Extracting X and Y from the ADCIRC grid GDF ..")
    xtemp, ytemp = gdf['geometry'].x, gdf['geometry'].y
    logger.debug(f"Extracted X and Y from the ADCIRC grid GDF")

    # Build Triangulate object for interpolating the input geopandas object
    logger.info(f"Computing X and Y triangulation object ...")
    t0 = time.time()
    tri = Tri.Triangulation(xtemp, ytemp, triangles=agdict['ele'])
    logger.debug(f"Tri interpolation object took {time.time()-t0:6.2f} secs to compute.")

    # Set up grid for geotiff data
    logger.debug(f"Execing compute_geotiff_grid ...")
    t0 = time.time()
    rasdict = compute_geotiff_grid(targetgrid, adcircepsg, targetepsg)
    logger.debug(f"compute_geotiff_grid took {time.time()-t0:6.2f} secs")
    
    # extract the raster pixel points
    xxm, yym = rasdict['xxm'], rasdict['yym']

    orig_filename = filename
    orig_png_filename = png_filename

    # filename = os.path.join(rootdir, orig_filename)
    filename = orig_filename
    png_filename = os.path.join(rootdir, orig_png_filename)

    nc = netCDF4.Dataset(url)
    if varname == 'inun_max':
        advardict = adcirc_utilities.get_adcirc_slice(nc, 'zeta_max')
        # compute inundation and replace advardict['data']
        advardict['data'] = computeInundation(advardict, agdict)
    else:
        advardict = adcirc_utilities.get_adcirc_slice(nc, varname)

    vmin = np.nanmin(advardict['data'])
    vmax = np.nanmax(advardict['data'])
    logger.info(f"Min/Max in ADCIRC Slice: {vmin:6.2f}, {vmax:6.2f}")

    # construct the interpolator
    logger.info(f"Computing linearInterpolator for variable ...")
    t0 = time.time()
    interp_lin = Tri.LinearTriInterpolator(tri, advardict['data'])
    logger.debug(f"Finished linearInterpolator in {time.time()-t0:6.2f} secs")

    # interpolate data to raster pixels
    zi_lin = interp_lin(xxm, yym)
    temp = pd.DataFrame({'lon': xxm.ravel(), 'lat': yym.ravel(), 'val': zi_lin.ravel(), })
    temp.to_csv('test.csv')

    #zi_lin = np.where(np.isnan(zi_lin) , 0, zi_lin)
    ivmin = np.nanmin(zi_lin)
    ivmax = np.nanmax(zi_lin)
    logger.info(f"Min/Max in interpolated grid: {ivmin:6.2f}, {ivmax:6.2f}")

    logger.info(f"Outputting tiff file {filename}")
    write_tif(rasdict, zi_lin, targetgrid, targetepsg,  args.datadir + '/' + filename)

    llpngfilename = filename.replace('tiff','png')
    logger.info(f"Create color png map {llpngfilename}")
    #createColorMap(args.datadir + '/' + filename,ivmin,ivmax)

    if main_config['S3']['SEND2AWS']:
        os.chdir(args.datadir)
        resp = s3_utilities.upload(thisBucket, args.s3path, filename)
        msg=f"Upload to s3://{thisBucket}/{args.s3path}/{args.filename} "
        if not resp:
            logger.info(f"{msg} failed.")
        else:
            logger.info(f"{msg} succeeded.")

        if os.path.exists(llpngfilename):
            resp = s3_utilities.upload(thisBucket, args.s3path, llpngfilename)
            msg=f"Upload to s3://{thisBucket}/{args.s3path}/{llpngfilename} "
            if not resp:
                logger.info(f"{msg} failed.")
            else:
                logger.info(f"{msg} succeeded.")

    #if png_filename is not None:
    #    logger.info('Outputting png file {}'.format(png_filename))
    #    write_png(filename, png_filename)

    logger.info("Finished")

if __name__ == '__main__':

    from argparse import ArgumentParser
    import sys
    parser = ArgumentParser()

    parser.add_argument('--datadir', action='store', dest='datadir', default='./',
                        help='String: Directory for image output file. ')
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
