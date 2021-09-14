#!/usr/bin/env python
# coding: utf-8

# If inputing filename :
# If simply inputting a raw url (not expected to be a common approach) then the utime value will be set to 0000

import os
import sys

import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl

import rasterio as rio
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio.warp import calculate_default_transform, reproject

from PIL import Image
from colour import Color

import cartopy.crs as ccrs
import contextily as cx

from utilities.utilities import utilities as utilities
import logging
from utilities.logging import LoggingUtil 

global logger

def utm2WGS84(infile):
    #Define new projection
    new_crs = CRS.from_epsg(4326)  # WGS84

    #Read UTM data and reproject it to new projection
    with rio.open(infile) as dataset:
        new_transform, width, height = calculate_default_transform(
            dataset.crs, new_crs, 
            dataset.width, dataset.height, *dataset.bounds)
        profile = dataset.profile.copy()
        profile.update(
            crs=new_crs,
            transform=new_transform
        )
        imgdata = np.array([dataset.read(i) for i in dataset.indexes])
        new_imgdata = np.zeros(imgdata.shape)
    
        reproject(source=imgdata,
                  destination=new_imgdata,
                  src_transform=dataset.transform,
                  src_crs=dataset.crs,
                  dst_transform=new_transform,
                  dst_crs=new_crs,
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

    #print(hexlist)
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

    #Get OSM image, and convert it to EPSG 4326
    w, s, e, n = (ds.bounds.left,ds.bounds.bottom,ds.bounds.right,ds.bounds.top)
    cximg, cxext = cx.bounds2img(w, s, e, n, ll=True, source=cx.providers.OpenStreetMap.Mapnik)
    wimg, wext = cx.warp_tiles(cximg, cxext, "EPSG:4326")

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

    #Create map image from OSM warped data
    ax.imshow(wimg, extent=wext)

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
    Entry point for ceate color geotiff files from raw geotiff files output from ADCIRC2geotiff.py 
    """
    log_level: int = int(os.getenv('LOG_LEVEL', logging.DEBUG))
    log_path: str = os.getenv('LOG_PATH', os.path.join(os.path.dirname(__file__), 'logs'))
    if not os.path.exists(log_path): os.mkdir(log_path)

    global logger
    logger = LoggingUtil.init_logging("APSVIZ.adras.hazus", level=log_level, line_format='long', log_file_path=log_path)
    logger.info(f"Start rawGeotiff2colorGeotiff")

    infilename = args.infilename
    logger.info(f"infilename={infilename}")

    vmin = 0.0
    vmax = 5.0

    createColorMap(infilename,vmin,vmax)

    logger.info("Finished")

if __name__ == '__main__':

    from argparse import ArgumentParser
    import sys
    parser = ArgumentParser()

    parser.add_argument('--tif_infilename', action='store', dest='infilename', default=None,
                        help='String: tiff input file name.')

    args = parser.parse_args()
    sys.exit(main(args))
