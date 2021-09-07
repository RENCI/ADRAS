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
from PIL import Image
from colour import Color

from utilities.utilities import utilities as utilities
import logging
from utilities.logging import LoggingUtil 

global logger

def get_discrete_cmap(valueList):
    # defind bottom and top data values 
    bottomvalue = 0.0
    topvalue =  5.0

    # create top and bottom color values
    bottomcolor = Color('#0000ff')
    topcolor = Color('#ff0000')

    # Create list of 32 discrete color values
    colorramp=list(bottomcolor.range_to(topcolor, 32))

    # Output list of color values as list of hex colors
    hexlist = []
    for color in colorramp:
        hexlist.append(color.hex_l)

    # Create color map from list of hex colors, and return from program
    cmp = mpl.colors.ListedColormap(hexlist)
    return cmp

def rotate_img(img_path, rt_degr):
    '''
    This function rotates the color bar image so it is horizontal
    '''
    img = Image.open(img_path)
    return img.rotate(rt_degr, expand=1)

def create_colorbar(cmap,values,barPathFile):
    """Create tick marks for values in meters"""
    valrange = abs(values[0] - values[-1])

    ticks = [values[0], valrange/4, valrange/2, valrange/1.33, values[-1]]

    tick1m = '<'+str("{:.2f}".format(ticks[0]))
    tick2m = str("{:.2f}".format(ticks[1]))
    tick3m = str("{:.2f}".format(ticks[2]))
    tick4m = str("{:.2f}".format(ticks[3]))
    tick5m = str("{:.2f}".format(ticks[4]))+'>'

    ticks_labels = [tick1m,tick2m,tick3m,tick4m,tick5m]

    """Get color map and plot range"""
    cmap = plt.cm.get_cmap(cmap)
    norm = mpl.colors.Normalize(vmin=values[0], vmax=values[-1])

    """Plot color bar and first axis"""
    unit = 'meters'
    fig, ax = plt.subplots(figsize=(1, 8))
    cbar = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, ticks=ticks, orientation='vertical')
    cbar.ax.yaxis.set_label_position("left")
    ax.tick_params(direction='out', length=10, width=2, labelsize=17, colors='black', grid_color='black', grid_alpha=0.5)
    cbar.set_label(unit, fontsize=17)
    cbar.ax.set_yticklabels(ticks_labels, rotation=90, va="center")

    """Create tick marks for values in feet"""
    econversionval = 3.28084
    eunit = 'feet'

    valrangeft = valrange * econversionval
    iticks = [(values[0] * econversionval), valrangeft/4, valrangeft/2, valrangeft/1.33, (values[-1] * econversionval)]

    tick1ft = '<'+str("{:.2f}".format(iticks[0]))
    tick2ft = str("{:.2f}".format(iticks[1]))
    tick3ft = str("{:.2f}".format(iticks[2]))
    tick4ft = str("{:.2f}".format(iticks[3]))
    tick5ft = str("{:.2f}".format(iticks[4]))+'>'

    iticks_labels = [tick1ft,tick2ft,tick3ft,tick4ft,tick5ft]

    """Plot second axis"""
    ax2 = ax.twinx()
    ax2.tick_params(direction='out', length=10, width=2, labelsize=17, colors='black', grid_color='black', grid_alpha=0.5)
    ax2.set_ylim([(values[0] * econversionval),(values[-1] * econversionval)])
    ax2.set_yticks(iticks)
    ax2.set_yticklabels(iticks_labels, rotation=90, va="center")
    ax2.set_ylabel(eunit, fontsize=17)

    """Save colorbar image and close plot"""
    fig.savefig(barPathFile, transparent=True, bbox_inches = 'tight', pad_inches = 0.25)
    plt.close()

    """Rotate colorbar so it is horizontal"""
    img_rt_270 = rotate_img(barPathFile, 270)
    img_rt_270.save(barPathFile)

def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return(tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3)))

def getColors():
    # define bottom and top data values
    bottomvalue = 0.0
    topvalue =  5.0

    # create top and bottom color values 
    bottomcolor = Color('#0000ff')
    topcolor = Color('#ff0000')

    # create list of 32 discrete color values
    colorramp=list(bottomcolor.range_to(topcolor, 32))

    # calculate range value between the bottom and top color values
    if bottomvalue < 0:
        vrange = topvalue + bottomvalue
    else:
        vrange = topvalue - bottomvalue 

    # create list of data values
    valueLst = np.append(np.arange(bottomvalue, topvalue, vrange/30), topvalue)

    # create list of data values and color hex values 
    valcolorlst = []    
    for i in range(len(valueLst)-1):
        valcolorlst.append([valueLst[i], valueLst[i+1], hex_to_rgb(colorramp[i+1].hex_l)])

    # Return data value list and data value/color value lists
    return(valueLst, valcolorlst)

def styleRaster(input_file, valcolorlst):
    # Open raw geotiff file, get profile information, and read data in as array
    with rio.open(input_file) as raster:
        profile = raster.profile.copy()
        profile.update(
            dtype=rio.uint8,
            count=4,
            nodata=None,
            compress='lzw')
        data = raster.read()
        data = np.reshape(data, (data.shape[1], data.shape[2]))

    # Create empy band arrays for input of color values
    outband1 = np.empty((data.shape[0], data.shape[1]), dtype='int')
    outband2 = np.empty((data.shape[0], data.shape[1]), dtype='int')
    outband3 = np.empty((data.shape[0], data.shape[1]), dtype='int')
    outband4 = np.empty((data.shape[0], data.shape[1]), dtype='int')

    # Define bottom color value, create index from data using value from valcolorlst, and input into outbands
    bottomcolor = hex_to_rgb(Color('#0000ff').hex_l)
    index = np.where(data <= valcolorlst[0][0])
    outband1[index] = bottomcolor[0]
    outband2[index] = bottomcolor[1]
    outband3[index] = bottomcolor[2]
    outband4[index] = 191

    # loop through valcolorlist, creating index from data, and input into outbands
    for i in range(len(valcolorlst)):
        index = np.where((data > valcolorlst[i][0]) & (data <= valcolorlst[i][1]))
        outband1[index] = valcolorlst[i][2][0]
        outband2[index] = valcolorlst[i][2][1]
        outband3[index] = valcolorlst[i][2][2]
        outband4[index] = 191

    # Define top color value, create index from data using value from valcolorlist, and input into outbands
    topcolor = hex_to_rgb(Color('#ff0000').hex_l)
    index = np.where(data > valcolorlst[-1][1])
    outband1[index] = topcolor[0]
    outband2[index] = topcolor[1]
    outband3[index] = topcolor[2]
    outband4[index] = 191

    # mergre outbands into arr_merge array
    arr_merge = np.array([outband1, outband2, outband3, outband4])

    # define output file name
    output_file = input_file[:-5]+'.color.tiff'

    # write arr_merg to tiff file using information from profile
    with rio.open(output_file, 'w', **profile) as dst:
        # merge happens in sequence of rasters
        dst.write(arr_merge.astype(rio.uint8))

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

    # Create list of 30 discrete colors and values
    valueList, valColorLst = getColors()

    # Create color tiff
    styleRaster(infilename, valColorLst)

    # Get color map for colorbar
    cmap = get_discrete_cmap(valueList)

    # Get colorbar file path
    barPathFile = infilename[:-5]+'.colorbar.png'

    # Create colorbar
    create_colorbar(cmap,valueList,barPathFile)

    logger.info("Finished")

if __name__ == '__main__':

    from argparse import ArgumentParser
    import sys
    parser = ArgumentParser()

    parser.add_argument('--tif_infilename', action='store', dest='infilename', default=None,
                        help='String: tiff input file name.')

    args = parser.parse_args()
    sys.exit(main(args))
