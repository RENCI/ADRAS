#!/usr/bin/env python

# SPDX-FileCopyrightText: 2022 Renaissance Computing Institute. All rights reserved.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-License-Identifier: LicenseRef-RENCI
# SPDX-License-Identifier: MIT

import numpy as np
import netCDF4 
import re

#class Utilities:
#    """
#
#    """
#    def __init__(self):
#        pass

def get_adcirc_grid(nc):
    agdict = {}
    agdict['lon'] = nc.variables['x'][:]
    agdict['lat'] = nc.variables['y'][:]
    agdict['ele'] = nc.variables['element'][:, :] - 1
    agdict['latmin'] = np.mean(nc.variables['y'][:])  # needed for scaling lon/lat plots
    agdict['depth'] = nc.variables['depth'][:]
    return agdict

def get_adcirc_slice(nc, v, it=None):
    advardict = {}
    var = nc.variables[v]
#    if re.search('max', v):
    if len(var.dimensions) == 1 :
        var_d = var[:]  # the actual data
    else:
        var_d = var[it, :]  # the actual data
    var_d[var_d.mask] = np.nan
    advardict['data'] = var_d
    return advardict

# get ADCIRC grid parts;  this need only be done once, as it can be time-consuming over the network
def extract_url_grid(url):
    """
    Extract ADCIRC grid parts from THREDDS dataset
    Results:
    """
    nc = netCDF4.Dataset(url)
    agdict = get_adcirc_grid(nc)
    return nc, agdict
