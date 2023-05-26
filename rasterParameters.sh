#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2022 Renaissance Computing Institute. All rights reserved.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-License-Identifier: LicenseRef-RENCI
# SPDX-License-Identifier: MIT

rasterParameters () {
  case $gridname in 
  "nc_inundation_v9.99_w_rivers" | "ncv999wr" | "nc_inundation_v9.99a_w_rivers" | "ncv9.99d" | "nc_inundation_v9.99d" )
    center_lo=-78.25
    center_la=37.5
    res=500 # resolution in m
    nx=300
    ny=1000
    theta=-55
    target_crs="epsg:32617"
    adcirc_crs="epsg:4326"
  ;;

  NCSC* | ncsc* )
    #center_lo=-82.45
    #center_la=33.5
    #res=500  # resolution in m
    #nx=1260
    #ny=90
    #theta=40
    center_lo=-81
    center_la=32
    res="1"  # resolution in m
    nx=10
    ny=10
    theta=0
    #target_crs="epsg:32617"
    target_crs="epsg:4326"
    adcirc_crs="epsg:4326"
  ;;

  "SABv20a" )
    center_lo=-82.375
    center_la=31.5
    res=50  # resolution in m
    nx=4200
    ny=720
    theta=65
    target_crs="epsg:32617"
    adcirc_crs="epsg:4326"
  ;;

  "hsofs" ) 
    # 2022, al09, GA coast
    center_lo=-81.4
    center_la=34.3
    theta=-25
    res=100  # resolution in meters
    nx=1500
    ny=5600
    # 2022, al09, WFL coast
    #center_lo=-83.5
    #center_la=27.0
    #theta=-60
    #res=50  # resolution in meters
    #nx=6000
    #ny=1500
    target_crs="epsg:32617"
    adcirc_crs="epsg:4326"
  ;;

  "uriv18" | "uriv18_cl" | "NAC2014" ) 
    center_lo=-72.7
    center_la=41.75
    # resolution in m
    res=50  
    nx=4000
    ny=1800
    theta=15
    target_crs="epsg:32619"
    adcirc_crs="epsg:4326"
  ;;

  "ec95d") 
    center_lo=-75 
    center_la=35  
    res=1000  # resolution in m
    theta=45
    nx=100
    ny=100
    adcirc_crs="epsg:4326"
    target_crs="epsg:32617"

  ;;

  "LA_v20a-WithUpperAtch_chk" | "LAv20a" | "LAv21a" | "hsofs" )
    center_lo=-94.25
    center_la=30.5
    res=50  # resolution in m
    nx=10000
    ny=3900
    theta=0
    target_crs="epsg:32614"
    adcirc_crs="epsg:4326"
  ;;

  EGOM* | egom* )
    center_lo=-81.97
    center_la=26.467
    theta=-60
    res=50  # resolution in meters
    nx=4000
    ny=1500
    target_crs="epsg:32617"
    adcirc_crs="epsg:4326"
  ;;

  NGOM* | ngom* )
    center_lo=-91.5
    center_la=31.0
    res=50  # resolution in m
    nx=10000
    ny=6000
    theta=0
    target_crs="epsg:32616"
    adcirc_crs="epsg:4326"
  ;;

  "TX2008")
    center_lo=-98.7
    center_la=29.5
    # resolution in m
    res=50  
    nx=11000
    ny=3000
    theta=35
    target_crs="epsg:32614"
    adcirc_crs="epsg:4326"
  ;;

  "PRVI15" | "prvi15" ) 
    center_lo=-67.30
    center_la=18.535
    res=5000  # resolution in m
    nx=36
    ny=15
    theta=0
    target_crs="epsg:3920"
    adcirc_crs="epsg:4326"
  ;;

  *)
    center_lo=-77.09
    center_la=35.7
    res=100  # resolution in m
    nx=1000
    ny=1000
    theta=0
    target_crs="epsg:6346"
    adcirc_crs="epsg:4326"
  ;;

  esac
}

