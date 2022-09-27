#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2022 Renaissance Computing Institute. All rights reserved.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-License-Identifier: LicenseRef-RENCI
# SPDX-License-Identifier: MIT

rasterParameters () {
  case $gridname in 
  "nc_inundation_v9.99_w_rivers" | "ncv999wr" | "nc_inundation_v9.99a_w_rivers" | "ncv9.99d" | "nc_inundation_v9.99d" )
    upperleft_lo=-78.25
    upperleft_la=37.5
    res=500 # resolution in m
    nx=300
    ny=1000
    theta=-55
    target_crs="epsg:32617"
    adcirc_crs="epsg:4326"
  ;;

  "NCSC_SAB_v1.23" | "ncsc123" )
    upperleft_lo=-82.45
    upperleft_la=33.5
    res=500  # resolution in m
    nx=1260
    ny=90
    theta=40
    target_crs="epsg:32617"
    adcirc_crs="epsg:4326"
  ;;

  "SABv20a" )
    upperleft_lo=-82.375
    upperleft_la=31.5
    res=50  # resolution in m
    nx=4200
    ny=720
    theta=65
    target_crs="epsg:32617"
    adcirc_crs="epsg:4326"
  ;;

  "hsofs" ) 
    # 2022, al09, GA coast
    upperleft_lo=-81.1
    upperleft_la=33.3     
    theta=-20
    res=50  # resolution in meters
    nx=2000
    ny=7000
    target_crs="epsg:32619"
    adcirc_crs="epsg:4326"
  ;;

  "uriv18" | "uriv18_cl" | "NAC2014" ) 
    upperleft_lo=-72.7
    upperleft_la=41.75
    # resolution in m
    res=50  
    nx=4000
    ny=1800
    theta=15
    target_crs="epsg:32619"
    adcirc_crs="epsg:4326"
  ;;

  "ec95d") 
    upperleft_lo=-85.25       
    upperleft_la=29.5
    theta=-50
    res=500  # resolution in meters
    nx=900
    ny=300
    #upperleft_lo=-85.2       
    #upperleft_la=29.939
    #theta=-40
    #res=50  # resolution in meters
    #nx=7000
    #ny=3000
    target_crs="epsg:32614"
    adcirc_crs="epsg:4326"

    #upperleft_lo=-82.9
    #upperleft_la=30.5
    #res=5000  # resolution in m
    #nx=20
    #ny=100
    #theta=27.5
    #target_crs="epsg:32614"
    #adcirc_crs="epsg:4326"
    #upperleft_lo=-72.7
    #upperleft_la=41.75
    ## resolution in m
    #res=500  
    #nx=400
    #ny=180
    #theta=15
    #target_crs="epsg:32619"
    #adcirc_crs="epsg:4326"
    #upperleft_lo=-99
    #upperleft_la=32
    #res=1000  # resolution in m
    #nx=2000
    #ny=2000
    #theta=0
    #target_crs="epsg:6346"
    #adcirc_crs="epsg:4326"
  ;;

  "LA_v20a-WithUpperAtch_chk" | "LAv20a" | "LAv21a" | "hsofs" )
    upperleft_lo=-94.25
    upperleft_la=30.5
    res=50  # resolution in m
    nx=10000
    ny=3900
    theta=0
    target_crs="epsg:32614"
    adcirc_crs="epsg:4326"
  ;;

  EGOM*)
    upperleft_lo=-86.25       
    upperleft_la=28.4
    theta=-60
    res=50  # resolution in meters
    nx=15000
    ny=3000
    #upperleft_lo=-85.25       
    #upperleft_la=29.5
    #theta=-50
    #res=50  # resolution in meters
    #nx=9000
    #ny=3000
    target_crs="epsg:32614"
    adcirc_crs="epsg:4326"
  ;;

  NGOM*)
    upperleft_lo=-91.5
    upperleft_la=31.0
    res=50  # resolution in m
    nx=10000
    ny=6000
    theta=0
    target_crs="epsg:32614"
    adcirc_crs="epsg:4326"
  ;;

  "NGOM_RT_v20a_chk")
    upperleft_lo=-85.8
    upperleft_la=30.75
    res=25  # resolution in m
    nx=10000
    ny=6000
	theta=0
    target_crs='epsg:32614'
    adcirc_crs='epsg:4326'
  ;;

  "TX2008")
    upperleft_lo=-98.7
    upperleft_la=29.5
    # resolution in m
    res=50  
    nx=11000
    ny=3000
    theta=35
    target_crs="epsg:32614"
    adcirc_crs="epsg:4326"
  ;;

  "PRVI15" | "prvi15" ) 
    upperleft_lo=-67.30
    upperleft_la=18.535
    res=5000  # resolution in m
    nx=36
    ny=15
    theta=0
    target_crs="epsg:3920"
    adcirc_crs="epsg:4326"
  ;;

  *)
    upperleft_lo=-77.09
    upperleft_la=35.7
    res=100  # resolution in m
    nx=1000
    ny=1000
    theta=0
    target_crs="epsg:6346"
    adcirc_crs="epsg:4326"
  ;;

  esac
}

