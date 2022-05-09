#!/usr/bin/env bash

rasterParameters () {
  case $gridname in 
  "nc_inundation_v9.99_w_rivers" | "ncv999wr" | "nc_inundation_v9.99_rivers" | "ncv9.99d" | "nc_inundation_v9.99d" )
    upperleft_lo=-77.625
    upperleft_la=36.75
    res=50  # resolution in m
    nx=1700
    ny=9500
    theta=-54
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

#  "hsofs" ) 
#    upperleft_lo=-74.25
#    upperleft_la=41.35
#    # resolution in m
#    res=100  
#    nx=2100
#    ny=900
#    theta=15
#    target_crs="epsg:32619"
#    adcirc_crs="epsg:4326"
#  ;;

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
    upperleft_lo=-72.7
    upperleft_la=41.75
    # resolution in m
    res=500  
    nx=400
    ny=180
    theta=15
    target_crs="epsg:32619"
    adcirc_crs="epsg:4326"
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

  "EGOMv20b")
    upperleft_lo=-82.9
    upperleft_la=30.5
    res=50  # resolution in m
    nx=2000
    ny=10000
    theta=27.5
    target_crs="epsg:32614"
    adcirc_crs="epsg:4326"
  ;;

  "NGOMv19b")
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

