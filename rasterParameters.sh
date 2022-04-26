#!/usr/bin/env bash

rasterParameters () {
  case $gridname in 
  "nc_inundation_v9.99_w_rivers" | "ncv999wr" | "nc_inundation_v9.99_rivers" | "ncv9.99d" | "nc_inundation_v9.99d" )
    upperleft_lo=-78.75
    upperleft_la=36.75
    res=100  # resolution in m
    nx=3900
    ny=3500
	theta=0
    target_crs='epsg:32617'
    adcirc_crs='epsg:4326'
  ;;

  "ec95d") 
    upperleft_lo=-82.5
    upperleft_la=36.0
    res=1000  # resolution in m
    nx=500
    ny=500
	theta=0
    target_crs='epsg:32617'
    adcirc_crs='epsg:4326'
  ;;

  "LA_v20a-WithUpperAtch_chk")
    upperleft_lo=-95.0
    upperleft_la=31.0
    res=50  # resolution in m
    nx=10000
    ny=6000
	theta=0
    target_crs='epsg:32614'
    adcirc_crs='epsg:4326'
  ;;

  "EGOMv20b")
    upperleft_lo=-84
    upperleft_la=30.5
    res=50  # resolution in m
    nx=6000
    ny=10000
	theta=0
    target_crs='epsg:32614'
    adcirc_crs='epsg:4326'
  ;;

  "NGOMv19b")
    upperleft_lo=-91.5
    upperleft_la=31.0
    res=50  # resolution in m
    nx=10000
    ny=6000
	theta=0
    target_crs='epsg:32614'
    adcirc_crs='epsg:4326'
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
    upperleft_lo=-98.0
    upperleft_la=30.5
    # resolution in m
    res=50  
    nx=12000
    ny=12000
	theta=0
    target_crs='epsg:32614'
    adcirc_crs='epsg:4326'
  ;;

  *)
    upperleft_lo=-77.09
    upperleft_la=35.7
    res=100  # resolution in m
    nx=1000
    ny=1000
	theta=0
    target_crs='epsg:6346'
    adcirc_crs='epsg:4326'
  ;;

  esac
}

