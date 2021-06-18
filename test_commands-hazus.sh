#!/bin/bash
source ./init-ADRAS.sh
echo $PYTHONPATH
echo $GDAL_DATA

adv=28
s3path="laura/$adv/"
fname="laura_adv$adv.tif"

url="http://fortytwo.cct.lsu.edu/thredds/dodsC/2020/laura/$adv/LA_v20a-WithUpperAtch_chk/qbc.loni.org/LAv20a_al132020_jgf/nhcConsensus/maxele.63.nc"
varname="inun_max"
gridname="LA_v20a-WithUpperAtch_chk"
echo $PYTHON

if [ ! -d "pklfiles" ]; then
        mkdir pklfiles
fi
ls $ADRASHOME/SRC/newADCIRC2Geotiff.py

$PYTHON $ADRASHOME/SRC/newADCIRC2Geotiff.py --varname=$varname --gridname=$gridname --url=$url --tif_filename=$fname --s3path=$s3path

