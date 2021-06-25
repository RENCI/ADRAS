# ADRAS

Rasterizer for ADCIRC.  Primarily for computing inundation from maxele and the grid's topo/bathy

run as 

bash compute_geotiffs.sh <downloadurl>

Example: 

url='http://fortytwo.cct.lsu.edu/thredds/fileServer/2021/al03/06/LA_v20a-WithUpperAtch_chk/supermic.hpc.lsu.edu/LAv20a_al032021_jgf_23kcms/nhcConsensus'
bash compute_geotiffs.sh $url

local output will be 2 geotiff files, which will have been uploaded to the AWS S3 bucket hazus

