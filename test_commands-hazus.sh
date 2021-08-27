#!/bin/bash
PYTHONPATH="${PYTHONPATH:-$HOME/GitHub/RENCI/ADRAS}"
export PYTHONPATH
echo $PYTHONPATH

#url='https://fortytwo.cct.lsu.edu:443/thredds/fileServer/2021/al03/06/LA_v20a-WithUpperAtch_chk/supermic.hpc.lsu.edu/LAv20a_al032021_jgf_23kcms/nhcConsensus/'
#url='http://tds.renci.org:8080/thredds/fileServer/2021/nam/2021060918/ec95d/hatteras.renci.org/ec95d-nam-bob-rptest/namforecast/'
#url='http://tds.renci.org:8080/thredds/fileServer/2021/nam/2021082006/uriv18/hatteras.renci.org/uriv18-nam-bob-2021/namforecast/'
#url='http://adcircvis.tacc.utexas.edu:8080/thredds/fileServer/asgs/2021/al09/01/LAv21a/frontera.tacc.utexas.edu/LAv21a_al092021_jgf_10kcms/nhcConsensus/'
url='http://adcircvis.tacc.utexas.edu:8080/thredds/fileServer/asgs/2021/al09/05/LAv20a/frontera.tacc.utexas.edu/LAv20a_al092021_jgf_10kcms/veerRight50/'
bash compute_geotiffs.sh  $url

