#!/bin/bash
PYTHONPATH="${PYTHONPATH:-$HOME/GitHub/RENCI/ADRAS}"
export PYTHONPATH
echo $PYTHONPATH

#url='https://fortytwo.cct.lsu.edu:443/thredds/fileServer/2020/laura/28/LA_v20a-WithUpperAtch_chk/qbc.loni.org/LAv20a_al132020_jgf/nhcConsensus/'
url='https://fortytwo.cct.lsu.edu:443/thredds/fileServer/2021/al03/06/LA_v20a-WithUpperAtch_chk/supermic.hpc.lsu.edu/LAv20a_al032021_jgf_23kcms/nhcConsensus/'
#url='http://tds.renci.org:8080/thredds/fileServer/2021/nam/2021060918/ec95d/hatteras.renci.org/ec95d-nam-bob-rptest/namforecast/'
#url='http://tds.renci.org:8080/thredds/fileServer/2021/nam/2021082006/uriv18/hatteras.renci.org/uriv18-nam-bob-2021/namforecast/'
#url='http://adcircvis.tacc.utexas.edu:8080/thredds/fileServer/asgs/2021/al09/01/LAv21a/frontera.tacc.utexas.edu/LAv21a_al092021_jgf_10kcms/nhcConsensus/'
#url='http://adcircvis.tacc.utexas.edu:8080/thredds/fileServer/asgs/2021/al09/05/LAv20a/frontera.tacc.utexas.edu/LAv20a_al092021_jgf_10kcms/veerRight50/'
#url='http://adcircvis.tacc.utexas.edu:8080/thredds/fileServer/asgs/2021/al09/05/LAv20a/stampede2.tacc.utexas.edu/LAv20a_al092021_bde/nhcConsensus/'
#url="http://tds.renci.org:8080/thredds/fileServer/2022/nam/2022011006/hsofs/hatteras.renci.org/hsofs-nam-bob-2021/namforecast/"
url="http://tds.renci.org:8080/thredds/fileServer/2021/ida/12/ec95d/hatteras.renci.org/ec95d-al09-bob/nhcOfcl/"
#url='http://tds.renci.org:8080/thredds/fileServer/2021/nam/2021092312/NCSC_SAB_v1.15/hatteras.renci.org/ncsc115-nam-2021/namforecast/'
#url='http://adcircvis.tacc.utexas.edu:8080/thredds/fileServer/asgs/2021/nam/2021092306/SABv20a/frontera.tacc.utexas.edu/SABv20a_nam_jgf_status/namforecast/'
#url='http://adcircvis.tacc.utexas.edu:8080/thredds/fileServer/asgs/2021/nam/2021092906/SABv20a/frontera.tacc.utexas.edu/SABv20a_nam_jgf_status/namforecast/'

bash compute_geotiffs.sh  $url

