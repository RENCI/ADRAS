#!/bin/bash

#url='https://fortytwo.cct.lsu.edu:443/thredds/fileServer/2021/al03/06/LA_v20a-WithUpperAtch_chk/supermic.hpc.lsu.edu/LAv20a_al032021_jgf_23kcms/nhcConsensus/'
url='http://tds.renci.org:8080/thredds/fileServer/2021/nam/2021060918/ec95d/hatteras.renci.org/ec95d-nam-bob-rptest/namforecast/'
bash compute_geotiffs.sh  $url

