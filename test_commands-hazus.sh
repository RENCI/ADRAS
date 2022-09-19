#!/bin/bash

# SPDX-FileCopyrightText: 2022 Renaissance Computing Institute. All rights reserved.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-License-Identifier: LicenseRef-RENCI
# SPDX-License-Identifier: MIT

PYTHONPATH="${PYTHONPATH:-$HOME/GitHub/RENCI/ADRAS}"
export PYTHONPATH
echo $PYTHONPATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/bblanton/lib
echo $LD_LIBRARY_PATH

#url='https://fortytwo.cct.lsu.edu:443/thredds/fileServer/2020/laura/28/LA_v20a-WithUpperAtch_chk/qbc.loni.org/LAv20a_al132020_jgf/nhcConsensus/'
#url='https://fortytwo.cct.lsu.edu:443/thredds/fileServer/2021/al03/06/LA_v20a-WithUpperAtch_chk/supermic.hpc.lsu.edu/LAv20a_al032021_jgf_23kcms/nhcConsensus/'
#url='https://fortytwo.cct.lsu.edu/thredds/fileServer/2020/laura/28/LA_v20a-WithUpperAtch_chk/qbc.loni.org/LAv20a_al132020_jgf/nhcConsensus'
#url='http://tds.renci.org:8080/thredds/fileServer/2021/nam/2021060918/ec95d/hatteras.renci.org/ec95d-nam-bob-rptest/namforecast/'
#url='http://tds.renci.org:8080/thredds/fileServer/2021/nam/2021082006/uriv18/hatteras.renci.org/uriv18-nam-bob-2021/namforecast/'
#url='http://adcircvis.tacc.utexas.edu:8080/thredds/fileServer/asgs/2021/al09/01/LAv21a/frontera.tacc.utexas.edu/LAv21a_al092021_jgf_10kcms/nhcConsensus/'
#url='http://adcircvis.tacc.utexas.edu:8080/thredds/fileServer/asgs/2021/al09/05/LAv20a/frontera.tacc.utexas.edu/LAv20a_al092021_jgf_10kcms/veerRight50/'
#url='http://adcircvis.tacc.utexas.edu:8080/thredds/fileServer/asgs/2021/al09/05/LAv20a/stampede2.tacc.utexas.edu/LAv20a_al092021_bde/nhcConsensus/'
#url="http://tds.renci.org:8080/thredds/fileServer/2022/nam/2022011006/hsofs/hatteras.renci.org/hsofs-nam-bob-2021/namforecast/"
#url="http://tds.renci.org:8080/thredds/fileServer/2021/ida/12/ec95d/hatteras.renci.org/ec95d-al09-bob/nhcOfcl/"
#url='http://tds.renci.org:8080/thredds/fileServer/2021/nam/2021092312/NCSC_SAB_v1.15/hatteras.renci.org/ncsc115-nam-2021/namforecast/'
#url='http://adcircvis.tacc.utexas.edu:8080/thredds/fileServer/asgs/2021/nam/2021092306/SABv20a/frontera.tacc.utexas.edu/SABv20a_nam_jgf_status/namforecast/'
#url='http://adcircvis.tacc.utexas.edu:8080/thredds/fileServer/asgs/2021/nam/2021092906/SABv20a/frontera.tacc.utexas.edu/SABv20a_nam_jgf_status/namforecast/'
#url='http://tds.renci.org/thredds/fileServer/2022/nam/2022042606/ec95d/bridges2.psc.edu/ec95d-nam-bob-psc/namforecast'
#url='http://tds.renci.org/thredds/fileServer/2022/nam/2022042606/ec95d/hatteras.renci.org/ec95d-nam-bob-da-nowcast/nowcast'
#url='https://fortytwo.cct.lsu.edu/thredds/fileServer/2020/laura/28/LA_v20a-WithUpperAtch_chk/qbc.loni.org/LAv20a_al132020_jgf/nhcConsensus'
#url='https://fortytwo.cct.lsu.edu:443/thredds/fileServer/2021/al03/06/LA_v20a-WithUpperAtch_chk/supermic.hpc.lsu.edu/LAv20a_al032021_jgf_23kcms/nhcConsensus/'
#url='http://tds.renci.org/thredds/fileServer/2021/nam/2021060918/ec95d/hatteras.renci.org/ec95d-nam-bob-rptest/namforecast/'
url='http://tds.renci.org/thredds/fileServer/2022/nam/2022042606/ec95d/bridges2.psc.edu/ec95d-nam-bob-psc/namforecast'
#url='http://tds.renci.org/thredds/fileServer/2022/nam/2022050912/NCSC_SAB_v1.23/hatteras.renci.org/ncsc123-nam-sb/namforecast'

bash compute_geotiffs.sh  $url

