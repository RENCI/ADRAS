#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2022 Renaissance Computing Institute. All rights reserved.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-License-Identifier: LicenseRef-RENCI
# SPDX-License-Identifier: MIT

#------------------------------------------------------------------------
# compute_geotiffs.sh: computes geotiffs of ADCIRC max files
#------------------------------------------------------------------------
#set -x
#set -u

####################################

DEBUG=true
log="log.hazus"
WGET='wget --no-check-certificate '

filenames=( "maxele.63.nc" "maxele.63.nc" "maxele.63.nc" "swan_HS_max.63.nc" "maxwvel.63.nc" )
varnames=( "depth" "inun_max" "zeta_max" ) # "swan_HS_max" "wind_max")  
prodvarnames=( "depth" "inunmax" "wlmax" "hsignmax" "windspdmax")
keynames=( "Maximum Water Surface Elevation File Name"
           "Maximum Water Surface Elevation File Name"
           "Maximum Water Surface Elevation File Name"
           "Maximum Significant Wave Height File Name" 
           "Maximum Wind Speed File Name") 

RasterPartameterFileUrl='https://raw.githubusercontent.com/RENCI/ADRAS/main/rasterParameters.sh'

Usage ()
{
    echo "Usage: compute_geotiff.sh --downloadurl <url to run.properties> --datadir <path to output dir>"
    exit 2
}

# must be v4+ of bash
if [  ${BASH_VERSION:0:1} -lt 4 ] ; then
	echo "Must run in Bash version >=4.\n"
	exit 1
fi

printf "******************************************\n"  > $log
echo "Compute_geotiffs started at " `date -u`  >> $log

if [ "$#" -eq 0 ] ; then
        echo "Must have at least --downloadurl <url> on command line"
        exit 0
fi

GETOPT='getopt'
if [[ `uname` == "Darwin" ]]; then
        #GETOPT='/usr/local/Cellar/gnu-getopt/1.1.6/bin/getopt'
        GETOPT='/usr/local/opt/gnu-getopt/bin/getopt'
fi

OPTS=$($GETOPT -o v -n compute_geotiff --long downloadurl:,datadir:,verbose -n 'parse-options' -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
    #echo "here"
    Usage
fi

eval set -- "$OPTS"

if [ $? != "0" ]
then
    echo "Failed to parse commandline."
    exit 1
fi
downloadurl="run.properties"
datadir="./"
VERBOSE=true

while true ; do
    case "$1" in
        -v | --verbose)  VERBOSE=true;    shift;;
        --downloadurl)   downloadurl=$2;  shift 2;;
        --datadir)       datadir=$2;      shift 2;;
        --) shift; break;;
        *) echo "Unexpected option: $1 - this should not happen."
           Usage ;;
    esac
done

if [ "$VERBOSE" == true ]; then
        echo VERBOSE = $VERBOSE
        echo runproperties = $downloadurl
        echo Datadir = $datadir
        echo "Remaining Args:"
        echo $@
fi

PKLDIR="${PKLDIR:-$HOME/GitHub/RENCI/ADRAS/pklfiles}"
if [ ! -d ${PKLDIR} ] ; then
	echo "PKLDIR ${PKLDIR} DNE.  Making it ..." | tee -a $log
	mkdir -p $PKLDIR
fi

PYTHONPATH="${PYTHONPATH:-$HOME/GitHub/RENCI/ADRAS}"
if [ ! -d ${PYTHONPATH} ] ; then
	echo "PYTHONPATH ${PYTHONPATH} DNE.  Terminal." | tee -a $log
	exit 1
fi
export ADRASHOME=$PYTHONPATH

temp=`which python`
PYTHON="${PYTHON:-$temp}"
if [ ! -e ${PYTHON} ] ; then
	echo "PYTHON ${PYTHON} DNE.  Terminal." | tee -a $log
	exit 1
fi
export PYTHON=$PYTHON

if [[ $DEBUG == "true" ]] ; then
    echo "\$PYTHONPATH =" $PYTHONPATH  | tee -a $log
    echo "\$PYTHON     =" $PYTHON  | tee -a $log
    echo "\$PKLDIR     =" $PKLDIR  | tee -a $log
fi

RASTERPARAMFILE="./rasterParameters.sh"
#$WGET $RasterPartameterFileUrl  -O realtimeparams.sh # 2> /dev/null
#echo $?
#if [ $? -eq 0 ] ; then
#	RASTERPARAMFILE="./realtimeparams.sh"
#else
#    echo "$WGET of rasterParameters failed. Using defaults."
#fi
source $RASTERPARAMFILE
source ./properties.sh

## main starts here

# parse the commandline argument
#if [[ $# -lt 1 ]] || [[  $# -gt 1 ]] ; then
#    echo "Need only path to run properties file or downloadurl file on command line" | tee -a $log
#    echo "Compute_geotiffs terminated at " `date -u`  | tee -a $log
#    echo "\n\******************************************\n"  | tee -a $log
#    exit 1
#fi
RUNPROPERTIES=$downloadurl
if [ ${RUNPROPERTIES:0:4} == "http" ] ; then
    # assume its a download url, without run.properties at the end
    $WGET "$RUNPROPERTIES/run.properties" --output-document="run.properties" # 2> /dev/null
    if [ $? -ne 0 ] ; then
        echo "$WGET of $RUNPROPERTIES/run.properties failed." | tee -a $log
        exit 1
    fi
    RUNPROPERTIES="run.properties"
fi
echo "\$RUNPROPERTIES=$RUNPROPERTIES" 

# load run.properties file into associative array
#if [[ $DEBUG == "true" ]] ; then
#   echo "\$RUNPROPERTIES=$RUNPROPERTIES" | tee -a $log
#   echo "Loading properties."  | tee -a $log
#fi
declare -A properties
loadProperties $RUNPROPERTIES

#####################################################
# get needed parameters out of run.properties array #
#####################################################
gridname=${properties['adcirc.gridname']}
case $gridname in
   "ec95d")
      gridnameabbrev=$gridname
   ;;
   "nc_inundation_v9.99_w_rivers")
      gridnameabbrev="ncv999wr"
   ;;
  *)
    gridnameabbrev=$(echo $gridname | sed 's/_//g' | sed 's/\.//g')
  ;;
esac

if [[ ! -z ${properties[forcing.metclass]} ]] ; then
   weathertype=${properties[forcing.metclass]}
else
   weathertype='tropical'
fi
temp=${properties['coupling.waves']} 
wavemodel=$([ "$temp" == 'on' ] && echo "swan" || echo "None")
if [ $weathertype == "synoptic" ] ; then
    datetime="20"${properties['currentdate']} 
    adv=${properties['currentcycle']}"Z"
    year=${properties['forcing.nwp.year']}
else
    datetime="al"${properties['forcing.tropicalcyclone.stormnumber']} 
    adv=${properties['advisory']}
    year=${properties['forcing.tropicalcyclone.year']}
fi
temp=${properties['operator']}
operator=$([ "$temp" == 'ncfs-dev' ] && echo "bob" || echo "$temp")
ensname=${properties['asgs.enstorm']}

if [[ ! -z ${properties[forcing.nwp.model]} ]] ; then
   windmodel=${properties[forcing.nwp.model]}
else
   windmodel=${properties[forcing.tropicalcyclone.vortexmodel]}
fi

machine=${properties['hpc.hpcenvshort']}
url=${properties['downloadurl']}
urlbase=${url/fileServer/dodsC}

if [[ $DEBUG == "true" ]] ; then
	printf "\n"                             | tee -a $log
	echo "gridname       = $gridname"       | tee -a $log
	echo "gridnameabbrev = $gridnameabbrev" | tee -a $log
	echo "wavemodel      = $wavemodel"      | tee -a $log
	echo "datetime       = $datetime"       | tee -a $log
	echo "adv            = $adv"            | tee -a $log
	echo "operator       = $operator"       | tee -a $log
	echo "ensname        = $ensname"        | tee -a $log
	echo "weathertype    = $weathertype"    | tee -a $log
	echo "windmodel      = $windmodel"      | tee -a $log
	echo "machine        = $machine"        | tee -a $log
	echo "urlbase        = $urlbase"        | tee -a $log
fi

##############################
# get raster file parameters #
##############################
rasterParameters $gridname

########################################################
# write a temporary yaml file of the raster parameters #
########################################################
RFILE="raster.yml"
rm -rf $RFILE
echo "REGRID: &regrid" > $RFILE
echo "  $gridname:" >> $RFILE
echo "    center_lo: $center_lo" >> $RFILE
echo "    center_la: $center_la" >> $RFILE
echo "    res: $res  # resolution in m" >> $RFILE
echo "    theta: $theta"                >> $RFILE
echo "    nx: $nx"                      >> $RFILE
echo "    ny: $ny"                      >> $RFILE
echo "    target_crs: '$target_crs'"    >> $RFILE
echo "    adcirc_crs: '$adcirc_crs'"    >> $RFILE

ullo=`echo "scale=0; $center_lo*10/1" | bc`
ulla=`echo "scale=0; $center_la*10/1" | bc`
rasterparams=`printf "%d.%06d.%06d.%d.%d" $res $ullo $ulla $nx $ny`

if [[ $DEBUG == "true" ]] ; then
    printf "\nRaster parameters are: \n" | tee -a $log
    cat $RFILE | tee -a $log
fi

s3path="$year/$weathertype/$datetime/$adv/$windmodel"
if [[ $DEBUG == "true" ]] ; then
    echo "s3 path = $s3path" | tee -a $log
fi

other="mMSL"
k=-1
for v in ${varnames[@]}; do

    k=$((k+1))
    prodvarname=${prodvarnames[$k]}
    filename=${filenames[$k]}
    keyname=${keynames[$k]}
    if [[ -z ${properties[$keyname]} ]]; then 
        echo "rp key for $filename DNE.  Skipping..."
        continue
    fi
	if [[ x"$prodvarname" == "xwindspdmax" ]] ; then
		other="mps"
	fi
	if [[ x"$prodvarname" == "xinunmax" ]] ; then
		other="ftMSL"
	fi
    productId=`printf "%s_%s_%s_%s_%s_%s_%s_%s_%s_%s_%s.tiff" $datetime $adv \
               $prodvarname $gridnameabbrev $windmodel $wavemodel $ensname \
               $operator $machine $other $rasterparams`
    printf "\nProcessing %s %s %s %s ... \n" $k $v $prodvarname $productId | tee -a $log
    echo "s3 product id = $productId" | tee -a $log

    url=$urlbase/"$filename"
    com="$PYTHON $ADRASHOME/SRC/ADCIRC2Geotiff.py --pkldir=$PKLDIR --varname=$v \
         --gridname=$gridname --url=$url --tif_filename=$productId --s3path=$s3path \
         --rasterconfigfile=$RFILE --datadir=$datadir"

    echo $com | tee -a $log

    $com  | tee -a $log 2>&1

    printf "\n" | tee -a $log

done

echo "Compute_geotiffs finished at " `date -u` | tee -a $log
printf "\n******************************************\n\n\n" | tee -a $log

