#!/usr/bin/env bash
#------------------------------------------------------------------------
# compute_geotiffs.sh: computes inundation geotiffs
#------------------------------------------------------------------------
# THIS="output/compute_geotiffs.sh"
#set -x
#set -u

if [  ${BASH_VERSION:0:1} -lt 4 ] ; then
        echo "Must run in Bash version >=4.\n"
        exit 1
fi

DEBUG=true

log="log.hazus"

printf "\n\n\n******************************************\n"  > $log
echo "Compute_geotiffs started at " `date -u`  >> $log

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
#PYTHON="/home/bblanton/miniconda2/envs/geotiff_p2/bin/python" 
export PYTHON=`which python`

if [[ $DEBUG == "true" ]] ; then
    echo "\$PYTHONPATH =" $PYTHONPATH  | tee -a $log
    echo "\$GDAL_DATA  =" $GDAL_DATA   | tee -a $log
    echo "\$PYTHON     =" $PYTHON  | tee -a $log
    echo "\$PKLDIR     =" $PKLDIR  | tee -a $log
fi

source ./rasterParameters.sh
source ./properties.sh

## main starts here

# parse the commandline argument
if [[ $# -lt 1 ]] || [[  $# -gt 1 ]] ; then
    echo "Need only path to run properties file or downloadurl file on command line" | tee -a $log
    echo "Compute_geotiffs terminated at " `date -u`  | tee -a $log
    echo "\n\******************************************\n"  | tee -a $log
    exit 1
fi
RUNPROPERTIES=$1
if [ ${RUNPROPERTIES:0:4} == "http" ] ; then
    # assume its a download url, without run.properties at the end
    wget "$RUNPROPERTIES/run.properties" --output-document="run.properties" 2> /dev/null
    if [ $? -ne 0 ] ; then
        echo "wget of $RUNPROPERTIES/run.properties failed." | tee -a $log
        exit 1
   fi
   RUNPROPERTIES="run.properties"
fi

# load run.properties file into associative array
if [[ $DEBUG == "true" ]] ; then
   echo "\$RUNPROPERTIES=$RUNPROPERTIES" | tee -a $log
   echo "Loading properties."  | tee -a $log
fi

declare -A properties
loadProperties $RUNPROPERTIES

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

rasterParameters $gridname

# get needed parameters out of run.properties array
if [[ ! -z ${properties[forcing.metclass]} ]] ; then
   weathertype=${properties[forcing.metclass]}
else
   weathertype='unknown'
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
url=${url/fileServer/dodsC}"/maxele.63.nc"

if [[ $DEBUG == "true" ]] ; then
   printf "\ngridname=$gridname" | tee -a $log
   echo "gridnameabbrev=$gridnameabbrev" | tee -a $log
   echo "wavemodel=$wavemodel" | tee -a $log
   echo "datetime,adv=$datetime, $adv" | tee -a $log
   echo "operator=$operator" | tee -a $log
   echo "ensname, weathertype, windmodel=$ensname, $weathertype, $windmodel" | tee -a $log
   echo "machine=$machine" | tee -a $log
   echo "url=$url" | tee -a $log
fi

# write a temporary yaml file of the raster parameters
RFILE="raster.yml"
rm -rf $RFILE
echo "REGRID: &regrid" > $RFILE
echo "  $gridname:" >> $RFILE
echo "    upperleft_lo: $upperleft_lo" >> $RFILE
echo "    upperleft_la: $upperleft_la" >> $RFILE
echo "    res: $res  # resolution in m" >> $RFILE
echo "    nx: $nx" >> $RFILE
echo "    ny: $ny" >> $RFILE
echo "    target_crs: '$target_crs'" >> $RFILE
echo "    adcirc_crs: '$adcirc_crs'" >> $RFILE

if [[ $DEBUG == "true" ]] ; then
    printf "\nRaster parameters are: \n" | tee -a $log
    cat $RFILE | tee -a $log
fi

ullo=`echo "scale=0; $upperleft_lo*10/1" | bc`
ulla=`echo "scale=0; $upperleft_la*10/1" | bc`
rasterparams=`printf "%d.%06d.%06d.%d.%d" $res $ullo $ulla $nx $ny`

s3path="$year/$weathertype/$datetime/$adv/$windmodel"

varnames=( "inun_max" "zeta_max" ) 
prodvarnames=( "inunmax" "wlmax" )

other="None"
k=-1
for v in ${varnames[@]}; do
    k=$((k+1))
    prodvarname=${prodvarnames[$k]}

    productId=`printf "%s_%s_%s_%s_%s_%s_%s_%s_%s_%s_%s.tiff" $datetime $adv \
           $prodvarname $gridnameabbrev $windmodel $wavemodel $ensname \
           $operator $machine $other $rasterparams`
    echo $k, $v, $prodvarname, $productId | tee -a $log

    if [[ $DEBUG == "true" ]] ; then
        echo "s3 path = $s3path" | tee -a $log
        echo "s3 product id = $productId" | tee -a $log
    fi

    com="$PYTHON $ADRASHOME/SRC/ADCIRC2Geotiff_p2.py --pkldir=$PKLDIR --varname=$v \
         --gridname=$gridname --url=$url --tif_filename=$productId --s3path=$s3path \
         --rasterconfigfile=$RFILE"
    echo $com | tee -a $log
    $com  | tee -a $log 2>&1

done
echo "Compute_geotiffs finished at " `date -u` | tee -a $log
printf "\n******************************************\n" | tee -a $log
