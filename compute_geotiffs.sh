#!/bin/bash
#------------------------------------------------------------------------
# compute_geotiffs.sh: computes inundation geotiffs
#------------------------------------------------------------------------
THIS="output/compute_geotiffs.sh"

DEBUG=true

#export GDAL_DATA=/usr/share/gdal/
#export GDAL_DATA=/opt/miniconda2/share/gdal/
#export PYTHONPATH=$HOME/GitHub/RENCI/ADRAS
export PYTHONPATH=/repo/ADRAS
#PKLDIR="$HOME/GitHub/RENCI/ADRAS/pklfiles"
PKLDIR="/repo/ADRAS/pklfiles"
export ADRASHOME=$PYTHONPATH
#export VENV="$HOME/geotiff_p2"
#source $VENV/bin/activate 
#source $HOME/miniconda2/etc/profile.d/conda.sh
#conda activate geotiff_p2
#PYTHON="/home/bblanton/miniconda2/envs/geotiff_p2/bin/python"  # `which python`
#PYTHON=`which python`
PYTHON="/opt/conda/bin/python"

log="log.hazus"

printf "\n\n\n******************************************\n"  >> $log
echo "Compute_geotiffs started at " `date -u`  >> $log

if [[ $DEBUG == "true" ]] ; then
    echo "\$PYTHONPATH =" $PYTHONPATH  >> $log
    echo "\$GDAL_DATA  =" $GDAL_DATA  >> $log
    echo "\$PYTHON     =" $PYTHON  >> $log
    echo "\$PKLDIR     =" $PKLDIR  >> $log
fi

source ./rasterParameters.sh
source ./properties.sh

## main starts here

other="None"

# parse the commandline argument
if [[ $# -lt 1 ]] || [[  $# -gt 1 ]] ; then
    echo "Need only path to run properties file or downloadurl file on command line" >> $log
    echo "Compute_geotiffs terminated at " `date -u`  >> $log
    echo "\n\******************************************\n"  >> $log
    exit 1
fi
temp=$1
if [ ${temp:0:4} == "http" ] ; then
    # assume its a download url, without run.properties at the end
    RUNPROPERTIES="run.properties"
    wget "$temp/$RUNPROPERTIES" --output-document="$RUNPROPERTIES"  2> /dev/null
    if [ $? -ne 0 ] ; then
        echo "wget of $temp/run.properties failed."
        exit 1
   fi
else
    RUNPROPERTIES="$temp"
fi

# load run.properties file into associative array
if [[ $DEBUG == "true" ]] ; then
   echo "\$RUNPROPERTIES=$RUNPROPERTIES"  >> $log
   echo "Loading properties."  >> $log
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

if [ ! -d "$PKLDIR" ]; then
    echo "\nmaking pkldir $PKLDIR " >> $log
    mkdir -p "$PKLDIR"
fi

if [[ $DEBUG == "true" ]] ; then
   printf "\ngridname=$gridname" >> $log
   echo "gridnameabbrev=$gridnameabbrev" >> $log
   echo "wavemodel=$wavemodel" >> $log
   echo "datetime,adv=$datetime, $adv" >> $log
   echo "operator=$operator" >> $log
   echo "ensname, weathertype, windmodel=$ensname, $weathertype, $windmodel" >> $log
   echo "machine=$machine" >> $log
   echo "url=$url" >> $log
fi

# flush a temporary yaml file of the raster parameters
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
    printf "\nRaster parameters are: " >> $log
    cat $RFILE >> $log
fi

ullo=`echo "scale=0; $upperleft_lo*10/1" | bc`
ulla=`echo "scale=0; $upperleft_la*10/1" | bc`
rasterparams=`printf "%d.%06d.%06d.%d.%d" $res $ullo $ulla $nx $ny`

s3path="$year/$weathertype/$datetime/$adv/$windmodel"

varnames=( "inun_max" "zeta_max" ) 
prodvarnames=( "inunmax" "wlmax" )

k=-1
for v in ${varnames[@]}; do
    k=$((k+1))
    prodvarname=${prodvarnames[$k]}

    productId=`printf "%s_%s_%s_%s_%s_%s_%s_%s_%s_%s_%s.tiff" $datetime $adv \
           $prodvarname $gridnameabbrev $windmodel $wavemodel $ensname \
           $operator $machine $other $rasterparams`
    echo $k, $v, $prodvarname, $productId >> $log

    if [[ $DEBUG == "true" ]] ; then
        echo "s3 path = $s3path" >> $log
        echo "s3 product id = $productId" >> $log
    fi

    com="$PYTHON $ADRASHOME/SRC/ADCIRC2Geotiff_p2.py --pkldir=$PKLDIR --varname=$v \
         --gridname=$gridname --url=$url --tif_filename=$productId --s3path=$s3path \
         --rasterconfigfile=$RFILE"
    echo $com >> $log
    $com  >> $log 2>&1

done
echo "Compute_geotiffs finished at " `date -u` >> $log
printf "\n******************************************\n" >> $log
