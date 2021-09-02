for s in `cat url.list`; do 
	if [ ${s:0:1} == h ] ; then 
		bash compute_geotiffs.sh $s
	else
		echo "not running $s"
	fi
done
