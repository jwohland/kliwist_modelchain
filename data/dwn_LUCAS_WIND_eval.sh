#!/bin/bash
# download remo data from archive
# 062010 - LUCAS EVAL
set -exu
module load slk  # need to run slk login at least once before using this script

# start year
syear=1986
# end year
eyear=2015
# user number
user=062
# experiment number
exp=010
# output directory
odir=/scratch/g/g300106/LUCAS_WIND/${user}${exp}


# set up months and years
months="01 02 03 04 05 06 07 08 09 10 11 12"
years=`seq $syear $eyear`
# make directory
mkdir -p $odir
mkdir -p $odir/preprocessed
lfs setstripe -E 1G -c 1 -S 1M -E 4G -c 4 -S 1M -E -1 -c 8 -S 1M ${odir}  # recommendation from https://docs.dkrz.de/doc/datastorage/hsm/retrievals.html to increase speed

##########################
# RETRIEVAL
##########################

# search for all t-files
search_id_raw=`slk_helpers search_limited '{"$and": [{"path": {"$gte": "/arch/ch0636/g300089/exp062010"}}, {"resources.name": {"$regex": "t.......tar$"}}]}'`
search_id=`echo $search_id_raw | tail -n 1 | sed 's/[^0-9]*//g'`
echo "The search ID is ${search_id}"
# retrieve
slk retrieve ${search_id} $odir
# '$?' captures the exit code of the previous command (you can put it in
# the next line after each slk command).
if [ $? -ne 0 ]; then
    >&2 echo "an error occurred in slk retrieve call"
else
    echo "retrieval successful"
fi


##########################
# PREPROCESSING
##########################

cd $odir
# loop over years
for year in $years ; do
# loop over months
    for month in $months ; do
    # t-file
    tfile=e${user}${exp}t${year}${month}
    # rename tar-file to avoid issues with wild cards
    if [[ ! -f ${tfile}_dwn.tar ]];then
       mv $tfile.tar ${tfile}_dwn.tar
    fi
      echo $year $month
    # extract t-files (every 6h one file)
    tar xvf ${tfile}_dwn.tar
    # merge t-files and convert to netcdf with remo variable names
    cdo -f nc -t remo mergetime e${user}${exp}t${year}${month}???? e${user}${exp}t${year}${month}.nc
    # select variables and copy to scratch
    cdo selcode,129,130,131,132,134,156 e${user}${exp}t${year}${month}.nc preprocessed/e${user}${exp}t${year}${month}_wnd.nc
    # remove temp-files
    rm e${user}${exp}t${year}${month}????
    rm e${user}${exp}t${year}${month}.nc
    done
done

# remove files that are no longer needed later
#rm ${tfile}_dwn.tar  # skipped for now