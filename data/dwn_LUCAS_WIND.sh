#!/bin/ksh
# download remo data from archieve
set -exu

#star year
syear=1986
# end year
eyear=2015
# user number
user=062
# experiment number
exp=008
# work directory
wdir=/work/ch0636/g300089/remo_results_$user$exp
# temporay directory
tdir=${wdir}/tmp
# output directory
odir=/scratch/g/g300089/LUCAS_WIND/${user}${exp}
# archive directory
archdir=/arch/ch0636/g300089/exp${user}${exp}/


# set up months and years
months="01 02 03 04 05 06 07 08 09 10 11 12"
years=`seq $syear $eyear`
# make directories
mkdir -p $store
mkdir -p tmp
mkdir -p out

cd $tdir

# loop over years
for year in $years ; do
# loop over monhts
    for month in $months ; do
    # t-file
    tfile=e${user}${exp}t${year}${month}
    # file directory
    fdir=$archdir/year${year}
    # retrieve data from archive with slk
    slk retrieve $fdir/$tfile.tar $tdir
          # rename tar-file to avoid issues with wild cards
          if [[ ! -f ${tfile}_dwn.tar ]];then
             mv $tfile.tar ${tfile}_dwn.tar
    fi
          echo $year $month
    # extract t-files (every 6h one file)
    tar xvf ${tfile}_dwn.tar
          rm ${tfile}_{file}.tar
          # merge t-files and convert to netcdf with remo variable names
          cdo -f nc -t remo mergetime e${user}${exp}t${year}${month}???? e${user}${exp}t${year}${month}.nc
    # select variables
    cdo selcode,129,130,131,132,134,156 e${user}${exp}t${year}${month}.nc $odir/e${user}${exp}t${year}${month}_wnd.nc
          # remove temp-files
    rm e${user}${exp}t${year}${month}????
    rm e${user}${exp}t${year}${month}.nc
    done
done
    
