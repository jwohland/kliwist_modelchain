# Download LUH1 data, unpack it and remove files that are not needed

# make directory if it does not exist
mkdir -p ../data/LUH
# download historic data
wget -P ../data/LUH https://luh.umd.edu/luh_data/LUHa.v1/LUHa.v1.tgz
# download future data
wget -P ../data/LUH https://luh.umd.edu/luh_data/LUHa.v1_future.v1.1/IMAGE/LUHa.v1_image.v1.1.tgz
wget -P ../data/LUH https://luh.umd.edu/luh_data/LUHa.v1_future.v1/MiniCAM/LUHa.v1_minicam.v1.tgz
wget -P ../data/LUH https://luh.umd.edu/luh_data/LUHa.v1_future.v1/MESSAGE/LUHa.v1_message.v1.tgz
# download readme
wget -P ../data/LUH https://luh.umd.edu/luh_data/LUHa.v1/readme.txt
# unpack all downloaded files
mkdir -p ../data/LUH/historic/
tar zxvf ../data/LUH/LUHa.v1.tgz -C ../data/LUH/historic/
mkdir -p ../data/LUH/IMAGE/
tar zxvf ../data/LUH/LUHa.v1_image.v1.1.tgz -C ../data/LUH/IMAGE/
mkdir -p ../data/LUH/MINICAM/
tar zxvf ../data/LUH/LUHa.v1_minicam.v1.tgz -C ../data/LUH/MINICAM/
mkdir -p ../data/LUH/MESSAGE/
tar zxvf ../data/LUH/LUHa.v1_message.v1.tgz -C ../data/LUH/MESSAGE/
# remove compressed directories
rm -r ../data/LUH/*.tgz
# remove folders lu which are not needed
rm -r ../data/LUH/historic/lu
rm -r ../data/LUH/IMAGE/lu
rm -r ../data/LUH/MINICAM/lu
rm -r ../data/LUH/MESSAGE/lu
# remove historic files pre-1985
rm ../data/LUH/historic/updated_states/*{1500..1984}.txt
# remove future files except 2080-2100
rm ../data/LUH/*/updated_states/*{2005..2079}.txt
# remaining dataset could be further reduced in size by only keeping "gothr" and "gsecd" but will keep other variables
# for potential future use. Dataset about 5GB in size currently.