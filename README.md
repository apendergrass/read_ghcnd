# read_ghcnd
Read GHCN-Daily data

Before starting, you'll need to download the data from: https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/
* ```ghcnd_gsn.tar.gz```
* ```ghcnd-stations.txt```
* ```ghcnd-inventory.txt```
  
Then, unzip the data, and go to the directory gsn/ghcnd_gsn/
On linux or mac,  make a text file with the names of all the data files in it with the following command:
```ls -1 > fnames.txt```

In the matlab script, change this to your local directory
```directory='/path/to/files/'```

