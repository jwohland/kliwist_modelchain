# kliwist_modelchain
Evaluating the GCM-RCM modeling chain for wind energy applications.



##  Data access

Most data is accessed using the `intake` package in the DKRZ ecosystem. Catalogue data is made available in the `output` folder to enable reproducibility at other institutions as well. 

### REMO LUCAS

Idealized REMO experiments from the LUCAS project (assuming all of Europe is covered in forest, grass or current land use) are retrieved from the DKRZ tape archive. 
Retrieval and preprocessing is performed by running:

`bash data/dwn_LUCAS_WIND_eval.sh`

`bash data/dwn_LUCAS_WIND_forest.sh`

`bash data/dwn_LUCAS_WIND_grass.sh`

### NorESM1-M

is not available at DKRZ. 
The model didn't make wind speeds (`sfcWind`) available but only wind components (`ua` and `va`). 
The lowest model level is substantially higher than 10m, namely around 80m. 
For lack of comparability, NorESM1-M is thus excluded here. 


### Shapefile of exclusive economic zones (EEZ)

The offshore assessment relies on the shapes of EEZ, in particular the World EEZ v11 (2019-11-18) shapefile that can be downloaded from https://www.marineregions.org/downloads.php 

Download and extract the data to `data/EEZ/` and remove everything except for the LICENSE and `eez_v11.shp`.
