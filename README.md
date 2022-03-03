# kliwist_modelchain
Evaluating the GCM-RCM modeling chain for wind energy applications.



##  Data access

Most data is accessed using the `intake` package in the DKRZ ecosystem. Catalogue data is made available in the `output` folder to enable reproducibility at other institutions as well. 

### NorESM1-M

is not available at DKRZ. The model also didn't make wind speeds (`sfcWind`) available but only wind components (`ua` and `va`). Components on 6h temporal resolutions are used to compute wind speeds which are then aggregated temporally. The 6h wind components are instantaneous. Downloads have to be done manually using the `wget` scripts provided in `input`. Those scripts were taken from https://aims2.llnl.gov/metagrid/search and then screened manually for `ua` and `va` only (filtering didn't work in the ESGF context for unclear reasons). 
