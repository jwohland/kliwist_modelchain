# kliwist_modelchain

Code underlying analysis performed in 

> Jan Wohland, Process-based climate change impact assessment for European winds using EURO-CORDEX and global models, Environmental Research Letters (provisionally accepted), 2022

If you use content of this repository or code derived from it in academic work, please cite the above publication. 

The intention of this repository is to document the analysis in an attempt to make scientific work more transparent and reproducible. 


## Access to raw input data 

### Climate model data

Most data is accessed using the `intake` package in the DKRZ ecosystem. Intermediate data is made available in the `output` folder to enable reproducibility at other institutions as well, see below for the link and instructions. 

### Land use change data

Land use change data is taken from LUH and can be retrieved by executing `download_LUH1.sh`. 

### Shapefile of exclusive economic zones (EEZ)

The offshore assessment relies on the shapes of EEZ, in particular the World EEZ v11 (2019-11-18) shapefile provided by the Flanders Marine Institute and available at https://doi.org/10.14284/386

Download and extract the data to `data/EEZ/` and remove everything except for the LICENSE and `eez_v11.shp`.

## Access to intermediate data

Intermediate data (i.e., results from running everything under *Calculations* in `run_all.py`) are provided in a zenodo data repository. Download the data from https://doi.org/10.5281/zenodo.7372998, extract the zipped `output.zip` folder and move it into the same folder as `code`.

## Anaconda environment

The anaconda environment can be constructed from the provided files using

`conda env create --file environment.yaml`

You can activate it using 

`conda activate kliwist_modelchain`

on a UNIX system. You might have to use `source` on other systems. 

## Figure overview
| Figure | Filename | Creating python function |
|---|---|---|
| Fig. 1| windchange_mean.png | `plot_s10_maps.make_joint_plots()` |
| Fig. 2 a-c | diff_landuse.png | `plot_lu_maps.make_LUH_maps()` |
| Fig. 2d | pattern_correlation.png | `compute_plot_pattern_correlation.make_plot()` |
| Fig. 3 | heatmap_mean_countries.jpeg | `plot_s10_country_heatmaps.make_s10_heatmaps(onshore=True)` |
| Fig. 4 | heatmap_mean_countries_offshore.jpeg | `plot_s10_country_heatmaps.make_s10_heatmaps(onshore=False)` |
| Fig. 5 a-c | scatter_diff_United Kingdom_offshore.jpeg | `plot_s10_scatter.make_s10_scatter(onshore=False)` |
| Fig. 5d | Scatter_plot_United Kingdom_offshore_True.jpeg | `plot_temperature_gradient.make_all_plots()` |
| Fig. 6a | Correlation_map_Europe.jpeg | `plot_temperature_gradient.make_all_plots()` |
| Fig. 6b | Amplitude_map_Europe.jpeg | `plot_temperature_gradient.make_all_plots()` |


All Figures can be created at once by executing the script `run_all.py` or the Notebook `run_all.ipynb`. Script and notebook execute the same code, so feel free to choose whichever route you prefer. 
