import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

def plot_wind_temp_gradient(country, offshore=True, metric="diff"):
    if offshore:
        path_sfc = "../output/sfcWind/country_aggregates/offshore/"
    else:
        path_sfc = "../output/sfcWind/country_aggregates/"
    path_tas = "/work/ch0636/g300106/projects/kliwist_modelchain/output/tas/monthly/"
    wind_df_list=[]
    tas_df_list=[]
    for experiment_id in ["rcp26", "rcp45", "rcp85"]:
        # wind
        name = metric + "_" + experiment_id
        ds_wind = xr.open_dataset(
            path_sfc + "country_cmip5_" + name + ".nc"
        )
        if country == "all":
            ds_wind = ds_wind.mean(dim="country").drop(["height"])
        else:
            ds_wind = ds_wind.sel(country=country).drop(["height", "lon", "lat", "country"])
        ds_wind["experiment_id"] = experiment_id
        wind_df_list.append(ds_wind.to_dataframe())

        # temperature gradient
        ds_tmp = xr.open_dataset(path_tas + "cmip5_" + name + ".nc")
        T_equator = ds_tmp.sel(lat=slice(-10, 10)).mean(dim=["lat", "lon", "month"])
        T_pole = ds_tmp.sel(lat=slice(70, 90)).mean(dim=["lat", "lon", "month"])
        T_gradient = T_equator - T_pole
        T_gradient = T_gradient.drop("height")
        T_gradient["experiment_id"] = experiment_id
        tas_df_list.append(T_gradient.to_dataframe())

    tas_df = pd.concat(tas_df_list)
    wind_df = pd.concat(wind_df_list)
    tas_df.index = pd.Index([".".join(x.split(".")[:2]) for x in tas_df.index], name="identifier")
    wind_df.index = pd.Index([".".join(x.split(".")[:2]) for x in wind_df.index], name="identifier")
    df = pd.merge(tas_df, wind_df, on=["identifier", "experiment_id"])

    sns.scatterplot(x="tas", y="sfcWind", data=df, hue="experiment_id")
    plt.title("Correlation coefficient is " + str(np.round(df["tas"].corr(df["sfcWind"]),2)))
    plt.ylabel(country + " Offshore Wind Speed change [m/s]")
    plt.xlabel("Equator minus Pole Temperature change [K]")


plot_wind_temp_gradient("Norway")
plot_wind_temp_gradient("Norway", offshore=False)


#todo save Figures