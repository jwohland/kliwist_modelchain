import sys

sys.path.append("../code")
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from plot_s10_maps import add_coast_boarders, SUBPLOT_KW
import cartopy.crs as ccrs


def load_winds_tempgradients(metric="diff"):
    """
    open wind speed change and tas change data across all scenarios and concatenates them together.

    temperatures are output as changes in the equator-to-pole gradient rather than absolute
    :param metric:
    :return:
    """
    path_sfc = "../output/sfcWind/monthly/"
    path_tas = "/work/ch0636/g300106/projects/kliwist_modelchain/output/tas/monthly/"

    wind_list, gradient_list = [], []
    for i, experiment_id in enumerate(["rcp26", "rcp45", "rcp85"]):
        name = metric + "_" + experiment_id
        ds_wind = xr.open_dataset(path_sfc + "cmip5_" + name + ".nc")
        ds_tmp = xr.open_dataset(path_tas + "cmip5_" + name + ".nc")
        T_equator = ds_tmp.sel(lat=slice(-10, 10)).mean(dim=["lat", "lon", "month"])
        T_pole = ds_tmp.sel(lat=slice(70, 90)).mean(dim=["lat", "lon", "month"])
        T_gradient = T_equator - T_pole
        for identifier in ds_wind.identifier.values:
            wind_list.append(
                ds_wind.mean(dim="month")
                .sel(identifier=identifier, drop=True)
                .drop("height")
            )
            gradient_list.append(T_gradient.sel(identifier=identifier, drop=True))
    return xr.concat(wind_list, dim="experiments"), xr.concat(
        gradient_list, dim="experiments"
    )


def prepare_figure(scope):
    """
    Prepare Figure with grid ready to plot either global or European data
    :param scope:
    :return:
    """
    if scope == "Globe":
        f, ax = plt.subplots(
            ncols=1, figsize=((5, 3)), subplot_kw={"projection": ccrs.PlateCarree()}
        )
    else:
        f, ax = plt.subplots(ncols=1, figsize=((5, 3)), **SUBPLOT_KW)
    return f, ax


def correlation_compute_plot(wind_ds, gradient_ds):
    """
    Calculates the correlation between changes in wind gradient and wind speeds and plots
    global correlation map and European correlation map
    :param wind_ds:
    :param gradient_ds:
    :return:
    """
    corr = xr.corr(wind_ds["sfcWind"], gradient_ds["tas"], dim="experiments")
    levels = [-1, -0.8, -0.6, 0.6, 0.8, 1]  # mask out correlation smaller than r=0.6
    # global correlation plot

    for scope in ["Globe", "Europe"]:
        f, ax = prepare_figure(scope)
        corr.plot(
            levels=levels,
            ax=ax,
            cbar_kwargs={
                "label": "Correlation temperature gradient change and wind speed change",
                "orientation": "horizontal",
            },
        )
        add_coast_boarders(ax)
        plt.savefig("../plots/tas_gradient/Correlation_map_" + scope + ".jpeg", dpi=300)


def amplitude_compute_plot(wind_ds, gradient_ds, method):
    assert method in ["mean_changes", "extreme"]
    if method == "mean_changes":
        # A) based on mean changes
        delta_wind = wind_ds["sfcWind"].mean(dim="experiments")
        delta_gradient = gradient_ds["tas"].mean(dim="experiments")
        slope_proxy = delta_wind / delta_gradient
        label = "Wind change per gradient change proxy  [m/s /K]"
    elif method == "extreme":
        # B) based on extreme cases
        strongest_change = gradient_ds["tas"].argmin()  # minimum is correct because gradient change negative
        slope_proxy = wind_ds["sfcWind"].sel(experiments=strongest_change)
        label = "Wind change in max gradient change experiment [m/s]"
    #elif method == "regression":
    #    print("Not implemented yet")
    #    fit = xr.polyfit(wind_ds["sfcWind"], gradient_ds["tas"], dim="experiments")
    #    return
    for scope in ["Globe", "Europe"]:
        f, ax = prepare_figure(scope)
        slope_proxy.plot(
            ax=ax,
            levels=np.linspace(-0.1, 0.1, 11),
            extend="both",
            cbar_kwargs={
                "label": label,
                "orientation": "horizontal",
            },
        )
        add_coast_boarders(ax)
        plt.savefig(
            "../plots/tas_gradient/Amplitude_map_" + method + "_" + scope + ".jpeg",
            dpi=300,
        )


def scatter_compute_plot_country(country, offshore=True, metric="diff"):
    """
    Open wind speed change information per country and plots scatter plot of wind speed change
    vs. equator to pole temperature gradient change.
    :param country:
    :param offshore:
    :param metric:
    :return:
    """
    """
    Country level scatter plots
    :param country:
    :param offshore:
    :param metric:
    :return:
    """
    if offshore:
        path_sfc = "../output/sfcWind/country_aggregates/offshore/"
    else:
        path_sfc = "../output/sfcWind/country_aggregates/"
    path_tas = "../output/tas/monthly/"
    wind_df_list, tas_df_list = [], []

    for experiment_id in ["rcp26", "rcp45", "rcp85"]:
        # wind
        name = metric + "_" + experiment_id
        ds_wind = xr.open_dataset(path_sfc + "country_cmip5_" + name + ".nc")
        if country == "all":
            ds_wind = ds_wind.mean(dim="country").drop(["height"])
        else:
            ds_wind = ds_wind.sel(country=country).drop(
                ["height", "lon", "lat", "country"]
            )
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
    tas_df.index = pd.Index(
        [".".join(x.split(".")[:2]) for x in tas_df.index], name="identifier"
    )
    wind_df.index = pd.Index(
        [".".join(x.split(".")[:2]) for x in wind_df.index], name="identifier"
    )
    df = pd.merge(tas_df, wind_df, on=["identifier", "experiment_id"])

    sns.scatterplot(x="tas", y="sfcWind", data=df, hue="experiment_id")
    plt.title(
        "Correlation coefficient is " + str(np.round(df["tas"].corr(df["sfcWind"]), 2))
    )
    plt.ylabel(country + " offshore wind speed change [m/s]")
    plt.xlabel("Equator minus pole temperature change [K]")

    plt.savefig(
        "../plots/tas_gradient/Scatter_plot_"
        + country
        + "_offshore_"
        + offshore
        + ".jpeg",
        dpi=300,
    )


def make_all_plots():
    wind_ds, gradient_ds = load_winds_tempgradients()

    # correlation maps
    correlation_compute_plot(wind_ds, gradient_ds)

    # Proxies for amplitude of change
    amplitude_compute_plot(wind_ds, gradient_ds, method="mean_changes")
    amplitude_compute_plot(wind_ds, gradient_ds, method="extremes")

    # Country level scatter plots
    for country in ["Norway", "United Kingdom", "Ireland", "Germany", "Portugal", "all"]:
        for offshore in [True, False]:
            scatter_compute_plot_country(country, offshore)
    # add couple of countries without coast
    for country in ["Czechia", "Hungary"]:
        scatter_compute_plot_country(country, offshore=False)