import matplotlib.pyplot as plt
from compute_country_aggregates import COUNTRIES, OFFSHORE_COUNTRIES
import xarray as xr
import seaborn as sns
import pandas as pd
from plot_utils import *


def replace_long_gcm(longname, df_cmip5):
    """
    Helper function to ensure that GCM names in the CMIP5 and EUROCORDEX
    data are identical. This function is needed because EUROCORDEX uses
    longer names that also include the institution while CMIP5 does not.
    :param longname:
    :param df_cmip5:
    :return:
    """
    shortnames = df_cmip5["GCM"].values
    try:
        new_name = [x for x in shortnames if x in longname][0]
    except IndexError:
        new_name = "-".join(
            longname.split("-")[1:]
        )  # if model doens't exist remove at least institute name
    return new_name


def make_scatter_plot(df, metric, experiment_id, ax):
    sns.set_theme(style="whitegrid", palette="muted")
    hue_order = [
        "CNRM-CM5",
        "EC-EARTH",
        "IPSL-CM5A-MR",
        "HadGEM2-ES",
        "MPI-ESM-LR",
        "MIROC5",
        "GFDL-ESM2G",
        "NorESM1-M",
        "IPSL-CM5A-LR",
    ]  # manually define hue order to ensure that GCMs are always mapped to the same color independent of ensemble size
    sns.swarmplot(
        data=df,
        x="experiment_family",
        y="sfcWind",
        hue="GCM",
        size=8 if metric == "diff" else 3,  # markers need to be smaller for means because they are closer
        ax=ax,
        hue_order=hue_order,
    )
    sns.boxplot(
        data=df,
        x="experiment_family",
        y="sfcWind",
        color="white",
        saturation=0.5,
        width=0.4,
        showmeans=True,
        meanline=True,
        medianprops={"visible": False},
        whiskerprops={"visible": False},
        meanprops={"color": "k", "ls": "-", "lw": 2, "zorder": 10},
        showbox=False,
        showcaps=False,
        ax=ax,
    )
    ax.axhline(y=0, color="grey", ls="--")
    if metric == "diff":
        ylabel = "Wind speed change [m/s]"
    else:
        ylabel = "Wind speeds [m/s] "
    ax.set(
        ylabel=ylabel,
        xlabel="",
        title=experiment_id,
    )
    sns.move_legend(
        ax,
        "lower center",
        bbox_to_anchor=(0.34, -0.3),
        ncol=5,
        title=None,
        frameon=True,
    )


def make_s10_scatter(onshore=True):
    path = "../output/sfcWind/country_aggregates/"
    if onshore:
        relevant_countries = COUNTRIES
    else:
        relevant_countries = OFFSHORE_COUNTRIES
        path += "offshore/"

    for country in relevant_countries:
        # prepare data as pandas dataframe for plotting with seaborn
        for metric in ["diff", "mean"]:
            if metric == "diff":
                experiments = ["rcp26", "rcp45", "rcp85"]
            else:
                experiments = ["historical", "rcp26", "rcp45", "rcp85"]
            f, axs = plt.subplots(ncols=len(experiments), figsize=(10, 5), sharey=True)
            for i, experiment_id in enumerate(experiments):
                name = metric + "_" + experiment_id
                diff_agg_CMIP5 = xr.open_dataset(path + "country_cmip5_" + name + ".nc")
                diff_agg_CORDEX = xr.open_dataset(
                    path + "country_cordex_" + name + ".nc"
                )

                # CMIP5
                df_cmip5 = diff_agg_CMIP5.sel({"country": country}).to_dataframe()
                df_cmip5["GCM"] = [x.split(".")[1] for x in df_cmip5.index]
                df_cmip5["experiment_family"] = "CMIP5"

                # CORDEX
                df_cordex = diff_agg_CORDEX.sel({"country": country}).to_dataframe()
                df_cordex["GCM"] = [x.split(".")[1] for x in df_cordex.index]
                df_cordex["GCM"] = [
                    replace_long_gcm(x, df_cmip5) for x in df_cordex["GCM"]
                ]
                df_cordex["RCM"] = [x.split(".")[-2] for x in df_cordex.index]
                df_cordex["experiment_family"] = "EURO-CORDEX"

                # combine both
                df_combined = pd.concat([df_cmip5, df_cordex])

                # now use seaborn to generate scatterplots
                filename = "scatter_" + metric + "_" + country
                if not onshore:
                    filename += "_offshore"
                make_scatter_plot(df_combined, metric, experiment_id, axs[i])

            plt.subplots_adjust(bottom=0.23, left=0.08, right=0.98, top=0.95)
            add_letters(axs, y=1.04, fs=12)
            for i in range(len(experiments)):
                if i > 0:
                    axs[i].set_ylabel("")  # keep y label only once
                if i != 1:
                    axs[
                        i
                    ].get_legend().remove()  # legends are identical. only keep central one
            plt.savefig(
                "../plots/countries/" + filename + ".jpeg",
                dpi=300,
            )
            plt.close("all")
