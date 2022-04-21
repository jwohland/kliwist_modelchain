import matplotlib.pyplot as plt
from compute_country_aggregates import COUNTRIES
import xarray as xr
import seaborn as sns
import pandas as pd


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


def make_scatter_plot(df, filename, metric, country, experiment_id):
    sns.set_theme(style="whitegrid", palette="muted")
    ax = sns.swarmplot(data=df, x="experiment_family", y="sfcWind", hue="GCM", size=10)
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
        meanprops={"color": "k", "ls": "-", "lw": 2, "zorder":10},
        showbox=False,
        showcaps=False,
    )
    ax.axhline(y=0, color="grey", ls="--")
    if metric == "diff":
        ylabel = "Wind speed change [m/s] \n 2080-2100 minus 1985-2005"
    else:
        ylabel = "Wind speeds [m/s] "
    ax.set(
        ylabel=ylabel,
        xlabel="",
        title=country + ", " + experiment_id,
    )
    sns.move_legend(
        ax, "lower center", bbox_to_anchor=(0.4, -0.5), ncol=3, title=None, frameon=True
    )
    plt.subplots_adjust(bottom=0.3, left=0.2)
    plt.savefig(
        "../plots/countries/" + filename + ".jpeg",
        dpi=300,
    )
    plt.close("all")

def make_s10_scatter():
    for metric in ["diff", "mean"]:
        if metric == "diff":
            experiments = ["rcp26", "rcp45", "rcp85"]
        else:
            experiments = ["historical", "rcp26", "rcp45", "rcp85"]
        for experiment_id in experiments:
            name = metric + "_" + experiment_id
            diff_agg_CMIP5 = xr.open_dataset(
                "../output/country_aggregates/country_cmip5_" + name + ".nc"
            )
            diff_agg_CORDEX = xr.open_dataset(
                "../output/country_aggregates/country_cordex_" + name + ".nc"
            )

            for country in COUNTRIES:
                # prepare data as pandas dataframe for plotting with seaborn

                # CMIP5
                df_cmip5 = diff_agg_CMIP5.sel({"country": country}).to_dataframe()
                df_cmip5["GCM"] = [x.split(".")[1] for x in df_cmip5.index]
                df_cmip5["experiment_family"] = "CMIP5"

                # CORDEX
                df_cordex = diff_agg_CORDEX.sel({"country": country}).to_dataframe()
                df_cordex["GCM"] = [x.split(".")[1] for x in df_cordex.index]
                df_cordex["GCM"] = [replace_long_gcm(x, df_cmip5) for x in df_cordex["GCM"]]
                df_cordex["RCM"] = [x.split(".")[-2] for x in df_cordex.index]
                df_cordex["experiment_family"] = "CORDEX"

                # combine both
                df_combined = pd.concat([df_cmip5, df_cordex])

                # now use seaborn to generate scatterplots
                filename = "scatter_" + metric + "_" + experiment_id + "_" + country
                make_scatter_plot(df_combined, filename, metric, country, experiment_id)