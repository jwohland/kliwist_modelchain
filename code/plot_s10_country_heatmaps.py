import xarray as xr
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from numpy import linspace
from plot_s10_maps import xarray_ttest

def make_s10_heatmaps(onshore):
    metric = "diff"
    df_list = []
    path_to_data = "../output/sfcWind/country_aggregates/"
    if not onshore:
        path_to_data += "offshore/"
    for experiment_id in ["rcp26", "rcp45", "rcp85"]:
        name = metric + "_" + experiment_id
        for experiment_family in ["CORDEX", "CMIP5"]:
            diff_agg = xr.open_dataset(
                path_to_data
                + "country_"
                + experiment_family.lower()
                + "_"
                + name
                + ".nc"
            )
            p_value = xarray_ttest(diff_agg)
            diff_df = (
                diff_agg.mean(dim="identifier")
                .drop(["height", "lat", "lon", "member"], errors="ignore")
                .to_dataframe()
            )  # errors=ignore means that no error is raised if subset of height, lat, lon does not exist
            diff_df["experiment_family"] = experiment_family
            diff_df["experiment_id"] = experiment_id
            diff_df["p_value"] = p_value
            df_list.append(diff_df)
    df = pd.concat(df_list)

    f, axs = plt.subplots(ncols=2, figsize=(12, 9))
    cbar_ax = f.add_axes([0.2, 0.05, 0.6, 0.02])

    for i, experiment_family in enumerate(["CORDEX", "CMIP5"]):
        df_tmp = df[df.experiment_family == experiment_family]
        df_wind = df_tmp.pivot(columns="experiment_id", values="sfcWind")
        df_p = df_tmp.pivot(columns="experiment_id", values="p_value")

        # plot significant changes in bold
        sns.heatmap(
            df_wind,
            ax=axs[i],
            annot=True,
            fmt=".2f",
            annot_kws={"weight": "bold"},
            cmap=plt.get_cmap("RdBu_r", 5),
            vmin=-.3, 
            vmax=.3,
            mask=df_p>0.05,  # mask countries where N0 not rejected
            cbar_kws={
                "orientation": "horizontal",
                "label": "Ensemble mean wind speed change [m/s]",
                "ticks": linspace(-.3, .3, 6),
            },
            cbar_ax=cbar_ax,
        )
        # add non-significant changes in normal font
        sns.heatmap(
            data=df_wind,
            ax=axs[i],
            annot=True,
            fmt=".2f",
            cmap=plt.get_cmap("RdBu_r", 5),
            vmin=-.3,
            vmax=.3,
            mask=df_p<0.05,  # mask countries where N0 rejected
            cbar_kws={
                "orientation": "horizontal",
                "label": "Ensemble mean wind speed change [m/s]",
                "ticks": linspace(-.3, .3, 6),
            },
            cbar_ax=cbar_ax,
        )
        title_dic = {"CORDEX": "EURO-CORDEX", "CMIP5": "CMIP5"}
        axs[i].set(xlabel="", ylabel="", title=title_dic[experiment_family])

    plt.subplots_adjust(left=0.11, right=0.98, top=0.97, wspace=0.3)
    figname = "heatmap_mean_countries"
    if not onshore:
        figname += "_offshore"
    plt.savefig(
        "../plots/countries/" + figname + ".jpeg",
        dpi=300,
    )
