import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import xarray as xr
from numpy import unique
from pandas import MultiIndex


def list_of_models(ds, model_type):
    index = -2
    if model_type == "GCM":
        index = 1
    return sorted(list(unique([x.split(".")[index] for x in ds.identifier.values])))


def reindex_per_model(ds):
    tmp = diff.copy()
    RCMs = [x.split(".")[-2] for x in ds.identifier.values]
    GCMs = [x.split(".")[1] for x in ds.identifier.values]
    new_index = MultiIndex.from_arrays([RCMs, GCMs], names=("RCMs", "GCMs"))
    return tmp.assign(identifier=new_index).unstack("identifier")


def plot_array(ds):
    # prepare plotting
    RCMs, GCMs = list_of_models(ds, "RCM"), list_of_models(ds, "GCM")
    f, axs = plt.subplots(
        ncols=len(GCMs),
        nrows=len(RCMs),
        subplot_kw={"projection": ccrs.PlateCarree(), "extent": [-15, 50, 35, 70]},
        figsize=(10, 8),
        dpi=300,
    )
    cbar_ax = f.add_axes([0.2, 0.05, 0.6, 0.01])
    label = "Wind speed change 2080-2100 minus 1985-2005 [m/s]"

    for ident in sorted(ds.identifier.values):
        GCM, RCM = ident.split(".")[1], ident.split(".")[-2]
        ax = axs[RCMs.index(RCM), GCMs.index(GCM)]
        # plot values
        ds["sfcWind"].sel(identifier=ident).plot(
            ax=ax,
            x="lon",
            y="lat",
            cbar_ax=cbar_ax,
            vmin=-0.8,
            vmax=0.8,
            extend="both",
            cmap=plt.get_cmap("coolwarm"),
            cbar_kwargs={"label": label, "orientation": "horizontal"},
        )
        # add hatching
        ds["sfcWind"].sel(identifier=ident).plot(
            ax=ax,
            x="lon",
            y="lat",
            levels=[-10, -0.1, 0.1, 10],
            colors=["none", "white", "none"],
            add_colorbar=False,
        )
        ax.add_feature(cf.COASTLINE)
        ax.add_feature(cf.BORDERS)
        ax.set_extent([-15, 50, 35, 70])
        ax.set_title("")
    for i, GCM in enumerate(GCMs):  # add GCM name as column headings
        axs[0, i].set_title(GCM, fontsize=6)
    for j, RCM in enumerate(RCMs):  # add RCM name as row y axis
        axs[j, 0].text(
            -0.1,
            0.5,
            RCM,
            rotation="vertical",
            fontsize=5,
            horizontalalignment="center",
            verticalalignment="center",
            transform=axs[j, 0].transAxes,
        )
    plt.subplots_adjust(0.04, 0.1, 0.95, 0.97, hspace=0.05, wspace=0.05)


def plot_array_CMIP5(
    ds,
):  # todo currently copied from plot_array need cleaner plotting functions!
    # prepare plotting
    f, axs = plt.subplots(
        ncols=ds.identifier.size,
        nrows=1,
        subplot_kw={"projection": ccrs.PlateCarree(), "extent": [-15, 50, 35, 70]},
        figsize=(8, 2),
        dpi=300,
    )
    cbar_ax = f.add_axes(
        [0.2, 0.3, 0.6, 0.05]
    )  # maybe use inset_locator instead? https://stackoverflow.com/questions/13310594/positioning-the-colorbar
    label = "Wind speed change 2080-2100 minus 1985-2005 [m/s]"

    for i, ident in enumerate(sorted(ds.identifier.values)):
        GCM = ident.split(".")[1]
        # plot values
        ds["sfcWind"].sel(identifier=ident).plot(
            ax=axs.flatten()[i],
            x="lon",
            y="lat",
            cbar_ax=cbar_ax,
            vmin=-0.8,
            vmax=0.8,
            extend="both",
            cmap=plt.get_cmap("coolwarm"),
            cbar_kwargs={"label": label, "orientation": "horizontal"},
        )
        # add hatching
        ds["sfcWind"].sel(identifier=ident).plot(
            ax=axs.flatten()[i],
            x="lon",
            y="lat",
            levels=[-10, -0.1, 0.1, 10],
            colors=["none", "white", "none"],
            add_colorbar=False,
        )

        axs.flatten()[i].set_title(GCM, fontsize=6)

    for ax in axs.flatten():
        ax.add_feature(cf.COASTLINE)
        ax.add_feature(cf.BORDERS)
    plt.subplots_adjust(0.02, 0.15, 0.98, 0.99)


### Individual plots of all models ###
for experiment_id in ["rcp26", "rcp45", "rcp85"]:
    # CORDEX
    diff = xr.open_dataset("../output/cordex_diff_" + experiment_id + ".nc")
    plot_array(diff)
    plt.savefig(
        "../plots/cordex_windchange_" + experiment_id + ".png",
        dpi=300,
        facecolor="w",
        transparent=False,
    )
    # CMIP5
    diff = xr.open_dataset("../output/cmip5_diff_" + experiment_id + ".nc")
    plot_array_CMIP5(diff)
    plt.savefig(
        "../plots/cmip5_windchange_" + experiment_id + ".png",
        dpi=300,
        facecolor="w",
        transparent=False,
    )


### Aggregate plots ###
# CORDEX
RCM_ensemble_size = {"rcp26": 4, "rcp45": 5, "rcp85": 10}  # number of RCMs with >1 downscaled GCM
for experiment_id in ["rcp26", "rcp45", "rcp85"]:
    diff = xr.open_dataset("../output/cordex_diff_" + experiment_id + ".nc")
    ds = reindex_per_model(diff)
    for metric in ["mean", "standard_deviation", "mean_per_std"]:
        label = metric + " wind speed change 2080-2100 minus 1985-2005 [m/s]"
        RCMs = unique(ds.RCMs)
        f, axs = plt.subplots(
            nrows=RCM_ensemble_size[experiment_id],
            subplot_kw={"projection": ccrs.PlateCarree(), "extent": [-15, 50, 35, 70]},
            figsize=(4, len(experiment_id) * 2),
        )
        plt.subplots_adjust(bottom=0.03, top=0.98)
        cbar_ax = f.add_axes([0.1, 0.015, 0.8, 0.01])
        i = 0
        for RCM in RCMs:
            N_models = (
                len(ds.GCMs)
                - ds.sel(RCMs=RCM)
                .mean(dim=["rlat", "rlon"])["sfcWind"]
                .isnull()
                .values.sum()
            )
            if N_models > 1:
                if metric == "mean":
                    plot_data = ds.sel(RCMs=RCM).mean(dim="GCMs", skipna=True)
                    vmin, vmax = -0.5, 0.5
                    cmap = plt.get_cmap("coolwarm")
                elif metric == "standard_deviation":
                    plot_data = ds.sel(RCMs=RCM).std(dim="GCMs", skipna=True)
                    vmin, vmax = 0, 0.25
                    cmap = plt.get_cmap("Greens")
                elif metric == "mean_per_std":
                    plot_data = ds.sel(RCMs=RCM).mean(dim="GCMs", skipna=True) / ds.sel(
                        RCMs=RCM
                    ).std(dim="GCMs", skipna=True)
                    vmin, vmax = -2, 2
                    cmap = plt.get_cmap("coolwarm")
                plot_data["sfcWind"].plot(
                    x="lon",
                    y="lat",
                    ax=axs[i],
                    vmin=vmin,
                    vmax=vmax,
                    extend="both",
                    cbar_ax=cbar_ax,
                    cmap=cmap,
                    cbar_kwargs={"label": label, "orientation": "horizontal"},
                )  # todo add how many non-nan
                axs[i].set_title(RCM + ", " + str(N_models), fontsize=8)
                axs[i].add_feature(cf.COASTLINE)
                axs[i].add_feature(cf.BORDERS)
                i += 1
                plt.savefig(
                    "../plots/cordex_windchange_"
                    + experiment_id
                    + "_"
                    + metric
                    + ".png",
                    dpi=300,
                    facecolor="w",
                    transparent=False,
                )

"""
GCMs = unique(ds.GCMs)
f, axs = plt.subplots(
    ncols=len(GCMs),
    nrows=2,
    figsize=(10, 4),
    subplot_kw={"projection": ccrs.PlateCarree(), "extent": [-15, 50, 35, 70]},
)
for i, GCM in enumerate(GCMs):
    ds.sel(GCMs=GCM).mean(dim="RCMs", skipna=True)["sfcWind"].plot(
        ax=axs[0, i]
    )  # todo add how many non-nan
    ds.sel(GCMs=GCM).std(dim="RCMs", skipna=True)["sfcWind"].plot(ax=axs[1, i])
"""
