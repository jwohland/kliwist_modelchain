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


def plot_aggregate(ds, model, metric, axs, i, aggregate_dimension):
    """
    Plot aggregate information by evaluating all RCMs/GCMs available for a specific GCM/RCM in terms of
    - mean change
    - standard deviation of mean change accross the ensemble
    - mean change divided by standard deviation as indicator of signal to noise ratio  # todo think about 1/N scaling or so
    :param ds:
    :param model:
    :param metric:
    :param axs:
    :param i:
    :param aggregate_dimension:
    :return:
    """
    # todo there must be a nicer way of doing this with fewer if and elifs
    if aggregate_dimension == "RCM":
        N_models = (
            len(ds.GCMs)
            - ds.sel(RCMs=model)
            .mean(dim=["rlat", "rlon"])["sfcWind"]
            .isnull()
            .values.sum()
        )
    elif aggregate_dimension == "GCM":
        N_models = (
            len(ds.RCMs)
            - ds.sel(GCMs=model)
            .mean(dim=["rlat", "rlon"])["sfcWind"]
            .isnull()
            .values.sum()
        )
    if N_models > 1:
        cmap = plt.get_cmap("coolwarm")
        if metric == "mean":
            if aggregate_dimension == "RCM":
                plot_data = ds.sel(RCMs=model).mean(dim="GCMs", skipna=True)
            elif aggregate_dimension == "GCM":
                plot_data = ds.sel(GCMs=model).mean(dim="RCMs", skipna=True)
            vmin, vmax = -0.5, 0.5
        elif metric == "standard_deviation":
            if aggregate_dimension == "RCM":
                plot_data = ds.sel(RCMs=model).std(dim="GCMs", skipna=True)
            elif aggregate_dimension == "GCM":
                plot_data = ds.sel(GCMs=model).std(dim="RCMs", skipna=True)
            vmin, vmax = 0, 0.25
            cmap = plt.get_cmap("Greens")
        elif metric == "mean_per_std":
            if aggregate_dimension == "RCM":
                plot_data = ds.sel(RCMs=model).mean(dim="GCMs", skipna=True) / ds.sel(
                    RCMs=model
                ).std(dim="GCMs", skipna=True)
            elif aggregate_dimension == "GCM":
                plot_data = ds.sel(GCMs=model).mean(dim="RCMs", skipna=True) / ds.sel(
                    GCMs=model
                ).std(dim="RCMs", skipna=True)
            vmin, vmax = -2, 2
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
        )
        if metric == "mean":
            # add hatching  # todo maybe make hatching differently, e.g., by masking out areas with low signal to noise ratio
            plot_data["sfcWind"].plot(
                ax=axs[i],
                x="lon",
                y="lat",
                levels=[-10, -0.1, 0.1, 10],
                colors=["none", "white", "none"],
                add_colorbar=False,
            )
        axs[i].set_title(model + ", " + str(N_models), fontsize=8)
        axs[i].add_feature(cf.COASTLINE)
        axs[i].add_feature(cf.BORDERS)
        i += 1
        plt.savefig(
            "../plots/aggregate/cordex_windchange_"
            + experiment_id
            + "_"
            + metric
            + "_aggdim_"
            + aggregate_dimension
            + ".png",
            dpi=300,
            facecolor="w",
            transparent=False,
        )
    return i


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
ensemble_size = {
    "RCM": {
        "rcp26": 4,
        "rcp45": 5,
        "rcp85": 10,
    },  # number of RCMs with >1 downscaled GCM
    "GCM": {
        "rcp26": 5,
        "rcp45": 5,
        "rcp85": 8,
    },  # number of GCM with >1 downscaling RCM
}
for experiment_id in ["rcp26", "rcp45", "rcp85"]:
    diff = xr.open_dataset("../output/cordex_diff_" + experiment_id + ".nc")
    ds = reindex_per_model(diff)
    for aggregate_dimension in ["RCM", "GCM"]:
        if aggregate_dimension == "RCM":  # RCMs per row as in matrix plot
            nrows, ncols = ensemble_size[aggregate_dimension][experiment_id], 1
            figsize = (4, len(experiment_id) * 2)
            subplot_params = {"bottom": 0.08, "top": 0.98}
            cbar_params = [0.1, 0.05, 0.8, 0.01]
        elif aggregate_dimension == "GCM":  # GCMs per column
            nrows, ncols = 1, ensemble_size[aggregate_dimension][experiment_id]
            figsize = (len(experiment_id) * 4, 3.5)
            subplot_params = {"bottom": 0.2, "top": 0.98, "right": 0.96, "left": 0.04}
            cbar_params = [0.1, 0.15, 0.8, 0.04]
        for metric in ["mean", "standard_deviation", "mean_per_std"]:
            # prep figure
            label = metric + " wind speed change 2080-2100 minus 1985-2005 [m/s]"
            f, axs = plt.subplots(
                nrows=nrows,
                ncols=ncols,
                subplot_kw={
                    "projection": ccrs.PlateCarree(),
                    "extent": [-15, 50, 35, 70],
                },
                figsize=figsize,
            )
            plt.subplots_adjust(**subplot_params)
            cbar_ax = f.add_axes(cbar_params)
            if aggregate_dimension == "RCM":
                models_of_interest = unique(ds.RCMs)
            elif aggregate_dimension == "GCM":
                models_of_interest = unique(ds.GCMs)
            i = 0  # can not use enumerate because sometimes loop is over insufficient data
            for model in models_of_interest:
                i = plot_aggregate(ds, model, metric, axs, i, aggregate_dimension)

# CMIP5
for experiment_id in ["rcp26", "rcp45", "rcp85"]:
    # CMIP5
    diff = xr.open_dataset("../output/cmip5_diff_" + experiment_id + ".nc")
    f, axs = plt.subplots(
        ncols=3,
        subplot_kw={"projection": ccrs.PlateCarree(), "extent": [-15, 50, 35, 70]},
        figsize=(12, 4),
        dpi=300,
    )
    diff.mean(dim="identifier")["sfcWind"].plot(
        ax=axs[0],
        vmin=-0.5,
        vmax=0.5,
        cmap=plt.get_cmap("coolwarm"),
        extend="both",
        cbar_kwargs={"label": "Wind speed change [m/s]", "orientation": "horizontal"},
    )
    diff.mean(dim="identifier")["sfcWind"].plot(
        ax=axs[0],
        levels=[-10, -0.1, 0.1, 10],
        colors=["none", "white", "none"],
        add_colorbar=False,
    )
    diff.std(dim="identifier")["sfcWind"].plot(
        ax=axs[1],
        vmin=0,
        vmax=0.25,
        cmap=plt.get_cmap("Greens"),
        extend="both",
        cbar_kwargs={
            "label": "Standard deviation of wind speed change [m/s]",
            "orientation": "horizontal",
        },
    )
    (diff.mean(dim="identifier") / diff.std(dim="identifier"))["sfcWind"].plot(
        ax=axs[2],
        vmin=-2,
        vmax=2,
        extend="both",
        cmap=plt.get_cmap("coolwarm"),
        cbar_kwargs={
            "label": "Mean wind change / standard deviation [1]",
            "orientation": "horizontal",
        },
    )
    for i, title in enumerate(["Mean", "Standard Deviation", "Signal to Noise"]):
        axs[i].set_title(title)
        axs[i].add_feature(cf.COASTLINE)
        axs[i].add_feature(cf.BORDERS)
    axs[0].text(
        -0.1,
        1.1,
        experiment_id,
        horizontalalignment="left",
        verticalalignment="center",
        transform=axs[0].transAxes,
        weight="bold",
    )
    plt.savefig(
        "../plots/aggregate/cmip5_windchange_" + experiment_id + "_mean.png",
        dpi=300,
        facecolor="w",
        transparent=False,
    )
