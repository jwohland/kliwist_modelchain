import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import xarray as xr
from numpy import unique, linspace
from pandas import MultiIndex
import xesmf as xe
from plot_utils import *
from scipy.stats import ttest_1samp


FIG_PARAMS = {
    "dpi": 300,
    "facecolor": "w",
    "transparent": False,
}

SUBPLOT_KW = {
    "subplot_kw": {"projection": ccrs.PlateCarree(), "extent": [-15, 50, 35, 70]}
}

TEXT_PARAMS = {
    "horizontalalignment": "left",
    "verticalalignment": "center",
}

MEAN_PLOT_PARAMS = {
    "x": "lon",
    "y": "lat",
    "levels": linspace(-0.7, 0.7, 8),
    "extend": "both",
    "cmap": plt.get_cmap("RdBu_r"),
}


def add_coast_boarders(ax):
    ax.add_feature(cf.COASTLINE)
    ax.add_feature(cf.BORDERS)


def list_of_models(ds, model_type):
    index = -2
    if model_type == "GCM":
        index = 1
    return sorted(list(unique([x.split(".")[index] for x in ds.identifier.values])))


def reindex_per_model(ds):
    tmp = ds.copy()
    RCMs = [x.split(".")[-2] for x in ds.identifier.values]
    GCMs = [x.split(".")[1] for x in ds.identifier.values]
    new_index = MultiIndex.from_arrays([RCMs, GCMs], names=("RCMs", "GCMs"))
    return tmp.assign(identifier=new_index).unstack("identifier")


def plot_array(ds):
    """
    Plot array of change signals in Eurocordex ensemble where each column represents one GCM and each row is one RCM.
    :param ds: xr.Dataset
    :return:
    """
    # prepare plotting
    RCMs, GCMs = list_of_models(ds, "RCM"), list_of_models(ds, "GCM")
    f, axs = plt.subplots(
        ncols=len(GCMs),
        nrows=len(RCMs),
        figsize=(len(RCMs) + 1, len(GCMs) * 1.5),
        **SUBPLOT_KW
    )
    plt.subplots_adjust(0.04, 0.1, 0.95, 0.97, hspace=0.05, wspace=0.05)
    cbar_ax = f.add_axes([0.2, 0.06, 0.6, 0.01])
    label = "Wind speed change 2080-2100 minus 1985-2005 [m/s]"

    for i_ident, ident in enumerate(sorted(ds.identifier.values)):
        GCM, RCM = ident.split(".")[1], ident.split(".")[-2]
        if (
            ident == "EUR-11.CNRM-CERFACS-CNRM-CM5.ICTP.ICTP-RegCM4-6.mon"
            and ds["sfcWind"].sel(identifier=ident).isnull().all()
        ):
            print(ident + " is ignored because nan everywhere")
        else:
            ax = axs[RCMs.index(RCM), GCMs.index(GCM)]
            if i_ident == 0:  # plot with colorbar once
                ds["sfcWind"].sel(identifier=ident).plot(
                    ax=ax,
                    cbar_ax=cbar_ax,
                    cbar_kwargs={"label": label, "orientation": "horizontal"},
                    **MEAN_PLOT_PARAMS
                )
            else:
                ds["sfcWind"].sel(identifier=ident).plot(
                    ax=ax, add_colorbar=False, **MEAN_PLOT_PARAMS
                )
            add_coast_boarders(ax)
            ax.set_title("")
    for i, GCM in enumerate(GCMs):  # add GCM name as column headings
        axs[0, i].set_title(GCM, fontsize=6)
    for j, RCM in enumerate(RCMs):  # add RCM name as row y axis
        axs[j, 0].text(
            -0.1,
            0.5,
            "-".join(RCM.split("-")[1:]),
            rotation="vertical",
            fontsize=5,
            transform=axs[j, 0].transAxes,
            **TEXT_PARAMS
        )


def plot_array_CMIP5(ds):
    """
    Plot array of CMIP5 change signals per model (in different columns) as provided in ds
    :param ds: xr.Dataset
    :return:
    """
    # prepare plotting
    f, axs = plt.subplots(ncols=ds.identifier.size, figsize=(11, 3), **SUBPLOT_KW)
    cbar_ax = f.add_axes([0.2, 0.3, 0.6, 0.05])
    label = "Wind speed change 2080-2100 minus 1985-2005 [m/s]"

    for i, ident in enumerate(sorted(ds.identifier.values)):
        GCM = ident.split(".")[1]
        # plot values
        if i == 0:  # once with colorbar
            ds["sfcWind"].sel(identifier=ident).plot(
                ax=axs.flatten()[i],
                cbar_ax=cbar_ax,
                cbar_kwargs={"label": label, "orientation": "horizontal"},
                **MEAN_PLOT_PARAMS
            )
        else:
            ds["sfcWind"].sel(identifier=ident).plot(
                ax=axs.flatten()[i], add_colorbar=False, **MEAN_PLOT_PARAMS
            )
        axs.flatten()[i].set_title(GCM, fontsize=6)
        add_coast_boarders(axs.flatten()[i])
    plt.subplots_adjust(0.02, 0.15, 0.98, 0.99)


def plot_aggregate(
    ds, model, metric, axs, i, aggregate_dimension, cbar_ax, label, experiment_id
):
    """
    Plot aggregate information by evaluating all RCMs/GCMs available for a specific GCM/RCM in terms of
    - mean change
    - standard deviation of mean change accross the ensemble
    - mean change divided by standard deviation as indicator of signal to noise ratio
    :param ds:
    :param model:
    :param metric:
    :param axs:
    :param i:
    :param aggregate_dimension:
    :return:
    """
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
        cmap = plt.get_cmap("RdBu_r")
        if metric == "mean":
            if aggregate_dimension == "RCM":
                plot_data = ds.sel(RCMs=model).mean(dim="GCMs", skipna=True)
            elif aggregate_dimension == "GCM":
                plot_data = ds.sel(GCMs=model).mean(dim="RCMs", skipna=True)
            levels = [x for x in linspace(-0.5, 0.5, 11) if x != 0]
        elif metric == "standard_deviation":
            if aggregate_dimension == "RCM":
                plot_data = ds.sel(RCMs=model).std(dim="GCMs", skipna=True)
            elif aggregate_dimension == "GCM":
                plot_data = ds.sel(GCMs=model).std(dim="RCMs", skipna=True)
            levels = linspace(0, 0.25, 6)
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
            levels = linspace(-1.9, 1.9, 20)
        plot_data["sfcWind"].plot(
            x="lon",
            y="lat",
            ax=axs[i],
            levels=levels,
            extend="both",
            cbar_ax=cbar_ax,
            cmap=cmap,
            cbar_kwargs={"label": label, "orientation": "horizontal"},
        )
        axs[i].set_title(model + ", " + str(N_models), fontsize=8)
        add_coast_boarders(axs[i])
        i += 1
        plt.savefig(
            "../plots/aggregate/cordex_windchange_"
            + experiment_id
            + "_"
            + metric
            + "_aggdim_"
            + aggregate_dimension
            + ".png",
            **FIG_PARAMS
        )
    return i


def make_individual_plots():
    ### Individual plots of all models ###
    for experiment_id in ["rcp26", "rcp45", "rcp85"]:
        # CORDEX
        diff = xr.open_dataset("../output/sfcWind/cordex_diff_" + experiment_id + ".nc")
        plot_array(diff)
        plt.savefig(
            "../plots/cordex_windchange_" + experiment_id + ".png", **FIG_PARAMS
        )
        # CMIP5
        diff = xr.open_dataset("../output/sfcWind/cmip5_diff_" + experiment_id + ".nc")
        plot_array_CMIP5(diff)
        plt.savefig("../plots/cmip5_windchange_" + experiment_id + ".png", **FIG_PARAMS)
        # todo add historical values. Should be as easy as:
        # ref = xr.open_dataset("../output/sfcWind/cordex_hist_" + experiment_id + ".nc")
        # plot_array(ref)
        # plt.savefig("../plots/cmip5_wind_" + experiment_id + ".png", **FIG_PARAMS)


def make_aggregate_plots():
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
        diff = xr.open_dataset("../output/sfcWind/cordex_diff_" + experiment_id + ".nc")
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
                subplot_params = {
                    "bottom": 0.2,
                    "top": 0.98,
                    "right": 0.96,
                    "left": 0.04,
                }
                cbar_params = [0.1, 0.15, 0.8, 0.04]
            for metric in ["mean", "standard_deviation", "mean_per_std"]:
                # prep figure
                label = metric + " wind speed change 2080-2100 minus 1985-2005 [m/s]"
                f, axs = plt.subplots(
                    nrows=nrows, ncols=ncols, figsize=figsize, **SUBPLOT_KW
                )
                plt.subplots_adjust(**subplot_params)
                cbar_ax = f.add_axes(cbar_params)
                if aggregate_dimension == "RCM":
                    models_of_interest = unique(ds.RCMs)
                elif aggregate_dimension == "GCM":
                    models_of_interest = unique(ds.GCMs)
                i = 0  # can not use enumerate because sometimes loop is over insufficient data
                for model in models_of_interest:
                    i = plot_aggregate(
                        ds,
                        model,
                        metric,
                        axs,
                        i,
                        aggregate_dimension,
                        cbar_ax,
                        label,
                        experiment_id,
                    )

    # CMIP5 and CORDEX average over GCMs and RCMs
    for experiment_family in ["CORDEX", "CMIP5"]:
        for experiment_id in ["rcp26", "rcp45", "rcp85"]:
            diff = xr.open_dataset(
                "../output/sfcWind/"
                + experiment_family.lower()
                + "_diff_"
                + experiment_id
                + ".nc"
            )
            f, axs = plt.subplots(ncols=3, figsize=(12, 4), **SUBPLOT_KW)
            plot_params = {"x": "lon", "y": "lat", "extend": "both"}
            diff.mean(dim="identifier")["sfcWind"].plot(
                ax=axs[0],
                levels=[x for x in linspace(-0.5, 0.5, 11) if x != 0],
                cbar_kwargs={
                    "label": "Wind speed change [m/s]",
                    "orientation": "horizontal",
                },
                **plot_params
            )
            diff.std(dim="identifier")["sfcWind"].plot(
                ax=axs[1],
                vmin=0,
                vmax=0.25,
                cmap=plt.get_cmap("Greens"),
                cbar_kwargs={
                    "label": "Standard deviation of wind speed change [m/s]",
                    "orientation": "horizontal",
                },
                **plot_params
            )
            (diff.mean(dim="identifier") / diff.std(dim="identifier"))["sfcWind"].plot(
                ax=axs[2],
                vmin=-2,
                vmax=2,
                cbar_kwargs={
                    "label": "Mean wind change / standard deviation [1]",
                    "orientation": "horizontal",
                },
                **plot_params
            )
            for i, title in enumerate(
                ["Mean", "Standard Deviation", "Signal to Noise"]
            ):
                axs[i].set_title(title)
                add_coast_boarders(axs[i])
            axs[0].text(
                -0.1,
                1.1,
                experiment_id,
                transform=axs[0].transAxes,
                weight="bold",
                **TEXT_PARAMS
            )
            plt.subplots_adjust(left=0.05, right=0.95, top=0.95)
            plt.savefig(
                "../plots/aggregate/"
                + experiment_family.lower()
                + "_windchange_"
                + experiment_id
                + "_mean.png",
                **FIG_PARAMS
            )


def xarray_ttest(ds):
    """
    Compute two sided T test for data in ds using the null hypothesis that the population mean is zero.
    :param ds: Xarray dataset containing changes in wind speeds with dimensions (r)lat, (r)lon, identifier
    :return:
    """
    p = ttest_1samp(ds["sfcWind"].values, popmean=0).pvalue
    da_p = xr.DataArray(
        p.transpose(),
        coords=ds.mean(dim="identifier").coords,
        dims=ds.mean(dim="identifier").dims,
    )
    return da_p


def make_joint_plots():
    """
    creates plots of the CORDEX and CMIP5 mean changes in all three scenarios and a plot showing the difference in changes
    :return:
    """
    # CORDEX vs. CMIP5 changes
    for cordex_vs_CMIP5_changes in [False, True]:
        plot_params = {"x": "lon", "y": "lat", "extend": "both"}
        if cordex_vs_CMIP5_changes:
            f, axs = plt.subplots(ncols=3, figsize=(12, 4), **SUBPLOT_KW)
        else:
            f, axs = plt.subplots(ncols=3, nrows=2, figsize=(10, 5), **SUBPLOT_KW)
            cbar_ax = f.add_axes([0.15, 0.1, 0.7, 0.05])
        for i, experiment_id in enumerate(["rcp26", "rcp45", "rcp85"]):
            # open data
            ds_cordex = xr.open_dataset(
                "../output/sfcWind/cordex_diff_" + experiment_id + ".nc"
            ).squeeze()
            diff_cordex = ds_cordex.mean(dim="identifier")
            p_cordex = xarray_ttest(ds_cordex)
            ds_cmip5 = xr.open_dataset(
                "../output/sfcWind/cmip5_diff_" + experiment_id + ".nc"
            ).squeeze()
            diff_cmip5 = ds_cmip5.mean(dim="identifier")
            p_cmip5 = xarray_ttest(ds_cmip5)
            if cordex_vs_CMIP5_changes:
                # regrid CORDEX to CMIP5 grid
                regridder = xe.Regridder(
                    diff_cordex, diff_cmip5, "bilinear", periodic=True
                )
                diff_cordex = regridder(diff_cordex).squeeze()
                # plot
                (diff_cordex - diff_cmip5)["sfcWind"].plot(
                    ax=axs[i],
                    levels=[-0.3, -0.2, -0.1, 0.1, 0.2, 0.3],
                    cbar_kwargs={
                        "label": "Wind speed change [m/s]",
                        "orientation": "horizontal",
                    },
                    **plot_params
                )
                axs[i].set_title(experiment_id + "; CORDEX - CMIP5")
                add_coast_boarders(axs[i])
                plt.savefig("../plots/aggregate/diff_windchange_mean.png", **FIG_PARAMS)
            else:
                plt.subplots_adjust(
                    hspace=0.05,
                    wspace=0.05,
                    bottom=0.18,
                    top=0.97,
                    left=0.04,
                    right=0.97,
                )
                for j, ds in enumerate([diff_cmip5, diff_cordex]):
                    if (j == 1) & (i == 1):
                        ds["sfcWind"].plot(
                            ax=axs[j, i],
                            levels=[x for x in linspace(-0.5, 0.5, 11) if x != 0],
                            cbar_ax=cbar_ax,
                            cbar_kwargs={
                                "label": "Ensemble mean wind speed change [m/s]",
                                "orientation": "horizontal",
                            },
                            **plot_params
                        )
                    else:
                        ds["sfcWind"].plot(
                            ax=axs[j, i],
                            levels=[x for x in linspace(-0.5, 0.5, 11) if x != 0],
                            add_colorbar=False,
                            **plot_params
                        )
                    add_coast_boarders(axs[j, i])
                    axs[j, 0].text(
                        -0.1,
                        0.5,
                        ["CMIP5", "EURO-CORDEX"][j],
                        rotation="vertical",
                        fontsize=10,
                        transform=axs[j, 0].transAxes,
                        **TEXT_PARAMS
                    )
                axs[0, i].set_title(experiment_id)
                axs[1, i].set_title("")
                add_letters(axs, x=-0.03, y=1.06, fs=12)
                # add significance hatching, masking areas that are not significant at the 95% level
                p_cmip5.plot.contourf(
                    ax=axs[0, i],
                    levels=[0, 0.05],
                    hatches=["..", ""],
                    alpha=0,
                    **plot_params,
                    add_colorbar=False
                )
                p_cordex.plot.contourf(
                    ax=axs[1, i],
                    levels=[0, 0.05],
                    hatches=["..", ""],
                    alpha=0,
                    **plot_params,
                    add_colorbar=False
                )
                plt.savefig("../plots/aggregate/windchange_mean.png", **FIG_PARAMS)


def make_aggregate_monthly_plots(
    variable_id="sfcWind", method="diff", experiment_ids=["rcp26", "rcp45", "rcp85"]
):
    """
    creates 3x12 subplots where
    - each row is a month
    - each column is a scenario
    Plotted values are changes in the monthly mean future minus historical
    :return:
    """
    units = {
        "sfcWind": "m/s",
        "sic": "fraction of grid cell",
        "tas": "K",
        "ts": "K",
        "tas-ts": "K",
    }
    for experiment_family in ["CORDEX", "CMIP5"]:
        f, axs = plt.subplots(ncols=3, nrows=12, figsize=(6, 18), **SUBPLOT_KW)
        plt.subplots_adjust(0.05, 0.1, 0.97, 0.97, hspace=0.05, wspace=0.05)
        cbar_ax = f.add_axes([0.1, 0.06, 0.8, 0.01])
        plot_params = {"x": "lon", "y": "lat", "extend": "both"}
        if variable_id in ["sfcWind", "sic", "tas-ts"]:
            levels = linspace(-0.9, 0.9, 10)
        elif variable_id in ["tas", "ts"]:
            levels = [x for x in linspace(-5, 5, 11) if x != 0]
        if (variable_id == "tas-ts") & (method == "mean"):
            levels = [x for x in linspace(-2.4, 2.4, 17) if x != 0]
        for i_col, experiment_id in enumerate(experiment_ids):
            diff = xr.open_dataset(
                "../output/"
                + variable_id
                + "/monthly/"
                + experiment_family.lower()
                + "_"
                + method
                + "_"
                + experiment_id
                + ".nc"
            )
            for i_month, month in enumerate(range(1, 13)):
                if i_month + i_col == 0:
                    # plot with colorbar
                    diff.sel(month=month).mean(dim="identifier")[
                        variable_id
                    ].squeeze().plot(
                        ax=axs[i_month, i_col],
                        levels=levels,
                        cbar_kwargs={
                            "label": experiment_family
                            + " "
                            + variable_id
                            + " "
                            + method
                            + " ["
                            + units[variable_id]
                            + "]",
                            "orientation": "horizontal",
                        },
                        cbar_ax=cbar_ax,
                        **plot_params
                    )
                else:
                    # plot without colorbar
                    diff.sel(month=month).mean(dim="identifier")[
                        variable_id
                    ].squeeze().plot(
                        ax=axs[i_month, i_col],
                        levels=levels,
                        add_colorbar=False,
                        **plot_params
                    )

                add_coast_boarders(axs[i_month, i_col])
                axs[i_month, 0].text(
                    -0.15,
                    0.5,
                    str(month),
                    transform=axs[i_month, 0].transAxes,
                    weight="bold",
                    **TEXT_PARAMS
                )
                axs[i_month, i_col].set_title("")
            axs[0, i_col].set_title(experiment_id)

        plt.savefig(
            "../plots/aggregate/monthly/"
            + experiment_family
            + "_"
            + method
            + "_"
            + variable_id
            + "_monthly.png",
            **FIG_PARAMS
        )


def make_s10_maps():
    make_individual_plots()  # wind speed change per GCM and RCM
    make_aggregate_plots()  # wind speed change aggregated over GCMs/RCMs
    make_joint_plots()  # wind speed changes for EURO-CORDEX and CMIP combined
    for variable in [
        "sfcWind",
        "tas",
        "ts",
    ]:
        make_aggregate_monthly_plots(variable)
    try:
        make_aggregate_monthly_plots("sic")
    except:
        print("SIC data incomplete because not available for CMIP5")
    make_aggregate_monthly_plots("tas-ts", "mean", ["historical", "rcp45", "rcp85"])
