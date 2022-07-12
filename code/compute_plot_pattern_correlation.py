import xarray as xr
import xesmf as xe
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from plot_s10_maps import SUBPLOT_KW, FIG_PARAMS
from scipy.stats import pearsonr, spearmanr
from compute_country_aggregates import add_CMIP5_bounds

SCENARIO_DICT = {"IMAGE": "rcp26", "MINICAM": "rcp45", "MESSAGE": "rcp85"}
EXTENT = SUBPLOT_KW["subplot_kw"]["extent"]


def load_LUH():
    """
    Load LUH datasets into single xarray dataset
    :return:
    """
    luh_list = []
    for luhmodel in ["IMAGE", "MINICAM", "MESSAGE"]:
        ds = xr.open_dataset("../output/LUH/diff_" + luhmodel + ".nc")
        ds = ds.assign_coords({"experiment_id": SCENARIO_DICT[luhmodel]})
        luh_list.append(ds)
    ds = xr.concat(luh_list, dim="experiment_id")
    # Following 3 lines are lat lon metadata additions needed for conservative remapping
    ds.cf.add_bounds(
        ["lon", "lat"]
    )  # add bounds in between grid points (correct here because grid is regular)
    ds.lon.attrs["standard_name"] = "longitude"
    ds.lat.attrs["standard_name"] = "latitude"
    return ds


def load_CMIP5():
    """
    Load CMIP5 datasets into single xarray dataset
    :return:
    """
    CMIP5_list = []
    for experiment_id in ["rcp26", "rcp45", "rcp85"]:
        diff = (
            xr.open_dataset("../output/sfcWind/cmip5_diff_" + experiment_id + ".nc")
            .mean(dim="identifier")
            .squeeze()
        )
        diff = diff.assign_coords({"experiment_id": experiment_id})
        diff = add_CMIP5_bounds(diff)
        CMIP5_list.append(diff)
    ds = xr.concat(CMIP5_list, dim="experiment_id")
    ds = ds.cf.add_bounds(
        ["lon", "lat"]
    )  # add bounds that are needed for conservative remapping
    return ds


def regrid_LUH_onto_CMIP5(ds_LUH, ds_CMIP5, method="conservative"):
    """
    LUH grid is finer than CMIP5. Need to
    :param ds_LUH:
    :param ds_CMIP5:
    :param method: defaults to conservative to conserve area when remapping grid cell fractions from finer LUH to coarser CMIP grid
    :return:
    """
    regridder = xe.Regridder(ds_LUH, ds_CMIP5, method=method)
    return regridder(ds_LUH)


def select_rectangle(ds):
    """
    Select rectangular area defined in EXTENT.
    :param ds: Input dataset with longitudes from 0 to 360 degrees
    :return:
    """
    lon_min, lon_max, lat_min, lat_max = EXTENT
    ds = xr.concat(
        [ds.sel(lon=slice(0, lon_max)), ds.sel(lon=slice(360 + lon_min, 360))],
        dim="lon",
    )  # slicing done in two steps because periodicity at 360° not understood. -15 degrees corresponds to 345 degree here.
    ds = ds.sel(lat=slice(lat_min, lat_max))
    return ds



def convert_to_dataframe(ds_LUH, ds_CMIP5, threshold=0.01):
    """
    Convert xarray datasets to pandas dataframes for plotting

    Only considers those grid boxes where absolute land use change is above the given threshold
    :param ds_LUH:
    :param ds_CMIP5:
    :param threshold: value between 0 (no change) and 1 (maximum change)
    :return:
    """
    # make and populate dataframe
    ds_both = xr.merge([ds_LUH, ds_CMIP5])
    ds_both = ds_both.drop(["lat_bounds", "lon_bounds", "height"])
    df = ds_both.to_dataframe()
    df = df.reset_index()
    df = df[np.abs(df.gothr+df.gsecd) > threshold]  # only consider grid cells with at least 1% change
    return df


def plot_scatter(df, experiment_family):
    # plotting
    R = pearsonr(df["gothr+gsecd"], df["sfcWind"])[0]
    rho = spearmanr(df["gothr+gsecd"], df["sfcWind"])[0]
    g = sns.scatterplot(x="gothr+gsecd", y="sfcWind", hue="experiment_id", data=df, alpha=0.7)
    g.legend_.set_title(None)
    plt.xlabel("Change in primary plus secondary land [fraction of grid cell]")
    plt.ylabel("Wind speed change [m/s]")
    plt.title(experiment_family + ", r=" + str(np.round(R,2))+ ", rho=" + str(np.round(rho,2)))
    plt.tight_layout()
    plt.savefig("../plots/LUH/pattern_correlation.png", **FIG_PARAMS)
    plt.clf()


def make_plot():
    # Load data
    ds_LUH = load_LUH()
    ds_CMIP5 = load_CMIP5()

    # crop and regrid data
    ds_LUH = select_rectangle(
        regrid_LUH_onto_CMIP5(ds_LUH, ds_CMIP5, method="conservative")
    )
    ds_CMIP5 = select_rectangle(ds_CMIP5)


    df = convert_to_dataframe(ds_LUH, ds_CMIP5, 0.01)
    plot_scatter(df, "CMIP5")
