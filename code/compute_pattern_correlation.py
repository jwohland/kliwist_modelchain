import xarray as xr
import xesmf as xe
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from plot_s10_maps import SUBPLOT_KW
from scipy.stats import pearsonr, spearmanr

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
            xr.open_dataset("../output/cmip5_diff_" + experiment_id + ".nc")
            .mean(dim="identifier")
            .squeeze()
        )
        diff = diff.assign_coords({"experiment_id": experiment_id})
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


# Load data
ds_LUH = load_LUH()
ds_CMIP5 = load_CMIP5()

# crop and regrid data
ds_LUH = select_rectangle(
    regrid_LUH_onto_CMIP5(ds_LUH, ds_CMIP5, method="conservative")
)
ds_CMIP5 = select_rectangle(ds_CMIP5)


f, axs = plt.subplots(ncols=3, sharey=True, figsize=((10, 6)))
xtot, ytot = [], []
for i, experiment_id in enumerate(["rcp26", "rcp45", "rcp85"]):
    ds_LUH_tmp = ds_LUH.sel(experiment_id=experiment_id)["gothr+gsecd"]
    ds_CMIP5_tmp = ds_CMIP5.sel(experiment_id=experiment_id)
    x = ds_LUH_tmp.where(np.abs(ds_LUH_tmp) > 0.0001).values.flatten()
    y = ds_CMIP5_tmp["sfcWind"].where(np.abs(ds_LUH_tmp) > 0.0001).values.flatten()
    mask = np.isfinite(x)
    x = x[mask]
    y = y[mask]
    xtot.extend(x)
    ytot.extend(y)
    R = pearsonr(x, y)[0]
    rho = spearmanr(x, y)[0]
    sns.scatterplot(
        x=x,
        y=y,
        ax=axs[i],
        label="R=" + str(np.round(R, 2)) + "; rho = " + str(np.round(rho, 2)),
    )
    # sns.set_theme(style="darkgrid")
    # sns.jointplot(x=x, y=y, ax=axs[i], kind="reg")
    axs[i].set_title(experiment_id)
axs[0].set_ylabel("Ensemble mean wind speed change [m/s]")
axs[1].set_xlabel("Change in primary plus secondary land [fraction of grid cell]")
plt.savefig("../plots/LUH/pattern_correlation.png")


# plot all 3 scenarios in one
xtot, ytot = np.asarray(xtot), np.asarray(ytot)

for th in [0, 1, 5, 1]:
    x = xtot[np.abs(xtot) > th / 100.0]
    y = ytot[np.abs(xtot) > th / 100.0]
    plt.clf()
    R = pearsonr(x, y)[0]
    rho = spearmanr(x, y)[0]
    # f, ax = plt.subplots()
    g = sns.jointplot(
        x=x,
        y=y,
        label="R=" + str(np.round(R, 2)) + "; rho = " + str(np.round(rho, 2)),
    )
    g.plot_joint(sns.kdeplot, color="r", zorder=10, levels=6)
    plt.title("Threshold = " + str(th))
    plt.xlabel("Change in primary plus secondary land")
    plt.ylabel("Wind speed change [m/s]")
    plt.savefig("../plots/LUH/pattern_correlation_all_th_" + str(th) + ".png", dpi=300)

# todo build up data in a pd Dataframe and make one scatterplot for all scenarios
