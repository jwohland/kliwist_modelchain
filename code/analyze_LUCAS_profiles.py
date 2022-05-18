import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import numpy as np
import cartopy.crs as ccrs

sys.path.append("../code/")
from plot_s10_maps import SUBPLOT_KW, FIG_PARAMS, add_coast_boarders, TEXT_PARAMS

data_path = "/work/ch0636/g300106/projects/kliwist_modelchain/data/LUCAS/"
experiment_dictionary = {"grass": "062008", "forest": "062009", "eval": "062010"}


def open_and_preprocess(experiment):
    """
    Opens 30 year period from 1980 until the end of 2009 of REMO-iMOVE simulations for a given experiment.

    Also performs some preprocessing namely:
        - harmonizing coordinates (winds shifted by half a grid box)
        - computing wind speeds S from components U and V
        - dropping variables that are not needed
        - calculate the temporal mean

    :param experiment: one of "grass", "forest", "eval"
    :return:
    """
    ds = xr.concat(
        xr.open_mfdataset(
            [
                data_path + experiment_dictionary[experiment] + filename
                for filename in ["/*198*.nc", "/*199*.nc", "/*200*.nc"]
            ]
        ),
        dim="time",
    )
    # remove variables that aren't needed for now
    ds = ds.drop(["PS", "hyai", "hybi", "hyam", "hybm", "T"])
    # bring wind variables on same grid. Their grids are shifted by half a cell east or northward
    ds["U"] = (
        ds["U"].assign_coords({"rlon_2": ds.rlon.values}).rename({"rlon_2": "rlon"})
    )
    ds["V"] = (
        ds["V"].assign_coords({"rlat_2": ds.rlat.values}).rename({"rlat_2": "rlat"})
    )
    # bring geopotential on same grid
    ds["FI"] = ds["FI"].rename(
        {"lev_2": "lev"}
    )  # todo check that this actually makes sense
    # drop old grid variables that aren't needed any more
    ds = ds.drop(["rlon_2", "rlat_2", "lev_2"])
    # add wind speed
    ds["S"] = (ds["U"] ** 2 + ds["V"] ** 2) ** (1.0 / 2)
    ds = ds.drop(["U", "V", "rotated_pole"])
    ds = ds.mean(dim=["time"])
    ds = ds.assign_coords({"experiment": experiment})
    return ds


def average_over_region(ds, r_lat_center, r_lon_center):
    # select region of interest
    ds = ds.sel(
        {
            "rlat": slice(-2 + r_lat_center, 2 + r_lat_center),
            "rlon": slice(-2 + r_lon_center, 2 + r_lon_center),
        }
    ).mean(
        dim=["rlat", "rlon"]
    )  # todo generalize
    return ds


def plot_vertical_profile(ds_all):
    """
    # plot vertical wind profile at example locations
    """
    for lat_offset in [-10, 0, 10]:
        for lon_offset in [-10, 0, 10]:
            ds_region = average_over_region(ds_all, lat_offset, lon_offset).compute()
            f, axs = plt.subplots(ncols=2, figsize=((12, 5)))
            df = ds_all.sel(
                {"lev": slice(0, 27)}
            ).to_dataframe()  # convert 8 lowest levels to dataframe for plotting
            sns.scatterplot(data=df, y="S", x="FI", ax=axs[0], hue="experiment")
            axs[0].axvline(x=df["FIB"].values[0], color="grey", ls="--")
            axs[0].set_ylabel(r"Wind speed [m/s]")
            axs[0].set_xlabel("Geopotential [m]")

            df = ds_all.sel(
                {"lev": slice(23, 27)}
            ).to_dataframe()  # convert 5 lowest levels to dataframe for plotting
            # axs[1].set(yscale="log")
            sns.scatterplot(data=df, y="S", x="FI", ax=axs[1], hue="experiment")
            axs[1].axvline(x=df["FIB"].values[0], color="grey", ls="--")
            axs[1].set_ylabel(r"Wind speed [m/s]")
            axs[1].set_xlabel("Geopotential [m]")

            plt.suptitle(
                r"1979 averaged over rlon"
                + str(lon_offset)
                + ", rlat "
                + str(lat_offset)
                + " $\pm$ 2 degree"
            )

            plt.savefig(
                "../plots/LUCAS/scatter_1980-2010_latoffset"
                + str(lat_offset)
                + "_lonoffset_"
                + str(lon_offset)
                + ".png",
                dpi=300,
                facecolor="white",
            )


def plot_mean_change(ds_grass_minus_eval, ds_eval_minus_forest):
    """
    # plot mean change maps at different height levels
    """
    rotated_pole = ccrs.RotatedPole(-162, 39.25)
    # pole location can be obtained from ds.rotated_pole.attrs["grid_north_pole_longitude"], ds.rotated_pole.attrs["grid_north_pole_latitude"]
    f, axs = plt.subplots(
        ncols=2, nrows=5, figsize=((10, 14)), subplot_kw={"projection": rotated_pole}
    )

    for j, lev in enumerate(list(ds_all.lev.values)[-5:]):
        ds_grass_minus_eval["S"].sel({"lev": lev}).plot(
            ax=axs[j, 0],
            # x="lon",
            # y="lat",
            cmap=plt.get_cmap("coolwarm"),
            levels=np.linspace(-2, 2, 9),
            extend="both",
            cbar_kwargs={
                "label": "Wind speed change [m/s]",
                "orientation": "horizontal",
            },
        )
        ds_eval_minus_forest["S"].sel({"lev": lev}).plot(
            ax=axs[j, 1],
            # x="lon",
            # y="lat",
            cmap=plt.get_cmap("coolwarm"),
            levels=np.linspace(-2, 2, 9),
            extend="both",
            cbar_kwargs={
                "label": "Wind speed change [m/s]",
                "orientation": "horizontal",
            },
        )

        for i in [0, 1]:
            add_coast_boarders(axs[j, i])
            axs[j, i].set_title("")
        axs[j, 0].text(
            -0.1,
            0.5,
            "level =" + str(lev),
            rotation="vertical",
            fontsize=10,
            transform=axs[j, 0].transAxes,
            **TEXT_PARAMS
        )
    axs[0, 0].set_title("GRASS - EVAL")
    axs[0, 1].set_title("FOREST - GRASS")

    plt.tight_layout()
    plt.savefig("../plots/LUCAS/REMO_1980-2010_absolute_change_winds_vertical.png", **FIG_PARAMS)


def plot_convergence_maps(ds_grass_minus_eval, ds_eval_minus_forest):
    """
    # plot convergence maps for different absolute and relative convergence criteria
    """
    # todo somewhat better approach but needs double checking whether correct
    rotated_pole = ccrs.RotatedPole(-162, 39.25)
    # pole location can be obtained from ds.rotated_pole.attrs["grid_north_pole_longitude"], ds.rotated_pole.attrs["grid_north_pole_latitude"]
    f, axs = plt.subplots(
        ncols=2, nrows=4, figsize=((6, 10)), subplot_kw={"projection": rotated_pole}
    )

    for i, threshold in enumerate([0.4, 0.2, 0.1, 0.05]):
        for j, ds in enumerate([ds_grass_minus_eval, ds_eval_minus_forest]):
            tmp = ds["S"] < threshold  # binary classifier
            tmp = tmp.cumsum(dim="lev").argmax(
                dim="lev"
            )  # levels are sorted from top of atmosphere to ground level. This locates the lowermost level where change is smaller than threshold
            ds["convergence_map"] = tmp
            ds["convergence_map"].plot(
                levels=np.arange(23, 27),
                cbar_kwargs={"label": "Convergence level"},
                ax=axs[i, j],
                extend="min",
            )
            add_coast_boarders(axs[i, j])
            axs[i, j].set_title("")

        axs[i, 0].text(
            -0.1,
            0.5,
            "threshold =" + str(threshold) + " m/s",
            rotation="vertical",
            fontsize=10,
            transform=axs[i, 0].transAxes,
            **TEXT_PARAMS
        )

        axs[i, 0].set_title("GRASS - EVAL")
        axs[i, 1].set_title("EVAL - FOREST")

    plt.tight_layout()
    plt.savefig(
        "../plots/LUCAS/REMO_1980-2010_convergence_maps.png",
        dpi=300,
        facecolor="white",
    )

def analyze_LUCAS():
    # open and prepare data
    ds_list = [open_and_preprocess(x) for x in experiment_dictionary.keys()]
    ds_all = xr.concat(ds_list, dim="experiment").compute()
    ds_grass_minus_eval = ds_all.sel({"experiment": "grass"}) - ds_all.sel(
        {"experiment": "eval"}
    )
    ds_eval_minus_forest = ds_all.sel({"experiment": "eval"}) - ds_all.sel(
        {"experiment": "forest"}
    )

    # plot vertical wind profiles at example locations
    plot_vertical_profile(ds_all)
    # plot mean change maps at different heights
    plot_mean_change(ds_grass_minus_eval, ds_eval_minus_forest)
    # plot convergence maps
    plot_convergence_maps(ds_grass_minus_eval, ds_eval_minus_forest)
