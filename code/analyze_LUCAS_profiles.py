import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns
from cordex import preprocessing as preproc

data_path = "/work/ch0636/g300106/projects/kliwist_modelchain/data/LUCAS/"
experiment_dictionary = {"grass": "062008", "forest": "062009", "eval": "062010"}


def open_and_preprocess(experiment, r_lat_center=0, r_lon_center=0):
    ds = xr.open_mfdataset(
        data_path + experiment_dictionary[experiment] + "/*1979*.nc"
    )  # todo generalize beyond 1979
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
    # select region of interest
    ds = ds.sel(
        {
            "rlat": slice(-2 + r_lat_center, 2 + r_lat_center),
            "rlon": slice(-2 + r_lon_center, 2 + r_lon_center),
        }
    )  # todo generalize
    # add wind speed
    ds["S"] = (ds["U"] ** 2 + ds["V"] ** 2) ** (1.0 / 2)
    ds = ds.drop(["U", "V", "rotated_pole"])
    ds = ds.mean(dim=["rlat", "rlon", "time"])
    ds = ds.assign_coords({"experiment": experiment})
    return ds.compute()


for lat_offset in [-10, 0, 10]:
    for lon_offset in [-10, 0, 10]:
        ds_list = [
            open_and_preprocess(x, lat_offset, lon_offset)
            for x in experiment_dictionary.keys()
        ]
        ds_all = xr.concat(ds_list, dim="experiment")

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
            "../plots/LUCAS/scatter_1979_latoffset"
            + str(lat_offset)
            + "_lonoffset_"
            + str(lon_offset)
            + ".png",
            dpi=300,
            facecolor="white"
        )
