import intake
import dask
from dask.distributed import Client
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
import cartopy.feature as cf

client = Client()  # start dask client


def get_dataset_dictionary(
    experiment_family,
    variable_id,
    frequency,
    CORDEX_domain,
    experiment_id,
    per_RCM=False,
    GCMs=None,
):
    """
    Check data availability and output dictionary of datasets
    :param experiment_family: e.g., CMIP5 or CORDEX
    :param variable_id: e.g., "sfcWind"
    :param frequency: e.g., "monthly"
    :param CORDEX_domain: e.g., "EUR-11"
    :param experiment_id: e.g., "rcp26", "rcp45, "historical"
    :param per_RCM: workaround to build up dictionary stepwise, needed because full loading during historical period fails
    :param GCMs: list of GCMs that are being searched for. Only implemented for CMIP5
    :return:
    """
    link_catalogue = "/pool/data/Catalogs/"  # path to cordex and cmip5 catalog on mistral cluster #todo generalize
    if experiment_family == "CMIP5":
        catalogue = "mistral-cmip5.json"
    elif experiment_family == "CORDEX":
        catalogue = "mistral-cordex.json"
    cat = intake.open_esm_datastore(link_catalogue + catalogue)
    if experiment_family == "CMIP5":
        subset = cat.search(
            variable=variable_id,
            frequency=frequency,
            experiment=experiment_id,
            ensemble_member="r1i1p1",  # todo this needs to be generalized. Chosen here to avoid duplicates
        )
        if GCMs:  # filter for those GCMs that are of interest
            ds_dict = {}
            for GCM in GCMs:
                try:
                    ds_dict.update(
                        subset.search(model=GCM).to_dataset_dict(
                            cdf_kwargs={"use_cftime": True, "chunks": {}}
                        )
                    )
                except:
                    print(GCM + " not found")
        else:
            ds_dict = subset.to_dataset_dict(
                cdf_kwargs={"use_cftime": True, "chunks": {}}
            )
    elif experiment_family == "CORDEX":
        subset = cat.search(
            variable_id=variable_id,
            frequency=frequency,
            CORDEX_domain=CORDEX_domain,
            experiment_id=experiment_id,
            member="r1i1p1",  # todo this needs to be generalized. Chosen here to avoid duplicates
        )
        if per_RCM:
            import numpy as np

            RCMs = list(np.unique(subset.df.model_id))
            ds_dict = {}
            for RCM in RCMs:
                try:
                    ds_dict.update(
                        subset.search(model_id=RCM).to_dataset_dict(
                            cdf_kwargs={"use_cftime": True, "chunks": {}}
                        )
                    )
                except:
                    print(RCM + " has a problem")
        else:
            ds_dict = subset.to_dataset_dict(
                cdf_kwargs={"use_cftime": True, "chunks": {}}
            )
    return ds_dict


def preprocess_cordex_dataset(ds, identifier):
    """

    :param ds:
    :param identifier:
    :param rotated: whether model uses rotated coords, True or False
    :return:
    """
    from cordex import preprocessing as preproc

    try:
        ds = preproc.replace_coords(ds)
        cont = True
    except ValueError:
        print("Probably conflicting sizes for rlat dimension. Ignoring " + identifier)
        cont = False

    # ignore those datasets that have x and y coordinates  #todo there must be a neater way
    if cont:
        try:
            lol = ds["x"]
            print(identifier + " does not use rotated coordinates and will be ignored")
        except:
            ds = (
                ds.drop(
                    [
                        "time_bnds",
                        "rotated_pole",
                        "bounds_lon",
                        "bounds_lat",
                        "time_bounds",
                        "Lambert_Conformal",
                        "height",
                        "rotated_latitude_longitude",
                        "lat_vertices",
                        "lon_vertices",
                    ],
                    errors="ignore",
                )
                .assign_coords({"identifier": identifier})
                .groupby("time.year")
                .mean("time")
            )
            return ds


def preprocess_cmip_dataset(ds, identifier):
    ds = ds.sel(lat=slice(35, 70), lon=slice(0, 60))
    ds = (
        ds.drop(["lat_bnds", "lon_bnds", "time_bnds"], errors="ignore")
        .assign_coords({"identifier": identifier})
        .groupby("time.year")
        .mean("time")
    )
    return ds


def dictionary_to_dataset(data_dict, cordex=True):
    """
    takes dictionary with keys of the form

    'EUR-11.MOHC-HadGEM2-ES.GERICS.GERICS-REMO2015.rcp45.mon'

    i.e.

    'CORDEX_domain.GCM.RCM_center_RCM.experiment_id.frequency'

    and datasets as values and converts them into a dataset with GCM-RCM as a coordinate
    :param data_dict:
    :return:
    """
    if cordex:
        list_ds = [
            preprocess_cordex_dataset(data_dict[identifier], identifier)
            for identifier in data_dict.keys()
        ]  # todo generalization needed. Maybe preprocess class with boolean cordex variable?
    else:
        list_ds = [
            preprocess_cmip_dataset(data_dict[identifier], identifier)
            for identifier in data_dict.keys()
        ]
    list_ds = [el for el in list_ds if el]  # remove None
    if not cordex:
        # regrid all CMIP5 results to grid of first model
        import xesmf as xe

        for i in range(1, len(list_ds)):
            regridder = xe.Regridder(list_ds[i], list_ds[0], "bilinear")
            list_ds[i] = regridder(list_ds[i])
    return xr.concat(list_ds, dim="identifier")


def plot_array(ds):
    # prepare plotting
    f, axs = plt.subplots(
        ncols=6,
        nrows=5,
        subplot_kw={"projection": ccrs.PlateCarree()},
        figsize=(10, 6),
        dpi=300,
    )
    cbar_ax = f.add_axes([0.2, 0.12, 0.6, 0.01])
    label = "Wind speed change 2080-2100 minus 1985-2005 [m/s]"

    i, j, GCM_old = -1, 0, "No"

    for ident in sorted(ds.identifier.values):
        GCM, RCM = ident.split(".")[1], ident.split(".")[3]
        if GCM != GCM_old:  # move one column to the right for each new GCM
            # turn off axes in plots without data
            [axs[j_e, i].set_axis_off for j_e in range(j + 1, 5)]
            i += 1
            j = 0
            GCM_old = GCM
        ax = axs[j, i]
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

        ax.set_title(GCM + ", " + RCM, fontsize=4)
        j += 1
        ax.add_feature(cf.COASTLINE)
        ax.add_feature(cf.BORDERS)
        ax.set_extent([-15, 50, 35, 70])
    plt.subplots_adjust(0.05, 0.15, 0.95, 0.99)


def plot_array_CMIP5(
    ds,
):  # todo currently copied from plot_array need cleaner plotting functions!
    # prepare plotting
    f, axs = plt.subplots(
        ncols=5,
        nrows=1,
        subplot_kw={"projection": ccrs.PlateCarree()},
        figsize=(8, 2),
        dpi=300,
    )
    cbar_ax = f.add_axes([0.2, 0.12, 0.6, 0.05])
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
    plt.subplots_adjust(0.05, 0.15, 0.95, 0.99)


def update_identifier(ds):
    """
    Create joint identifier that is identical for historical and rcp period to be able to subtract them
    :param ds:
    :return:
    """
    ds["identifier"] = [x.replace(".historical", "") for x in ds.identifier.values]
    ds["identifier"] = [
        x.replace(".rcp45", "") for x in ds.identifier.values
    ]  # todo generalize


cordex_dict_rcp45 = get_dataset_dictionary(
    "CORDEX", "sfcWind", "mon", "EUR-11", "rcp45"
)
cordex_dict_hist = get_dataset_dictionary(
    "CORDEX", "sfcWind", "mon", "EUR-11", "historical", per_RCM=True
)  # doesn't work out of the box. Have to exclude by looping over all RCMs and ignoring those with failing imports

ds_rcp45 = dictionary_to_dataset(cordex_dict_rcp45)
ds_rcp45 = ds_rcp45.sel(year=slice("2080", "2100")).mean("year")

ds_hist = dictionary_to_dataset(
    cordex_dict_hist
)  # some have 413 rlats etc. even though they should have 412
ds_hist = ds_hist.sel(year=slice("1985", "2005")).mean("year")

# calculate difference where it exists, ignore time
update_identifier(ds_hist)
update_identifier(ds_rcp45)
diff = ds_rcp45 - ds_hist

# plotting
plot_array(diff)
plt.savefig(
    "../plots/RCM-GCM_windchange.png", dpi=300, facecolor="w", transparent=False
)

# load CMIP5 models
import numpy as np

GCMs = np.unique(
    [x.split(".")[1] for x in diff.identifier.values]
)  # this gives a combination of institute and model id that is difficult to seperate because "-" is used as separator and as part of the model and institute name
# manual list
GCMs = ["CNRM-CM5", "EC-EARTH", "IPSL-CM5A-MR", "HadGEM2-ES", "MPI-ESM-LR", "NorESM1-M"]

# load RCP45
CMIP5_dict_rcp45 = get_dataset_dictionary(
    "CMIP5", "sfcWind", "mon", "", "rcp45", GCMs=GCMs
)  # todo: NorESM1-M not found
ds_CMIP5_rcp45 = dictionary_to_dataset(CMIP5_dict_rcp45, cordex=False)

# load historical
CMIP5_dict_hist = get_dataset_dictionary(
    "CMIP5", "sfcWind", "mon", "", "historical", GCMs=GCMs
)  # todo: NorESM1-M not found
ds_CMIP5_hist = dictionary_to_dataset(CMIP5_dict_hist, cordex=False)

ds_CMIP5_rcp45 = ds_CMIP5_rcp45.sel(year=slice(2080, 2100)).mean(dim="year")
ds_CMIP5_hist = ds_CMIP5_hist.sel(year=slice(1985, 2005)).mean(dim="year")

update_identifier(ds_CMIP5_hist)
update_identifier(ds_CMIP5_rcp45)


diff = ds_CMIP5_rcp45 - ds_CMIP5_hist

# plotting
plot_array_CMIP5(diff)  # todo CMIP5 models use different grids in their output!
plt.savefig("../plots/CMIP5_windchange.png", dpi=300, facecolor="w", transparent=False)
