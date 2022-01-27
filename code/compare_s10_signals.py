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
):
    """
    Check data availability and output dictionary of datasets
    :param experiment_family: e.g., CMIP5 or CORDEX
    :param variable_id: e.g., "sfcWind"
    :param frequency: e.g., "monthly"
    :param CORDEX_domain: e.g., "EUR-11"
    :param experiment_id: e.g., "rcp26", "rcp45, "historical"
    :param per_RCM: workaround to build up dictionary stepwise, needed because full loading during historical period fails
    :return:
    """
    link_catalogue = "/pool/data/Catalogs/"  # path to cordex and cmip5 catalog on mistral cluster #todo generalize
    if experiment_family == "CMIP5":
        catalogue = "mistral-cmip5.json"
    elif experiment_family == "CORDEX":
        catalogue = "mistral-cordex.json"
    cat = intake.open_esm_datastore(link_catalogue + catalogue)
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
        ds_dict = subset.to_dataset_dict(cdf_kwargs={"use_cftime": True, "chunks": {}})
    return ds_dict


def preprocess_cordex_dataset(ds, identifier, rotated=True):
    """

    :param ds:
    :param identifier:
    :param rotated: whether model uses rotated coords, True or False
    :return:
    """
    from cordex import preprocessing as preproc

    try:
        ds = preproc.replace_coords(ds)
        # ignore those datasets that have x and y coordinates  #todo there mus be a neater way
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
                    ],
                    errors="ignore",
                )
                    .assign_coords({"identifier": identifier})
                    .groupby("time.year")
                    .mean("time")
            )
            return ds
    except ValueError:
        print("Probably conflicting sizes for rlat dimension. Ignoring " + identifier)




def dictionary_to_dataset(data_dict):
    """
    takes dictionary with keys of the form

    'EUR-11.MOHC-HadGEM2-ES.GERICS.GERICS-REMO2015.rcp45.mon'

    i.e.

    'CORDEX_domain.GCM.RCM_center_RCM.experiment_id.frequency'

    and datasets as values and converts them into a dataset with GCM-RCM as a coordinate
    :param data_dict:
    :return:
    """
    list_ds = [
        preprocess_cordex_dataset(data_dict[identifier], identifier)
        for identifier in data_dict.keys()
    ]
    list_ds = [el for el in list_ds if el]  # remove None
    return xr.concat(list_ds, dim="identifier")


def plot_array(ds):
    # prepare plotting
    f, axs = plt.subplots(
        ncols=4,
        nrows=5,
        subplot_kw={"projection": ccrs.PlateCarree()},
        figsize=(10, 10),
        dpi=300,
    )
    plt.subplots_adjust(bottom=0.15)
    cbar_ax = f.add_axes([0.2, 0.12, 0.6, 0.01])
    label = "wind speed in 2080 - 2100"

    for i, ident in enumerate(ds.identifier.values):
        GCM, RCM = ident.split(".")[1], ident.split(".")[3]
        ds["sfcWind"].sel(identifier=ident).plot(
            ax=axs.flatten()[i],
            cbar_ax=cbar_ax,
            vmin=0,
            vmax=12,
            extend="both",
            cbar_kwargs={"label": label, "orientation": "horizontal"},
        )
        axs.flatten()[i].set_title(GCM + ", " + RCM, fontsize=6)


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
