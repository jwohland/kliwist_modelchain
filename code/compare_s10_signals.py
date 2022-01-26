import intake
import dask
from dask.distributed import Client
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf

client = Client()  # start dask client


def get_dataset_dictionary(
    experiment_family, variable_id, frequency, CORDEX_domain, experiment_id
):
    """
    Check data availability and output dictionary of datasets
    :param experiment_family: e.g., CMIP5 or CORDEX
    :param variable_id: e.g., "sfcWind"
    :param frequency: e.g., "monthly"
    :param CORDEX_domain: e.g., "EUR-11"
    :param experiment_id: e.g., "rcp26", "rcp45, "historical"
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
    return subset.to_dataset_dict(cdf_kwargs={"use_cftime": True, "chunks": {}})


def preprocess_cordex_dataset(ds, identifier):
    from cordex import preprocessing as preproc

    return preproc.replace_coords(
        ds.drop(
            [
                "time_bnds",
                "rotated_pole",
                "bounds_lon",
                "bounds_lat",
                "time_bounds",
                "Lambert_Conformal",
                "height",
                "x",
                "y"
            ],
            errors="ignore",
        )
        .assign_coords({"identifier": identifier})
        .groupby("time.year")
        .mean("time")
    )


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


cordex_dict = get_dataset_dictionary("CORDEX", "sfcWind", "mon", "EUR-11", "rcp45")

N = len(dset_dict.keys())  # 22


# prepare plotting
f, axs = plt.subplots(
    ncols=4,
    nrows=6,
    subplot_kw={"projection": ccrs.PlateCarree()},
    figsize=(10, 10),
    dpi=300,
)
plt.subplots_adjust(bottom=0.15)
cbar_ax = f.add_axes([0.2, 0.12, 0.6, 0.01])
label = "wind speed in 2080 - 2100"

for i, simulation_name in enumerate(dset_dict.keys()):
    print(simulation_name)
    GCM, RCM = simulation_name.split(".")[1], simulation_name.split(".")[3]
    ds_tmp = dset_dict[simulation_name]
    ds_tmp = ds_tmp.sel(time=slice("2080", "2100"))
    da_tmp = ds_tmp[
        "sfcWind"
    ]  # this had to happen before mean calculation because some variables don't have time dimension I think
    da_tmp = da_tmp.mean(dim="time")
    try:
        da_tmp.plot(
            ax=axs.flatten()[i],
            cbar_ax=cbar_ax,
            vmin=0,
            vmax=12,
            cbar_kwargs={"label": label, "orientation": "horizontal"},
        )
    except:
        # multi realization simulations: pick first ensemble member
        da_tmp = da_tmp.sel(member=da_tmp.member[0])
        da_tmp.plot(
            ax=axs.flatten()[i],
            cbar_ax=cbar_ax,
            vmin=0,
            vmax=12,
            cbar_kwargs={"label": label, "orientation": "horizontal"},
        )
        print("picked first member" + simulation_name)
    axs.flatten()[i].set_title(GCM + ", " + RCM, fontsize=6)
