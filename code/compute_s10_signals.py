import intake
from dask.distributed import Client
import xarray as xr
import numpy as np

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
            ensemble_member="r1i1p1",
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
            member="r1i1p1",
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

    # set coordinates to standard values
    try:
        ds = preproc.replace_coords(
            ds
        )  # fails of size of rlat/rlon dimension doesn't match
    except ValueError:
        if "MOHC-HadREM3-GA7-05.historical" not in identifier:
            print(
                "Unexpected mismatch in " + identifier
            )  # problem known for the mentioned model and irrelevant for this study because affects only historical runs
        return

    # Remap those datasets that have x and y coordinates
    if "x" in ds.keys():
        ds = preproc.remap_lambert_conformal(ds, domain="EUR-11")
        print(identifier + " does not use rotated coordinates and is remapped")
        ds = ds.drop_dims(["x", "y"])
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
        .groupby("time.year")
        .mean("time")
        .assign_coords({"identifier": identifier})
    )
    return ds


def preprocess_cmip_dataset(ds, identifier):
    ds = (
        ds.drop(["lat_bnds", "lon_bnds", "time_bnds"], errors="ignore")
        .assign_coords({"identifier": identifier})
        .groupby("time.year")
        .mean("time")
    )
    return ds


def dictionary_to_dataset(data_dict, experiment_family):
    """
    takes dictionary with keys of the form

    'EUR-11.MOHC-HadGEM2-ES.GERICS.GERICS-REMO2015.rcp45.mon'

    i.e.

    'CORDEX_domain.GCM.RCM_center_RCM.experiment_id.frequency'

    and datasets as values and converts them into a dataset with GCM-RCM as a coordinate
    :param data_dict:
    :return:
    """
    if experiment_family.lower() == "cordex":
        list_ds = [
            preprocess_cordex_dataset(data_dict[identifier], identifier)
            for identifier in data_dict.keys()
        ]
    else:
        list_ds = [
            preprocess_cmip_dataset(data_dict[identifier], identifier)
            for identifier in data_dict.keys()
        ]
    list_ds = [el for el in list_ds if el]  # remove None
    if experiment_family.lower() == "cmip5":
        # regrid all CMIP5 results to grid of first model
        import xesmf as xe

        for i in range(1, len(list_ds)):
            regridder = xe.Regridder(list_ds[i], list_ds[0], "bilinear")
            list_ds[i] = regridder(list_ds[i])
    return xr.concat(list_ds, dim="identifier")


def update_identifier(ds, experiment_id):
    """
    Create joint identifier that is identical for historical and rcp period to be able to subtract them
    :param ds:
    :return:
    """
    ds["identifier"] = [
        x.replace("." + experiment_id, "") for x in ds.identifier.values
    ]


def calculate_mean(
    experiment_family,
    variable_id,
    frequency,
    CORDEX_domain,
    experiment_id,
    per_RCM=False,
    GCMs=None,
):
    ds_dict = get_dataset_dictionary(
        experiment_family,
        variable_id,
        frequency,
        CORDEX_domain,
        experiment_id,
        per_RCM,
        GCMs=GCMs,
    )
    ds = dictionary_to_dataset(ds_dict, experiment_family)
    if experiment_id == "historical":
        years = slice("1985", "2005")
    else:
        years = slice("2080", "2100")
    ds = ds.sel(year=years).mean("year")
    return ds


for experiment_family in ["CORDEX", "CMIP5"]:
    if experiment_family == "CMIP5":
        GCMs = [
            "CNRM-CM5",
            "EC-EARTH",
            "IPSL-CM5A-MR",
            "HadGEM2-ES",
            "MPI-ESM-LR",
            "NorESM1-M",
        ]
        per_RCM = False
    else:
        GCMs, per_RCM = None, True
    ds_ref = calculate_mean(
        experiment_family, "sfcWind", "mon", "EUR-11", "historical", per_RCM, GCMs
    )
    ds_ref.to_netcdf(
        "../output/" + experiment_family.lower() + "_mean_historical.nc"
    )  # save mean historical
    update_identifier(ds_ref, "historical")
    for experiment_id in ["rcp45"]:  # extend to ["rcp26", "rcp45", "rcp85"]
        ds_future = calculate_mean(
            experiment_family, "sfcWind", "mon", "EUR-11", experiment_id, per_RCM, GCMs
        )
        ds_future.to_netcdf(
            "../output/" + experiment_family.lower() + "_mean_" + experiment_id + ".nc"
        )  # save mean future
        update_identifier(ds_future, experiment_id)

        # calculate and save difference
        diff = ds_future - ds_ref
        diff.to_netcdf(
            "../output/" + experiment_family.lower() + "_diff_" + experiment_id + ".nc"
        )