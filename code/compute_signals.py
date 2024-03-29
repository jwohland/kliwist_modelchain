import warnings

import intake
import xarray as xr
from cordex import preprocessing as preproc


def get_gcm_list(experiment_id):
    """
    return list of GCMs that were downscaled in EUROCORDEX-11 for a given experiment_id.

    This list is manually determined (by evaluating all EU-11 simulations). Automization was difficult because naming
    conventions seem to vary from GCM to GCM.

    :param experiment_id: rcp26, rcp45, rcp85
    :return:
    """
    assert experiment_id in ["rcp26", "rcp45", "rcp85", "all"]
    gcm_list = [
        "CNRM-CM5",
        "EC-EARTH",
        "IPSL-CM5A-MR",
        "HadGEM2-ES",
        "MPI-ESM-LR",
        "NorESM1-M",
    ]  # those models are included in all 3 rcps and the list is complete for rcp45
    if experiment_id == "rcp26":
        gcm_list.extend(
            [
                "MIROC5",
                "GFDL-ESM2G",
            ]
        )
    elif experiment_id == "rcp85":
        gcm_list.extend(
            [
                "CanESM2",
                "MIROC5",
            ]
        )
    elif (
        experiment_id == "all"
    ):  # full list of GCMs that were downscaled in at least 1 scenario
        gcm_list.extend(
            [
                "MIROC5",
                "GFDL-ESM2G",
                "CanESM2",
            ]
        )
    return gcm_list


def get_dataset_dictionary(
    experiment_family,
    variable_id,
    frequency,
    CORDEX_domain,
    experiment_id,
    per_RCM=False,
    GCMs=None,
    standard_ensemble_member="r1i1p1",
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
    link_catalogue = (
        "/pool/data/Catalogs/"  # path to cordex and cmip5 catalog on mistral cluster
    )
    if experiment_family == "CMIP5":
        catalogue = "dkrz_cmip5_disk.json"
    elif experiment_family == "CMIP6":
        catalogue = "dkrz_cmip6_disk.json"
    elif experiment_family == "CORDEX":
        catalogue = "dkrz_cordex_disk.json"
    cat = intake.open_esm_datastore(link_catalogue + catalogue)
    if experiment_family == "CMIP5":
        subset = cat.search(
            variable=variable_id,
            frequency=frequency,
            experiment=experiment_id,
        )
        subset.df.to_csv(
            "../input/CMIP5_" + experiment_id + "_" + variable_id + ".csv"
        )  # dump relevant parts of underlying catalogue
        if GCMs:  # filter for those GCMs that are of interest
            ds_dict = {}
            for GCM in GCMs:
                ensemble_member = standard_ensemble_member
                if GCM == "EC-EARTH":
                    ensemble_member = "r8i1p1"  # EC-Earth doesn't provide all scenarios in realization r1i1p1. Use r8i1p1 instead where all scenarios are provided.
                try:
                    ds_dict.update(
                        subset.search(
                            model=GCM, ensemble_member=ensemble_member
                        ).to_dataset_dict(cdf_kwargs={"use_cftime": True, "chunks": {}})
                    )
                except:
                    print(GCM + " not found")
        else:
            subset = subset.search(ensemble_member=standard_ensemble_member)
            ds_dict = subset.to_dataset_dict(
                cdf_kwargs={"use_cftime": True, "chunks": {}}
            )
    elif experiment_family == "CMIP6":
        subset = cat.search(
            variable_id=variable_id,
            frequency=frequency,
            experiment_id=experiment_id,
            member_id="r1i1p1f1",
            table_id="Amon",
        )
        if experiment_id == "historical":
            subset = subset.search(activity_id="CMIP")
        else:
            subset = subset.search(activity_id="ScenarioMIP")
        subset.df.to_csv(
            "../input/CMIP6_" + experiment_id + "_" + variable_id + ".csv"
        )  # dump relevant parts of underlying catalogue
        ds_dict = subset.to_dataset_dict(cdf_kwargs={"use_cftime": True, "chunks": {}})
    elif experiment_family == "CORDEX":
        subset = cat.search(
            variable_id=variable_id,
            frequency=frequency,
            CORDEX_domain=CORDEX_domain,
            experiment_id=experiment_id,
            member=standard_ensemble_member,
        )
        subset.df.to_csv(
            "../input/CORDEX_" + experiment_id + "_" + variable_id + ".csv"
        )  # dump relevant parts of underlying catalogue
        if per_RCM:
            import numpy as np

            RCMs = list(np.unique(subset.df.model_id))
            ds_dict = {}
            for RCM in RCMs:
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter(
                            "ignore"
                        )  # can be safely ignored here because model will be output
                        ds_dict.update(
                            subset.search(model_id=RCM).to_dataset_dict(
                                preprocess=preproc.rename_cordex,
                                cdf_kwargs={"use_cftime": True, "chunks": {}},
                            )
                        )
                except:
                    print(RCM + " has a problem")
        else:
            ds_dict = subset.to_dataset_dict(
                cdf_kwargs={"use_cftime": True, "chunks": {}}
            )
    return ds_dict


def remap_hadrem3(ds, CORDEX_domain="EUR-11", method="bilinear"):
    """
    Remaps MOHC-HadREM3-GA7-05 onto the normal rotated EUR-11 grid.
    The model uses 413 instead of 412 rlats which are shifted by about
    one half grid box
    :param ds:
    :param CORDEX_domain: name of the CORDEX domain, e.g., "EUR-11"
    :param method: method for interpolation, default "bilinear"
    :return:
    """
    from cordex import core as core
    import xesmf as xe

    ds_cordex = core.domain.cordex_domain(CORDEX_domain)
    regridder = xe.Regridder(ds, ds_cordex, method=method)
    ds = regridder(ds)
    return preproc.replace_coords(ds, "EUR-11")


def preprocess_cordex_dataset(ds, identifier):
    """

    :param ds:
    :param identifier:
    :param rotated: whether model uses rotated coords, True or False
    :return:
    """
    from cordex import preprocessing as preproc

    # one WRF downscaling introduces grid related problems, only for rcp85, so is manually excluded here
    # running preproc.rename_cordex didn't solve the issue (all values E+33 afterwards for this particular simulation)
    if identifier == "EUR-11.MIROC-MIROC5.UHOH.UHOH-WRF361H.rcp85.mon":
        return

    # fix HadREM3 grid (too large by 1 grid box)
    if "MOHC-HadREM3-GA7-05" in identifier:
        ds = remap_hadrem3(ds)
    # Remap those datasets that have x and y coordinates
    elif "x" in ds.keys():
        ds = preproc.remap_lambert_conformal(ds, domain="EUR-11")
        print(identifier + " does not use rotated coordinates and is remapped")
        ds = ds.drop_dims(["x", "y"], errors="ignore")
    elif (ds["rlat"].size == 412) & (ds["rlon"].size == 424):
        ds = preproc.replace_coords(ds)
    else:
        print("Grid not understood " + identifier)
        return None
    ds = ds.drop(
        [
            "time_bnds",
            "bnds",
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
    ).assign_coords({"identifier": identifier})
    return ds


def preprocess_cmip_dataset(ds, identifier):
    ds = ds.drop(
        ["lat_bnds", "lon_bnds", "time_bnds", "bnds"], errors="ignore"
    ).assign_coords({"identifier": identifier})
    try:  # CMIP5
        ds = ds.sel(ensemble_member=ds.ensemble_member, drop=True).squeeze()
        ds = ds.assign_attrs(
            ensemble_information="normally r1i1p1. Exception EC-EARTH r7i1p1"
        )
    except:  # CMIP6
        ds = ds.sel(member_id=ds.member_id, drop=True).squeeze()
        ds = ds.assign_attrs(ensemble_information="r1i1p1f1")
    return ds


def aggregate_temporally(ds, experiment_id, time_aggregation):
    """
    Aggregates the time dimension from input resolution to resolution specified in time_aggretation.
    Also selects the period based on experiment_id
    :param ds:
    :param experiment_id:
    :param time_aggregation:
    :return:
    """
    assert time_aggregation in [
        "annual",
        "monthly",
    ]  # monthly and annual resampling supported
    if experiment_id == "historical":
        years = slice("1985", "2005")
    else:
        years = slice("2080", "2100")
    if time_aggregation == "annual":
        try:
            ds = ds.sel(time=years).mean("time")
        except:  # some datasets do not cover the correct period
            ds = None
            """"""
    elif time_aggregation == "monthly":
        ds = ds.sel(time=years).groupby("time.month").mean("time")
    return ds


def dictionary_to_dataset(
    data_dict, experiment_family, experiment_id, time_aggregation
):
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
    list_ds = [ds for ds in list_ds if ds]  # remove None
    # aggregate temporally
    list_ds = [
        aggregate_temporally(ds, experiment_id, time_aggregation) for ds in list_ds
    ]
    list_ds = [ds for ds in list_ds if ds]  # remove None
    if experiment_family.lower() in ["cmip5", "cmip6"]:
        # regrid all CMIP5 results to a fixed grid
        import xesmf as xe
        import numpy as np

        ds_typical = xr.Dataset(
            {
                "lat": (["lat"], np.arange(-89.25, 90, 1.5)),
                "lon": (["lon"], np.arange(0, 360, 1.75)),
            }
        )  # this resolution is representative for CMIP5 horizontal resolutions and avoids grid points at the pole
        for i in range(0, len(list_ds)):
            try:
                regridder = xe.Regridder(
                    list_ds[i],
                    ds_typical,
                    "bilinear",
                    periodic=True,
                    ignore_degenerate=True,
                )
                list_ds[i] = regridder(list_ds[i])
            except:
                list_ds[i] = None
                if experiment_family == "CMIP6":
                    print("One simulation can not be regridded")  # this occurs once
                else:
                    print("Unexpected error")
        list_ds = [x.drop("height", errors="ignore") for x in list_ds if x]
    return xr.concat(list_ds, dim="identifier")


def update_identifier(ds, experiment_id, experiment_family="CMIP5"):
    """
    Create joint identifier that is identical for historical and rcp period to be able to subtract them
    :param ds:
    :return:
    """
    if experiment_family in ["CMIP5", "CORDEX"]:
        ds["identifier"] = [
            x.replace("." + experiment_id, "") for x in ds.identifier.values
        ]
    elif experiment_family == "CMIP6":
        ds["identifier"] = [x.split(".")[1] for x in ds.identifier.values]


def calculate_mean(
    experiment_family,
    variable_id,
    frequency,
    CORDEX_domain,
    experiment_id,
    per_RCM=False,
    GCMs=None,
    time_aggregation="annual",
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
    # manual fixes for some datasets
    # ICHEC-EC-EARTH throws a key Error for tas. Remove manually for now.
    if variable_id == "tas":
        for ignore_model in [
            "ICHEC.EC-EARTH.historical.Amon",
            "LASG-CESS.FGOALS-g2.historical.Amon",
            "NCAR.CCSM4.rcp26.Amon",
            "NCAR.CCSM4.rcp45.Amon",
            "NOAA-GFDL.GFDL-CM2p1.rcp45.Amon",
            "ICHEC.EC-EARTH.rcp45.Amon",
        ]:
            try:
                del ds_dict[ignore_model]
            except:
                """"""

    # in the large CMIP5 ensemble some tas are reported at 1.5m and some at 2m
    if variable_id == "tas":
        for model in list(ds_dict):
            ds_dict[model] = ds_dict[model].drop("height", errors="ignore")

    ds = dictionary_to_dataset(
        ds_dict, experiment_family, experiment_id, time_aggregation
    )
    return ds


def calculate_signals(
    time_aggregation="annual", variable_id="sfcWind", full_ensemble=False
):
    """
    Calculates 20y mean wind speed at the end of the historical and future periods.

    :param time_aggregation: "annual" or "monthly": Whether annual means or monthly means are requested
    :param variable_id: name of the variable, for example, "sfcWind" or "tas"
    :param full_ensemble: if false, only use those models that were downscaled in EURO-CORDEX. If true, use full ensemble
    :return:
    """
    out_path = "../output/" + variable_id + "/"
    if time_aggregation == "monthly":
        out_path += "monthly/"
    if full_ensemble:
        experiment_families = ["CMIP5", "CMIP6"]
        out_path += "all/"
    else:
        experiment_families = ["CORDEX", "CMIP5"]
    for experiment_family in experiment_families:
        if full_ensemble:
            GCMs, per_RCM = None, False
        elif (
            experiment_family == "CMIP5"
        ):  # only those CMIP5 models that were downscaled in EURO-CORDEX
            GCMs = get_gcm_list("all")
            per_RCM = False
        else:
            GCMs, per_RCM = None, True  # EURO-CORDEX

        ds_ref = calculate_mean(
            experiment_family,
            variable_id,
            "mon",
            "EUR-11",
            "historical",
            per_RCM,
            GCMs,
            time_aggregation,
        )
        ds_ref.to_netcdf(
            out_path + experiment_family.lower() + "_mean_historical.nc"
        )  # save mean historical
        update_identifier(ds_ref, "historical", experiment_family=experiment_family)
        if experiment_family == "CMIP6":
            experiment_ids = ["ssp126", "ssp245", "ssp370", "ssp585"]
        else:
            experiment_ids = ["rcp85", "rcp45", "rcp26"]

        for experiment_id in experiment_ids:
            if experiment_family == "CMIP5" and not full_ensemble:
                GCMs = get_gcm_list(experiment_id)
            ds_future = calculate_mean(
                experiment_family,
                variable_id,
                "mon",
                "EUR-11",
                experiment_id,
                per_RCM,
                GCMs,
                time_aggregation,
            )
            ds_future.to_netcdf(
                out_path + experiment_family.lower() + "_mean_" + experiment_id + ".nc"
            )  # save mean future
            update_identifier(ds_future, experiment_id, experiment_family)

            # calculate and save difference
            diff = ds_future - ds_ref
            diff.to_netcdf(
                out_path + experiment_family.lower() + "_diff_" + experiment_id + ".nc"
            )


def compute_monthly_stability_change():
    for experiment_family in ["CORDEX", "CMIP5"]:
        for method in ["mean", "diff"]:
            experiment_ids = ["rcp85", "rcp45", "rcp26"]
            if method == "mean":
                experiment_ids.append("historical")
            for experiment_id in experiment_ids:
                ds_tas = xr.open_dataset(
                    "../output/tas/monthly/"
                    + experiment_family.lower()
                    + "_"
                    + method
                    + "_"
                    + experiment_id
                    + ".nc"
                )
                ds_ts = xr.open_dataset(
                    "../output/ts/monthly/"
                    + experiment_family.lower()
                    + "_"
                    + method
                    + "_"
                    + experiment_id
                    + ".nc"
                )
                ds_stability = (
                    (ds_tas["tas"] - ds_ts["ts"])
                    .drop(["identifier", "member"], errors="ignore")
                    .to_dataset(name="tas-ts")
                    .squeeze()
                )
                ds_stability.to_netcdf(
                    "../output/tas-ts/monthly/"
                    + experiment_family.lower()
                    + "_"
                    + method
                    + "_"
                    + experiment_id
                    + ".nc"
                )
