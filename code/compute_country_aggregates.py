import geopandas as gpd
import cordex as cx
import xesmf as xe
import xarray as xr


COUNTRIES = [
    "Austria",
    "Belgium",
    "Bulgaria",
    "Croatia",
    "Cyprus",
    "Czechia",
    "Denmark",
    "Estonia",
    "Finland",
    "France",
    "Germany",
    "Greece",
    "Hungary",
    "Ireland",
    "Italy",
    "Latvia",
    "Lithuania",
    "Luxembourg",
    "Malta",
    "Netherlands",
    "Poland",
    "Portugal",
    "Romania",
    "Slovakia",
    "Slovenia",
    "Spain",
    "Sweden",
]  # EU
COUNTRIES.extend(["Norway", "Switzerland", "United Kingdom"])  # non-EU


def prepare_country_geometries():
    """
    Loads shapefile of all countries, subsets those countries defined in
    variable COUNTRIES and outputs their geometries
    :return:
    """
    regs = gpd.read_file(
        "https://cdn.jsdelivr.net/npm/world-atlas@2/countries-10m.json"
    )
    regs["geometry"] = regs.simplify(
        tolerance=0.05, preserve_topology=True
    )  # coarser from 10m resolution to half the resolution of the input data
    regs = regs.set_index("name")  # sort by country name
    return regs.loc[COUNTRIES].geometry


def make_averager_instance():
    """
    Uses clean grid of EUROCORDEX domain with all the metadata to
    instantiate spatial averaging instance
    :return:
    """
    eur11 = cx.cordex_domain("EUR-11", add_vertices=True)
    country_geometries = prepare_country_geometries()
    savg = xe.SpatialAverager(eur11, country_geometries, geom_dim_name="country")
    return savg


def calculate_aggregate(ds):
    savg = make_averager_instance()
    ds_agg = savg(ds)
    ds_agg = ds_agg.assign_coords(country=xr.DataArray(COUNTRIES, dims=("country",)))
    return ds_agg


def compute_aggregates():
    for experiment_family in ["CORDEX", "CMIP5"]:
        for experiment_id in ["rcp85", "rcp45", "rcp26"]:
            for metric in ["mean", "diff"]:
                name = experiment_family.lower() + "_" + metric + "_" + experiment_id
                diff = xr.open_dataset("../output/" + name + ".nc").squeeze()
                diff_agg = calculate_aggregate(diff)
                diff_agg.to_netcdf(
                    "../output/country_aggregates/country_" + name + ".nc"
                )
