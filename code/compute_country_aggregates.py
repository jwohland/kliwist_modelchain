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
COUNTRIES.extend(
    ["Norway", "Switzerland", "United Kingdom", "Ukraine", "Turkey"]
)  # non-EU

# Exclude unsuitable countries in offshore assessment
COUNTRIES_EXCLUDE = [
    "Austria",
    "Czechia",
    "Hungary",
    "Luxembourg",
    "Slovakia",
    "Switzerland",
]  # no coast
COUNTRIES_EXCLUDE.extend(
    ["Belgium", "Lithuania", "Slovenia"]
)  # EEZ smaller than 100*100km
OFFSHORE_COUNTRIES = [x for x in COUNTRIES if x not in COUNTRIES_EXCLUDE]


def prepare_country_geometries(onshore=True):
    """
    Loads shapefile of all countries, subsets those countries defined in
    variable COUNTRIES and outputs their geometries
    :return:
    onshore: boolean True if onshore assessment, False if offshore
    countries: list of countries to be evaluated
    """
    if onshore:
        regs = gpd.read_file(
            "https://cdn.jsdelivr.net/npm/world-atlas@2/countries-10m.json"
        )
        regs = regs.set_index("name")  # sort by country name
        countries = COUNTRIES
    else:
        regs = gpd.read_file("../data/EEZ/eez_v11.shp")
        regs = regs.set_index("TERRITORY1")
        regs = regs[regs.GEONAME.str.contains("Exclusive Economic Zone")]  # database also contains five small areas (smaller than 100x100km) that are not exclusive to individual countries
        countries = OFFSHORE_COUNTRIES
    regs["geometry"] = regs.simplify(
        tolerance=0.05, preserve_topology=True
    )  # coarser from 10m resolution to half the resolution of the input data
    return regs.loc[countries].geometry


def make_averager_instance(onshore):
    """
    Uses clean grid of EUROCORDEX domain with all the metadata to
    instantiate spatial averaging instance
    :return:
    """
    eur11 = cx.cordex_domain("EUR-11", add_vertices=True)
    country_geometries = prepare_country_geometries(onshore)
    savg = xe.SpatialAverager(eur11, country_geometries, geom_dim_name="country")
    return savg


def calculate_aggregate(ds, experiment_family, onshore):
    if experiment_family == "CORDEX":
        savg = make_averager_instance(onshore)
    elif experiment_family == "CMIP5":
        savg = xe.SpatialAverager(
            ds, prepare_country_geometries(onshore), geom_dim_name="country"
        )
    ds_agg = savg(ds)
    if onshore:
        countries = COUNTRIES
    else:
        countries = OFFSHORE_COUNTRIES
    ds_agg = ds_agg.assign_coords(country=xr.DataArray(countries, dims=("country",)))
    return ds_agg


def compute_all(onshore):
    for experiment_family in ["CMIP5", "CORDEX"]:
        for experiment_id in ["rcp85", "rcp45", "rcp26", "historical"]:
            for metric in ["mean", "diff"]:
                if metric == "diff" and experiment_id == "historical":
                    print("diffs do not exist for historical by definition")
                else:
                    name = (
                        experiment_family.lower() + "_" + metric + "_" + experiment_id
                    )
                    diff = xr.open_dataset("../output/" + name + ".nc").squeeze()
                    diff_agg = calculate_aggregate(diff, experiment_family, onshore)
                    if onshore:
                        diff_agg.to_netcdf(
                            "../output/country_aggregates/country_" + name + ".nc"
                        )
                    else:
                        diff_agg.to_netcdf(
                            "../output/country_aggregates/offshore/country_" + name + ".nc"
                        )