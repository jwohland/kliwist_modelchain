import xarray as xr
import pandas as pd


def open_LUH(filename, year=None, name=None):
    """Open LUH1 (land-use change) raster and add data on year and name (if given)"""
    da = (
        xr.open_rasterio(filename)
        .rename({"x": "lon", "y": "lat"})
        .drop("band")
        .squeeze()
    )
    if year:
        da = da.assign_coords({"time": pd.Timestamp(str(year) + "-01-01")})
    if name:
        da.name = name
    return da


def open_LUH_period(luhmodel, ts, te):
    """
    open time varying primary and secondary vegeation maps and combine them in a single xarray Dataset
    :param path: path to data
    :param ts: start year of respective experiment, e.g. 1850 for historical
    :param te: end year for respective experiment, e.g. 2000 for historical
    :return:
    """
    path_to_data = f"../data/LUH/{luhmodel}/updated_states/{{indicator}}.{{year}}.txt"

    ds = xr.merge(
        [
            xr.concat(
                [
                    open_LUH(
                        path_to_data.format(indicator=indicator, year=year),
                        year,
                        indicator,
                    )
                    for year in range(ts, te)
                ],
                dim="time",
            )
            for indicator in ["gothr", "gsecd"]
        ]
    )
    ds["gothr+gsecd"] = ds["gothr"] + ds["gsecd"]
    return ds

def compute_lu_differences():
    ds_historic = open_LUH_period("historic", 1985, 2005).mean(dim="time")

    for luhmodel in ["IMAGE", "MESSAGE", "MINICAM"]:
        ds_future = open_LUH_period(luhmodel, 2080, 2100).mean(dim="time")
        lu_diff = ds_future - ds_historic
        lu_diff.to_netcdf("../output/LUH/diff_" + luhmodel + ".nc")
