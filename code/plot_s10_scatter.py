import matplotlib.pyplot as plt

from compute_country_aggregates import COUNTRIES
import xarray as xr

metric = "diff"
for experiment_id in ["rcp26", "rcp45", "rcp85"]:
    name = metric + "_" + experiment_id
    diff_agg_CMIP5 = xr.open_dataset(
        "../output/country_aggregates/country_cmip5_" + name + ".nc"
    )
    diff_agg_CORDEX = xr.open_dataset(
        "../output/country_aggregates/country_cordex_" + name + ".nc"
    )

    for country in COUNTRIES:
        f, ax = plt.subplots(figsize=(4, 4))
        y_CMIP5 = diff_agg_CMIP5.sel({"country": country}).sfcWind
        x_CMIP5 = [0.2] * y_CMIP5.size
        y_CORDEX = diff_agg_CORDEX.sel({"country": country}).sfcWind
        x_CORDEX = [0.8] * y_CORDEX.size
        ax.scatter(x=x_CMIP5, y=y_CMIP5)
        ax.scatter(x=x_CORDEX, y=y_CORDEX)
        ax.set_xticks([0.2, 0.8], ["CMIP5", "EURO-CORDEX"])
        ax.set_ylabel("Wind speed change [m/s] 2080-2100 minus 1985-2005")
        ax.axhline(y=0, ls="--", color="grey")
        ax.set_title(country + " " + experiment_id)
        ax.set_xlim(0, 1)
        plt.tight_layout()
        plt.savefig(
            "../plots/countries/scatter_"
            + metric
            + "_"
            + experiment_id
            + "_"
            + country
            + ".jpeg"
        )
        # todo make nicer with seaborn https://seaborn.pydata.org/examples/scatterplot_categorical.html
        plt.close('all')