import matplotlib.pyplot as plt
import xarray as xr
from plot_s10_maps import add_coast_boarders, SUBPLOT_KW, FIG_PARAMS
from plot_utils import *

PLOT_PARAMS = {
    "levels": [-0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35],
    "extend": "both",
    "cmap": plt.get_cmap("RdBu"),
}


def make_LUH_maps():
    f, axs = plt.subplots(ncols=3, figsize=(10, 3), **SUBPLOT_KW)
    cbar_ax = f.add_axes([0.2, 0.16, 0.6, 0.08])

    for i, luhmodel in enumerate(["IMAGE", "MINICAM", "MESSAGE"]):
        lu_diff = xr.open_dataset("../output/LUH/diff_" + luhmodel + ".nc")
        if i == 0:
            lu_diff["gothr+gsecd"].plot(
                ax=axs[i],
                cbar_ax=cbar_ax,
                cbar_kwargs={
                    "label": "Change in primary plus secondary land [fraction of grid cell]",
                    "orientation": "horizontal",
                },
                **PLOT_PARAMS
            )
        else:
            lu_diff["gothr+gsecd"].plot(ax=axs[i], add_colorbar=False, **PLOT_PARAMS)
        add_coast_boarders(axs[i])
        axs[i].set_title(
            luhmodel + ", rcp" + ["2.6", "4.5", "8.5"][i]
        )  # mapping from IAM to RCP based on Hurtt et al. 2011 (DOI 10.1007/s10584-011-0153-2)Table 1
    add_letters(axs, fs=14, y=1.06, x=-0.02)
    plt.subplots_adjust(left=0.02, right=0.98, top=0.97, bottom=0.2, hspace=0.05, wspace=0.05)
    plt.savefig("../plots/LUH/diff_landuse.png", **FIG_PARAMS)
