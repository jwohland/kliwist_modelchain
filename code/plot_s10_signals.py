import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import xarray as xr

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

# CORDEX
diff = xr.open_dataset("../output/cordex_diff_rcp45.nc")
plot_array(diff)
plt.savefig(
    "../plots/RCM-GCM_windchange.png", dpi=300, facecolor="w", transparent=False
)

# CMIP5

diff = xr.open_dataset("../output/cmip5_diff_rcp45.nc")
plot_array_CMIP5(diff)  # todo CMIP5 models use different grids in their output!
plt.savefig("../plots/CMIP5_windchange.png", dpi=300, facecolor="w", transparent=False)
