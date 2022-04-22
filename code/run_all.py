# Execute the entire analysis and make all plots
import compute_country_aggregates
import plot_lu_maps
import plot_s10_maps
import plot_s10_scatter
import plot_s10_country_heatmaps
"""
Calculations
"""
# Compute country aggregates
compute_country_aggregates.compute_all()

"""
Plots
"""
# Maps of 10m wind speeds
plot_s10_maps.make_s10_maps()
# Scatter plots of wind change per country
plot_s10_scatter.make_s10_scatter()
# Heatmaps of wind per country
plot_s10_country_heatmaps.make_s10_heatmaps()
# Maps of land use change forcing data
plot_lu_maps.make_LUH_maps()
