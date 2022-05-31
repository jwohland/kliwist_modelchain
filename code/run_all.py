# Execute the entire analysis and make all plots
import compute_country_aggregates
import compute_s10_signals
import plot_lu_maps
import plot_s10_maps
import plot_s10_scatter
import plot_s10_country_heatmaps
import analyze_LUCAS_profiles
"""
Calculations
"""
# Compute changes in 10m wind speeds
compute_s10_signals.calculate_signals(time_aggregation="annual")
compute_s10_signals.calculate_signals(time_aggregation="monthly")
# Compute country aggregates
compute_country_aggregates.compute_all(onshore=True)  # onshore
compute_country_aggregates.compute_all(onshore=False)  #offshore

"""
Plots
"""
# Maps of 10m wind speeds
plot_s10_maps.make_s10_maps()
# Scatter plots of wind change per country
plot_s10_scatter.make_s10_scatter(onshore=True)
plot_s10_scatter.make_s10_scatter(onshore=False)
# Heatmaps of wind per country
plot_s10_country_heatmaps.make_s10_heatmaps(onshore=True)
plot_s10_country_heatmaps.make_s10_heatmaps(onshore=False)
# Maps of land use change forcing data
plot_lu_maps.make_LUH_maps()


# Plots of REMO LUCAS simulations
analyze_LUCAS_profiles.analyze_LUCAS()