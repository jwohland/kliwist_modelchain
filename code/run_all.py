# Execute the entire analysis and make all plots
import compute_country_aggregates
import compute_signals
import plot_lu_maps
import plot_s10_maps
import plot_s10_scatter
import plot_s10_country_heatmaps
import analyze_LUCAS_profiles
import plot_temperature_gradient
import compute_plot_pattern_correlation
"""
Calculations
"""
# Compute changes in 10m wind speeds
compute_signals.calculate_signals(time_aggregation="annual", variable_id="sfcWind")
compute_signals.calculate_signals(time_aggregation="monthly", variable_id="sfcWind")
compute_signals.calculate_signals("annual", "sfcWind", full_ensemble=True)
# compute changes in near-surface air temperature (tas), skin temperature (ts) and sea ice cover (sic)
for variable in ["tas", "ts", "sic"]:
    print(variable)
    compute_signals.calculate_signals(time_aggregation="monthly", variable_id=variable)
# compute tas changes for full ensemble
compute_signals.calculate_signals("annual", "tas", full_ensemble=True)
# compute changes in tas - ts, a proxy for stability change
compute_signals.compute_monthly_stability_change()
# Compute country aggregates
compute_country_aggregates.compute_all(onshore=True)  # onshore
compute_country_aggregates.compute_all(onshore=False)  #offshore

"""
Plots
"""
# Maps of 10m wind speeds
plot_s10_maps.make_s10_maps()
# Aggregate plots
plot_s10_maps.make_aggregate_monthly_plots("tas-ts", "mean", ["historical", "rcp45", "rcp85"])
# Scatter plots of wind change per country
plot_s10_scatter.make_s10_scatter(onshore=True)
plot_s10_scatter.make_s10_scatter(onshore=False)
# Heatmaps of wind per country
plot_s10_country_heatmaps.make_s10_heatmaps(onshore=True)
plot_s10_country_heatmaps.make_s10_heatmaps(onshore=False)
# Maps of land use change forcing data
plot_lu_maps.make_LUH_maps()
# Scatter plot of land use forcing vs wind speed change
compute_plot_pattern_correlation.make_plot()
# Temperature gradient plots
plot_temperature_gradient.make_all_plots()

# Plots of REMO LUCAS simulations
analyze_LUCAS_profiles.analyze_LUCAS()