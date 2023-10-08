# -----------------------------------------------------------
# Script Name: Run all julia scripts
# Purpose: Recreate the figures 2, S2, S6, S7, S8, S9
# Author: Felix Nößler
# Date: 2023-10-08
# -----------------------------------------------------------

# download julia (simulations were tested with julia 1.9.2 and 1.9.3)
# start julia in the subfolder julia_scripts/

# the simulations will be faster by using threads
# you can start julia with 8 threads: julia --threads 8
# see: https://docs.julialang.org/en/v1.9.3/manual/multi-threading/

import Pkg
Pkg.activate(".")
Pkg.instantiate()

# You can rerun the simulations,
# or skip the simulation and just create the plots.
# If you skip the simulation,
# saved files with simulation results will used to create the plots.

# resolution of the simulation
# 300 means:
#  - 300 ⋅ 300 parameter values per heatmap/filled contour plot
#  - 300 parameter values for the bifucation diagram
nvals = 300

################### figure 2: dmaxS vs dmaxI
## run simulation
include("scripts/02_dmax_dmax_env_het.jl")
include("scripts/02_dmax_dmax_pattern.jl")

## create the plot
include("scripts/02_dmax_dmax_plot.jl")
## → the striped patterns were added afterwards with Inkscape

###################  figure S2: bifurcation of a in local system
## run simulation
include("scripts/S2_bifurc_a.jl")

## create the plot
include("scripts/S2_bifurc_a_plot.jl")

################### figure S6: influence of kS=2
## run simulation
include("scripts/S6_kS_kI.jl")

## create the plot
include("scripts/S6_kS_kI_plot.jl")

################### figure S7: influence of kI and aI
## run the simulation
include("scripts/S7_aI_kI.jl")

## create the plot
include("scripts/S7_aI_kI_plot.jl")

################### figure S8: dmaxS vs dmaxI with HS as the invader
## run the simulation
include("scripts/S8_dmax_dmax_pattern_S_invader.jl")

## create the plot
include("scripts/S8_dmax_dmax_S_invader_plot.jl")

################### figure S9: timeseries with HS as the invader
include("scripts/S9_superior_invader_timeseries.jl")
