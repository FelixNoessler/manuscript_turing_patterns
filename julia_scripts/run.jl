# -----------------------------------------------------------
# Script Name: Run all julia scripts
# Purpose: Recreate the figures 2, S2, S6, S7, S8, S10, S11
# Author: Felix Nößler
# Date: 2023-11-30
# -----------------------------------------------------------

# download julia (simulations were tested with julia 1.9.2 and 1.9.3)
# start julia in the parent folder or in the subfolder julia_scripts/

# Reproducibility of the code and the same version number of the
# Julia packages is facilitated by the Manifest.toml.

# the simulations will be faster by using threads
# you can start julia with 8 threads: julia --threads 8
# see: https://docs.julialang.org/en/v1.9.3/manual/multi-threading/

import Pkg

if isdir("julia_scripts")
    cd("julia_scripts")
elseif "julia_scripts" != basename(pwd())
    error("Please start julia either in the parent folder or in julia_scripts")
end

Pkg.activate(".")
Pkg.instantiate()

# You can rerun the simulations (set only_plot to false),
# or skip the simulation and just create the plots.
# If you skip the simulation,
# saved files with simulation results will used to create the plots.
only_plot = false

# resolution of the simulation
# 300 means:
#  - 300 ⋅ 300 parameter values per heatmap/filled contour plot
#  - 300 parameter values for the bifucation diagram
nvals = 300

################### Figure 3
if !only_plot
    ## run simulation
    include("scripts/03_dS_dI_env_het.jl")
    include("scripts/03_dS_dI_pattern.jl")
end

## create the plot
include("scripts/03_dS_dI_plot.jl")
## → the striped patterns were added afterwards with Inkscape

###################  Figure S6
if !only_plot
    ## run simulation
    include("scripts/S6_aI_dI_panelA.jl")
    include("scripts/S6_aI_dI_panelB.jl")
end

## create the plot
include("scripts/S6_aI_dI_plot.jl")

################### Figure S7
if !only_plot
    ## run simulation
    include("scripts/S7_I_plastic_dispersal_env_het.jl")
    include("scripts/S7_I_plastic_dispersal_pattern.jl")
end

## create the plot
include("scripts/S7_I_plastic_dispersal_plot.jl")

################### Figure S8
if !only_plot
    ## run simulation
    include("scripts/S8_S_invader.jl")
end

## create the plot
include("scripts/S8_S_invader_plot.jl")
