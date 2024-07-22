# Manuscript code for: "Self-organised pattern formation promotes consumer coexistence by fluctuation-dependent mechanisms"
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10233860.svg)](https://doi.org/10.5281/zenodo.10233860)

Authors: Christian Guill, Felix Nößler, and Toni Klauschies

This repository contains scripts to reproduce the results in the manuscript.

## Julia scripts

The [`julia_scripts/run.jl`](julia_scripts/run.jl) file can be used to recreate the figures 3, S6, S7 and S8  of the manuscript. This script calls several smaller scripts to run simulations and create the figures. All information how to run the code can be found inside the header of the [`julia_scripts/run.jl`](julia_scripts/run.jl) file.

Additionally, time series can be generated with the script [`julia_scripts/scripts/timeseries.jl`](julia_scripts/scripts/timeseries.jl).

Main author: Felix Nößler

## Python Jupyter notebook

The [`supporting_figures.ipynb`](supporting_figures.ipynb) contains the Python code to recreate all other figures.

Main author: Christian Guill

--- 

### License

The code in this repository is licensed under the GNU General Public License v3.0 - see the [`COPYING`](COPYING) file in this directory for details.
