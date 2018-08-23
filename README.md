# SMARTSLAB

This project aims to parameterize mixed-layer depth (MLD) in a way that is efficient for climate models to simulate ocean mixed-layer response in timescale of decades.

## General information

Notes, working progress, references... etc. are in `notes`.
Codes are in `src` with each experiment described as sections below.
The current coupled model that I am using is [GFDL ESM2G](https://www.gfdl.noaa.gov/earth-system-model/). The experiment they run is the historical (145 years) run for [CMIP5](https://cmip.llnl.gov/cmip5/) where the data is downloaded from [ESGF](https://esgf-node.llnl.gov/projects/cmip5/).

## Target Models

- [GFDL ESM2G](https://www.gfdl.noaa.gov/earth-system-model/)
- [EC-Earth](http://www.ec-earth.org/)
- [NCAR CESM](http://www.cesm.ucar.edu/experiments/cesm1.0/)


## 01\_hQ\_fit
This part is to try fitting mixed-layer depth (MLD, or h) and Q-flux (Q). 

Usage:
```bash

julia ./src/01_hQ_fit/case1_h.jl      # Calculate a single value h for the entire ocean. But
                                      # only a number will be displayed.
                                      
julia ./src/01_hQ_fit/case2_h.jl      # Calculate h(x,y) but no time dependency
julia ./src/01_hQ_fit/case2_hQ.jl     # Calculate h(x,y) and Q(x,y) but no time dependency
julia ./src/01_hQ_fit/case3_hQ.jl     # Calculate h(x,y,t) and Q(x,y,t) with t = Jan, Feb, ... Dec. 


```
It should produce NetCDF files in directory `data`.

To plot (Currently the code is in Python3 (matplotlib+numpy):
```bash
python3 ./src/01_hQ_fit/plot_case2_all.py     # Should produce "img/case2_all.png"
python3 ./src/01_hQ_fit/plot_case3.py         # Should produce 12 pictures "img/case3_hQ_month_[month].png"

# python3 ./src/01_hQ_fit/plot_single_point.py [rlat] [rlon] [Description in title]
# [rlat] and [rlon] are rotated latitude and longitude. Currently the transformation between real and rotated coordinates is not clear. Need to ask Professor Keith Moore. This program will look for the nearest grid point to the input coordinate to plot the timeseries of averaged F_TOT, SST and fitted h, Q.

# Example:

python3 ./src/01_hQ_fit/plot_single_point.py   0  -89.5  "EasternPacific_Equator"

# An already written script for selected points
./src/01_hQ_fit/plot_single_points.sh
```

Note that all produced pictures are in `img`.


## 02\_simulation
Code here is to simulate the ML of a single point using the governing equation in [note](https://www.sharelatex.com/read/ffhwmpjxwbht). Some points I choosed are in the file `notes/11-selected_points.md`.

Usage: 
```bash
# julia simulate.jl [lat_index] [lon_index] [Description in title]
#
# Example:

julia ./src/02_simulation/simulate.jl   7  155  "SouthernOcean_Seaice_1 where sea ice exists"
```

To plot:
```bash
# An already written script for selected points
./src/02_simulation/plot_single_points.sh
```
The graphs are output to directory `img`.

