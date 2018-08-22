# SMARTSLAB

This project aims to parameterize mixed-layer depth (MLD) in a way that is efficient for climate models to simulate ocean mixed-layer response in timescale of decades.


## 01\_hQ\_fit
This part is to try fitting mixed-layer depth (MLD, or h) and Q-flux (Q). 

Usage:
```bash

julia ./src/01_hQ_fit/case1_h.jl      # calculate a single value h for the entire ocean. Only a number will be displayed.
julia ./src/01_hQ_fit/case2_h.jl      # calculate h(x,y) but no time dependency
julia ./src/01_hQ_fit/case2_hQ.jl     # calculate h(x,y) and Q(x,y) but no time dependency
julia ./src/01_hQ_fit/case3_hQ.jl     # calculate h(x,y,t) and Q(x,y,t) with t = Jan, Feb, ... Dec. 


```
It should produce NetCDF files in directory `data`.

To plot (Currently the code is in Python3 (matplotlib+numpy):
```bash
python3 ./src/01_hQ_fit/plot_case2_all.py   # should produce "img/case2_all.png"
python3 src/01_hQ_fit/plot_case3.py         # should produce 12 pictures "img/case3_hQ_month_[month].png"
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

# An already written script for selected points
./src/02_simulation/plot_single_points.sh
```
The graphs are output to directory `img`.

