# SMARTSLAB

## Terms definition
- root => top level directory

## Code in 01\_hQ\_fit
This part is to try fitting mixed-layer depth (MLD, or h) and Q-flux (Q). It should produce NC files in $root/data and pictures in $root/img by executing 

```bash
./src/01_hQ_fit/main.sh

```


## Code 02\_simulation
Code here is to simulate the ML of a single point using the governing equation in [note] (https://www.sharelatex.com/read/ffhwmpjxwbht). Some points I choosed are in the file `notes/11-selected_points.md`

Current usage: 
```bash
# simulate.jl [lat_index] [lon_index] [Description in title]
# Example:
julia src/02_simulation/simulate.jl   7  155  "SouthernOcean_Seaice_1 where sea ice exists"

```


