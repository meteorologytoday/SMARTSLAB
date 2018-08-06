import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np


fh = nc4.Dataset("mixed-layer-depth.nc", "r")

h_spat      = fh.variables['h_spat'][:,:]
h_temp_spat = fh.variables['h_temp_spat'][:,:,:]
rlats      = fh.variables['rlat'][:]
rlons      = fh.variables['rlon'][:]


levs = np.arange(-100, 110, 20)
cmap = plt.get_cmap("RdBu")

fig, ax = plt.subplots(1, 1, figsize=(12, 8))
   
ax.set_title("h_spat")
 
mappable = ax.contourf(rlons, rlats, h_spat, levs, cmap=cmap, extend="both")
fig.colorbar(mappable, ax=ax, ticks=levs)
    

fig.savefig("h_spat.png", dpi=200)



for i in range(12):
    month = i+1

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
   
    ax.set_title("h_temp_spat: %02d" % month)
 
    mappable = ax.contourf(rlons, rlats, h_temp_spat[i], levs, cmap=cmap, extend="both")
    fig.colorbar(mappable, ax=ax, ticks=levs)
    

    fig.savefig("h_temp_spat_%02d.png" % month, dpi=200)


