import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np

import sys

print("The first argument is %s" % (sys.argv[1],))
method = int(sys.argv[1])

filename = "./data/hQ_method_%d.nc" % method
print("Going to read %s" % filename)

fh = nc4.Dataset(filename, "r")

rlats      = fh.variables['rlat'][:]
rlons      = fh.variables['rlon'][:]


levs_h     = np.arange(-200, 210, 40)
levs_h_std = np.arange(0, 210, 40)

levs_Q     = np.arange(-200, 210, 40)
levs_Q_std = np.arange(0, 210, 40)

cmap_h     = plt.get_cmap("RdBu")
cmap_h_std = plt.get_cmap("hot_r")
cmap_Q     = plt.get_cmap("RdBu")
cmap_Q_std = plt.get_cmap("hot_r")


fig, ax = plt.subplots(2, 2, figsize=(12, 8))
for i in range(12):
    month = i+1
    
    print("Plotting month %02d... " % month, end="")


    fig, ax = plt.subplots(2, 2, figsize=(12, 8))

    ax[0][0].set_title("h best fit")
    ax[0][1].set_title("h std")
    ax[1][0].set_title("Q best fit")
    ax[1][1].set_title("Q std")
  
    fig.suptitle("Month : %02d" % month)
 
    mappable_0 = ax[0][0].contourf(rlons, rlats, fh.variables["h"][i,:,:], levs_h, cmap=cmap_h, extend="both")
    mappable_1 = ax[0][1].contourf(rlons, rlats, fh.variables["h_std"][i,:,:], levs_h_std, cmap=cmap_h_std, extend="max")
    mappable_2 = ax[1][0].contourf(rlons, rlats, fh.variables["Q"][i,:,:], levs_Q, cmap=cmap_Q, extend="both")
    mappable_3 = ax[1][1].contourf(rlons, rlats, fh.variables["Q_std"][i,:,:], levs_Q_std, cmap=cmap_Q_std, extend="max")

    fig.colorbar(mappable_0, ax=ax[0,0], ticks=levs_h)
    fig.colorbar(mappable_1, ax=ax[0,1], ticks=levs_h_std)
    fig.colorbar(mappable_2, ax=ax[1,0], ticks=levs_Q)
    fig.colorbar(mappable_3, ax=ax[1,1], ticks=levs_Q_std)
    

    fig.savefig("img/hQ_method_%d_month_%02d.png" % (method, month), dpi=200)
    print("done.")

