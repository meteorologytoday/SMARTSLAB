import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np

import sys

filenames = {
    'h'  : "./data/case2_h.jl.nc",
    'hQ' : "./data/case2_hQ.jl.nc",
    'TF' : "./data/dT_star_dt-TOT_F.nc",
}
fh = {}
for key in filenames.keys():
    print("Going to read %s" % (filenames[key],))
    fh[key] = nc4.Dataset(filenames[key], "r")

rlats      = fh['h'].variables['rlat'][:]
rlons      = fh['h'].variables['rlon'][:]

levs = {}
cmap = {}

levs["h"]     = np.arange(-200, 210, 40)
levs["h_std"] = np.arange(0, 210, 40)
levs["Q"]     = np.arange(-500, 500, 50) 
levs["Q_std"] = np.arange(0, 2000, 250)

cmap["h"]     = plt.get_cmap("RdBu")
cmap["h_std"] = plt.get_cmap("jet")
cmap["Q"]     = plt.get_cmap("RdBu_r")
cmap["Q_std"] = plt.get_cmap("jet")


fig, ax = plt.subplots(2, 3, figsize=(18, 12))
fig.subplots_adjust(**{
    'left' : 0.05,
    'right': 0.95
})

fig.suptitle("Case 2", size=30)
ind = (slice(None), slice(None))

# h
mappable_0 = ax[0, 0].contourf(rlons, rlats, fh['h'].variables['h'][ind],     levs['h'],     cmap=cmap['h'],     extend="both")
mappable_1 = ax[1, 0].contourf(rlons, rlats, fh['h'].variables['h_std'][ind], levs['h_std'], cmap=cmap['h_std'], extend="max")

fig.colorbar(mappable_0, ax=ax[0, 0], ticks=levs['h'])
fig.colorbar(mappable_1, ax=ax[1, 0], ticks=levs['h_std'])


# hQ
mappable_0 = ax[0, 1].contourf(rlons, rlats, fh['hQ'].variables['h'][ind],     levs['h'],     cmap=cmap['h'],     extend="both")
mappable_1 = ax[1, 1].contourf(rlons, rlats, fh['hQ'].variables['h_std'][ind], levs['h_std'], cmap=cmap['h_std'], extend="max")
mappable_2 = ax[0, 2].contourf(rlons, rlats, fh['hQ'].variables['Q'][ind],     levs['Q'],     cmap=cmap['Q'],     extend="both")
mappable_3 = ax[1, 2].contourf(rlons, rlats, fh['hQ'].variables['Q_std'][ind], levs['Q_std'], cmap=cmap['Q_std'], extend="max")

fig.colorbar(mappable_0, ax=ax[0, 1], ticks=levs['h'])
fig.colorbar(mappable_1, ax=ax[1, 1], ticks=levs['h_std'])
fig.colorbar(mappable_2, ax=ax[0, 2], ticks=levs['Q'])
fig.colorbar(mappable_3, ax=ax[1, 2], ticks=levs['Q_std'])

# titles

ax[0, 0].set_title(r"$h = h(\lambda, \phi), Q = 0$")
ax[1, 0].set_title(r"$\sigma_h$")
ax[0, 1].set_title(r"$h = h(\lambda, \phi)$")
ax[1, 1].set_title(r"$\sigma_h$")
ax[0, 2].set_title(r"$Q = Q(\lambda, \phi)$")
ax[1, 2].set_title(r"$\sigma_Q$")





fig.savefig("img/case_2.png", dpi=200)
plt.show()

