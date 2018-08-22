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

levs["h"]     = np.arange(-100, 110, 10)
levs["h_std"] = np.arange(0, 210, 40)
levs["Q"]     = np.arange(-100, 110, 10) 
levs["Q_std"] = np.arange(0, 2000, 250)
levs["TOT_F"]     = np.arange(-100, 110, 10) 

cmap["h"]     = plt.get_cmap("RdBu_r")
cmap["h_std"] = plt.get_cmap("jet")
cmap["Q"]     = plt.get_cmap("RdBu_r")
cmap["Q_std"] = plt.get_cmap("jet")
cmap["TOT_F"] = plt.get_cmap("RdBu_r")


fig, ax = plt.subplots(2, 4, figsize=(24, 8))
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

# TOT_F
mappable_0 = ax[0, 3].contourf(rlons, rlats, np.mean(fh['TF'].variables['TOT_F'][:,:,:], axis=0), levs['TOT_F'], cmap=cmap['TOT_F'], extend="both")
fig.colorbar(mappable_0, ax=ax[0, 3], ticks=levs['TOT_F'])

ax[1, 3].axis('off')

# titles

ax[0, 0].set_title(r"$h = h(\lambda, \phi)$ [$\mathrm{m}$], $Q = 0$")
ax[1, 0].set_title(r"$\sigma_h$")
ax[0, 1].set_title(r"$h = h(\lambda, \phi)$ [$\mathrm{m}$]")
ax[1, 1].set_title(r"$\sigma_h$")
ax[0, 2].set_title(r"$Q = Q(\lambda, \phi)$ [$\mathrm{W} / \mathrm{m}^2$]")
ax[1, 2].set_title(r"$\sigma_Q$")
ax[0, 3].set_title(r"$F_{tot} = F_{tot}(\lambda, \phi)$ [$\mathrm{W} / \mathrm{m}^2$]")

fig.savefig("img/case2_h_hQ.png", dpi=200)

