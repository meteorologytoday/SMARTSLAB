import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc4

import sys

rlat = float(sys.argv[1])
rlon = float(sys.argv[2])

place_name = None
if len(sys.argv) >= 4:
    place_name = sys.argv[3]


filename_hQ = "./data/case3_hQ.jl.nc"
filename_TF = "./data/dT_star_dt-TOT_F.nc"

print("Going to read %s and %s" % (filename_hQ, filename_TF))


fh_hQ = nc4.Dataset(filename_hQ, "r")
fh_TF = nc4.Dataset(filename_TF, "r")

rlats      = fh_hQ.variables['rlat'][:]
rlons      = fh_hQ.variables['rlon'][:]
times      = fh_hQ.variables['time'][:]

lats      = fh_hQ.variables['lat'][:, :]
lons      = fh_hQ.variables['lon'][:, :]

rlat_i = None
for i in range(len(rlats)-1):
    if rlats[i] <= rlat and rlat <= rlats[i+1]:
        rlat_i = i if (rlat - rlats[i]) < (rlats[i+1] - rlat) else i+1
        break

if rlat_i is None:
    raise Exception("Cannot find suitable latitude index!")

rlon_i = None
for i in range(len(rlons)-1):
    if rlons[i] <= rlon and rlon <= rlons[i+1]:
        rlon_i = i if (rlon - rlons[i]) < (rlons[i+1] - rlon) else i+1
        break

if rlon_i is None:
    raise Exception("Cannot find suitable longitude index!")



rlat = rlats[rlat_i]
rlon = rlons[rlon_i]
lat  = lats[rlat_i, rlon_i]
lon  = lons[rlat_i, rlon_i]
ind = (slice(None), rlat_i, rlon_i)

print("The nearest grid point to where you pick is rlat=%.2f, rlon=%.2f at grid point(%d, %d). Real lat=%.2f, lon=%.2f" % (rlat, rlon, rlat_i, rlon_i, lat, lon))

fig, (ax, ax2, ax3) = plt.subplots(3,1,figsize=(12,8))

ax.grid(True)
ax2.grid(True)
ax3.grid(True)

ax.set_ylabel(r"[$\mathrm{W} \, \mathrm{m}^{-2}$]")
ax2.set_ylabel(r"[$\mathrm{W} \, \mathrm{m}^{-2}$]")

title = "rlat = %.1f, rlon = %.1f" % (lat, lon)
if place_name is not None:
    title = "%s (%s)" % (title, place_name)

ax.set_title(title)
ax.plot(times, fh_hQ.variables["Q"][ind], color="r", label="Q")
ax.plot(times, fh_hQ.variables["dh_dt"][ind] * fh_TF.variables["T_star"][ind], dashes=(10,5,2,5), color="k", label=r"$T^*\partial_t h$")

ax2.plot(times, fh_TF.variables["TOT_F"][ind], color="b", label=r"$F_{tot}$")
ax2.plot(times, fh_hQ.variables["h"][ind] *  fh_TF.variables["dT_star_dt"][ind], dashes=(10,5), color="k", label=r"$h\partial_t T^* $")

h = fh_hQ.variables["h"][ind]
ax3.plot(times, h, color="k", label=r"$h$")
ax3.fill_between(times, h, y2=0, where=h < 0, hatch="//", facecolor="none", interpolate=True)
ax.legend(loc="lower left")
ax2.legend(loc="lower left")
ax3.legend(loc="lower left")

fig.savefig("./img/SinglePoint_[%.2f][%.2f].png" % (lat, lon) , dpi=200)

