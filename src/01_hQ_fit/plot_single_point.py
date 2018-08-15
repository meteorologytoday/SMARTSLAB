import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc4

import sys

print("The first argument is %s" % (sys.argv[1],))
method = int(sys.argv[1])

print("The second argument is %s" % (sys.argv[2],))
print("The third argument is %s" % (sys.argv[3],))
lat = float(sys.argv[2])
lon = float(sys.argv[3])

place_name = None
if len(sys.argv) >= 5:
    place_name = sys.argv[4]
    print("The fourth argument is %s" % (sys.argv[4],))


filename_hQ = "./data/hQ_fit_method_%d.jl.nc" % method
filename_TF = "./data/dT_star_dt-TOT_F.nc"

print("Going to read %s and %s" % (filename_hQ, filename_TF))


fh_hQ = nc4.Dataset(filename_hQ, "r")
fh_TF = nc4.Dataset(filename_TF, "r")

rlats      = fh_hQ.variables['rlat'][:]
rlons      = fh_TF.variables['rlon'][:]
times      = fh_TF.variables['time'][:]

lat_i = None
for i in range(len(rlats)-1):
    if rlats[i] <= lat and lat <= rlats[i+1]:
        lat_i = i if (lat - rlats[i]) < (rlats[i+1] - lat) else i+1
        break

if lat_i is None:
    raise Exception("Cannot find suitable latitude index!")

lon_i = None
for i in range(len(rlons)-1):
    if rlons[i] <= lon and lon <= rlons[i+1]:
        lon_i = i if (lon - rlons[i]) < (rlons[i+1] - lon) else i+1
        break

if lon_i is None:
    raise Exception("Cannot find suitable longitude index!")



lat = rlats[lat_i]
lon = rlons[lon_i]
ind = (slice(None), lat_i, lon_i)

print("The nearest grid point to where you pick is lat=%.1f, lon=%.1f" % (lat, lon))

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

plt.show()



