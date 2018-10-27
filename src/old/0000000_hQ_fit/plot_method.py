import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np

import sys

print("The first argument is %s" % (sys.argv[1],))
method = int(sys.argv[1])

filename_hQ = "./data/hQ_fit_method_%d.jl.nc" % method
filename_TF = "./data/dT_star_dt-TOT_F.nc"

print("Going to read %s and %s" % (filename_hQ, filename_TF))

fh_hQ = nc4.Dataset(filename_hQ, "r")
fh_TF = nc4.Dataset(filename_TF, "r")

rlats      = fh_hQ.variables['rlat'][:]
rlons      = fh_TF.variables['rlon'][:]

levs = {}
cmap = {}
if method == 1:

    levs["h"]     = np.arange(-200, 210, 40)
    levs["h_std"] = np.arange(0, 210, 40)

    levs["Q"]     = np.arange(-200, 210, 20)
    levs["Q_std"] = np.arange(0, 210, 20)

    r_cmap  = plt.get_cmap("RdBu_r")
    r_levs  = np.arange(-50, 55, 10)
    r2_levs = np.arange(-200, 210, 20)


elif method == 2:

    levs["h"]     = np.arange(-50, 55, 10) * 2
    levs["h_std"] = np.arange(0, 55, 10) * 2

    levs["Q"]     = np.arange(-50, 55, 10) 
    levs["Q_std"] = np.arange(0, 55, 10) /2

    r_cmap  = plt.get_cmap("RdBu_r")
    r_levs  = np.arange(-50, 55, 10)
    r2_levs = np.arange(-20, 21, 2) / 10


elif method == 3:

    levs["h"]     = np.arange(-50, 55, 10) * 2
    levs["h_std"] = np.arange(0, 55, 5) * 2

    levs["Q"]     = np.arange(-5000, 5500, 500) 
    levs["Q_std"] = np.arange(0, 5000, 250)

    r_cmap  = plt.get_cmap("RdBu_r")
    r_levs  = np.arange(-50, 55, 10)
    r2_levs = np.arange(-1000, 1100, 200)




levs["TOT_F"]     = np.arange(-150, 151, 30) 
levs["TOT_F_std"] = np.arange(0, 151, 15)

levs["dT_star_dt"]     = np.arange(-5, 5.5, .5) 
levs["dT_star_dt_std"] = np.arange(0,  2.1, .2)

cmap["h"]     = plt.get_cmap("RdBu")
cmap["h_std"] = plt.get_cmap("jet")

cmap["Q"]     = plt.get_cmap("RdBu_r")
cmap["Q_std"] = plt.get_cmap("jet")

cmap["TOT_F"]     = plt.get_cmap("RdBu_r")
cmap["TOT_F_std"] = plt.get_cmap("jet")

cmap["dT_star_dt"]     = plt.get_cmap("RdBu_r")
cmap["dT_star_dt_std"] = plt.get_cmap("jet")






for i in range(12):
    month = i+1
    
    print("Plotting month %02d... " % month, end="")


    fig, ax = plt.subplots(2, 5, figsize=(30, 12))
    fig.subplots_adjust(**{
        'left' : 0.05,
        'right': 0.95
    })


    fig.suptitle("Month : %02d (%s)" % (month, ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"][i]), size=30)

    ind = (i, slice(None), slice(None))

    for j, (varname, unit, fh) in enumerate([
        ("dT_star_dt", r"$\mathrm{W} / \mathrm{m}^3$", fh_TF),
        ("h",          r"$\mathrm{m}$",                fh_hQ),
        ("Q",          r"$\mathrm{W} / \mathrm{m}^2$", fh_hQ),
        ("TOT_F",      r"$\mathrm{W} / \mathrm{m}^2$", fh_TF),
    ]):

        ax[0, j].set_title("%s [%s]" % (varname, unit))
        ax[1, j].set_title("%s std [%s]" % (varname, unit))
 
        varname_std = varname + "_std"

        mappable_0 = ax[0, j].contourf(rlons, rlats, fh.variables[varname][ind],     levs[varname],     cmap=cmap[varname]    , extend="both")
        mappable_1 = ax[1, j].contourf(rlons, rlats, fh.variables[varname_std][ind], levs[varname_std], cmap=cmap[varname_std], extend="max")

        fig.colorbar(mappable_0, ax=ax[0, j], ticks=levs[varname])
        fig.colorbar(mappable_1, ax=ax[1, j], ticks=levs[varname_std])



    dh_dt = (fh_hQ.variables["h"][(i+1) % 12, :, :] - fh_hQ.variables["h"][(i-1) % 12, :, :]) / (2 * 365.0 * 86400.0 / 12.0)

    # Residue including dh_dt
    residue = fh_hQ.variables["h"][ind] * fh_TF.variables["dT_star_dt"][ind] \
            + dh_dt                     * fh_TF.variables["T_star"][ind]   \
            - fh_TF.variables["TOT_F"][ind] - fh_hQ.variables["Q"][ind]

    ax[0, 4].set_title(r"Residue [$\mathrm{W} / \mathrm{m}^2$]" + "\ndT_star_dt * h + dh_dt * T_star - Q - TOT_F")
    mappable = ax[0, 4].contourf(rlons, rlats, residue, r_levs, cmap=r_cmap, extend="both")
    fig.colorbar(mappable, ax=ax[0, 4], ticks=r_levs)

    # Residue including dh_dt
    residue = fh_hQ.variables["h"][ind] * fh_TF.variables["dT_star_dt"][ind] \
            - fh_TF.variables["TOT_F"][ind] - fh_hQ.variables["Q"][ind]

    ax[1, 4].set_title(r"Residue [$\mathrm{W} / \mathrm{m}^2$]" + "\ndT_star_dt * h - Q - TOT_F")
    mappable = ax[1, 4].contourf(rlons, rlats, residue, r2_levs, cmap=r_cmap, extend="both")
    fig.colorbar(mappable, ax=ax[1, 4], ticks=r2_levs)

   

    fig.savefig("img/hQ_method_%d_month_%02d.png" % (method, month), dpi=200)
    #fig.savefig("img/hQ_method_%d_month_%02d.eps" % (method, month), dpi=200)
    print("done.")

    plt.close(fig)

