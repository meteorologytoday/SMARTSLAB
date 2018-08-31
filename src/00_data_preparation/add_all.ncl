begin


;varnames = (/"rsntds", "rlds", "hfss", "hfls", "hfibthermds", "hfrunoffds", "hfsifrazil"/)

varnames = (/                        \
    "total_downward_heat_flux",      \
    "S",                             \
    "B",                             \
    "tos",                           \
    "hfds"                           \
/)


var_comps_list = NewList("lifo")
ListAppend(var_comps_list, (/"rsntds", "rlds", "hfss", "hfls"/)) 
ListAppend(var_comps_list, (/"rsntds"/)) 
ListAppend(var_comps_list, (/"rlds", "hfss", "hfls"/)) 
ListAppend(var_comps_list, (/"tos"/)) 
ListAppend(var_comps_list, (/"hfds"/)) 


year_interval = 5
beg_years = ispan(1861, 2001, year_interval)
;beg_years := ispan(1861, 1891, year_interval)
;
eg_years := ispan(1861, 1866, year_interval)

signature = "_Omon_GFDL-ESM2G_historical_r1i1p1_"


ref_fh = addfile("hfss" + signature + "186101-186512.nc", "r")
dim = dimsizes(ref_fh->hfss)
tmp = new(dim, "float")
copy_VarCoords(ref_fh->hfss, tmp)


time = todouble(ispan(1, year_interval * 12 * dimsizes(beg_years), 1))

dim := (/dimsizes(varnames), dimsizes(time), dim(1), dim(2)/)
output_vars = new(dim, "float")

print("Begin loop")
do i = 0, dimsizes(beg_years) - 1
    beg_year = beg_years(i)
    end_year = beg_year + 4


    do j = 0, dimsizes(varnames) - 1
        var_comps := var_comps_list[j]
        tmp = 0 ; clear all value
        
        do k = 0, dimsizes(var_comps) - 1
        
            var_comp = var_comps(k)
            filename = var_comp + signature + beg_year + "01-" + end_year + "12.nc"
            print("Reading file:" + filename)
            fh = addfile(filename, "r")

            tmp = tmp + fh->$var_comp$

            delete(fh)

        end do

        output_vars(j, year_interval * 12 * i : year_interval * 12 * (i+1) - 1, :, :) = tmp

    end do
    
end do




; output data




fout_name = "SMART" + signature + beg_years(0) + "01-" + beg_years(dimsizes(beg_years)-1) + "12.nc"
system("rm -f " + fout_name)
fout = addfile(fout_name, "c")
fout->lat = ref_fh->lat
fout->lon = ref_fh->lon
fout->vertices = ref_fh->vertices
fout->lon_vertices = ref_fh->lon_vertices
fout->lat_vertices = ref_fh->lat_vertices
fout->time = time

do j = 0, dimsizes(varnames)-1
    varname = varnames(j)
    print("Writing " + varname)
    x = output_vars(j, :, :, :)
    x!0 = "time"
    x&time = time
    fout->$varname$ = x 
end do

;fout->test = output_vars[0] - output_vars[1]

delete(fout)

delete(ref_fh)
end 
