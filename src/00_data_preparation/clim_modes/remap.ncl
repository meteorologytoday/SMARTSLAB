begin

undef("remap_5deg")
; Dimension assumed: time, lat, lon
function remap_5deg(v)
local old_lat, old_lon, old_dim, new_lat, new_lon
begin

old_lat = v&$v!1$
old_lon = v&$v!2$

old_dim = dimsizes(v)
old_lat := conform_dims(old_dim(1:2), old_lat, 0)
old_lon := conform_dims(old_dim(1:2), old_lon, 1)


new_lat = fspan(-90.0, 90.0, 36+1)
new_lon = fspan(0, 360, 72+1)
new_lon := new_lon(0:dimsizes(new_lon)-2) + 2.5

new_lat!0 = "lat"
new_lat&lat = new_lat
new_lat@unit = "degree_north"
new_lat@long_name = "latitude coordinate"

new_lon!0 = "lon"
new_lon&lon = new_lon
new_lon@unit = "degree_east"
new_lon@long_name = "longitude coordinate"

new_v = rcm2rgrid(old_lat, old_lon, v, new_lat, new_lon, 1)
new_v!0 = v!0
new_v!1 = v!1
new_v!2 = v!2

new_v&$v!0$ = v&$v!0$
new_v&$v!1$ = new_lat
new_v&$v!2$ = new_lon

return new_v
end


end
