
p=$(dirname "$0")

python3 $p/plot_single_point.py -73.2 -164.5 "SouthernOcean_Seaice_1 where sea ice exists" 
python3 $p/plot_single_point.py -72.1  -39.5 "SouthernOcean_Seaice_2 where sea ice exists"
python3 $p/plot_single_point.py  30.5 -154.5 "NorthPacific_NorthOfMaunaLoa"
python3 $p/plot_single_point.py -29.5 -154.5 "SouthPacific_MirrorOfNorthOfMaunaLoa"
python3 $p/plot_single_point.py   0.0  -89.5 "EasternPacific_Equator"
python3 $p/plot_single_point.py  75.2  -59.5 "DavisStrait"
