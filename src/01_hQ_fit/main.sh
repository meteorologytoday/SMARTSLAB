#!/bin/bash
set -x

path=$(dirname $0)

julia $path/hQ_fit_method_1.jl
julia $path/hQ_fit_method_2.jl

python3 $path/plot_method.py 1
python3 $path/plot_method.py 2
