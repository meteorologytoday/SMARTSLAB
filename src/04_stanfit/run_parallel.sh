#!/bin/bash

julia ./src/04_stanfit/stanfit_KT_a_lon.jl sh_test 1 &
julia ./src/04_stanfit/stanfit_KT_a_lon.jl sh_test 2 &
julia ./src/04_stanfit/stanfit_KT_a_lon.jl sh_test 3 &
julia ./src/04_stanfit/stanfit_KT_a_lon.jl sh_test 4 &
