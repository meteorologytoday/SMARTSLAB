#!/bin/bash


for i in {1..20}
do
    echo "Welcome $i times abusing CPU... hehehe"
    julia test_cores.jl &
done

echo "Done."
