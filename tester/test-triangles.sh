#!/bin/bash

pwd

mkdir -p scratch

n=1000
m=100000


time ./gen/gen-triangles $n $m ./scratch/graph.txt ./scratch/actual-triangles.txt
echo "Generated test data and correct result"

./gen/graph-to-bin ./scratch/header.dat ./scratch/graph.gra < ./scratch/graph.txt

threads=6
time mpirun -n $threads ./a3 --startk 1 --endk 25 --inputpath "./scratch/graph.gra" --headerpath "./scratch/header.dat" --outputpath ./scratch/_out > ./scratch/computed-triangles.txt

diff --brief ./scratch/computed-triangles.txt ./scratch/actual-triangles.txt
comp_value=$?

if [[ $comp_value -eq 1 ]]
then
    echo "Testcase differs"
else
    echo "Passed test"
    rm ./scratch/*
fi
