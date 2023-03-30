#!/bin/bash

pwd

mkdir -p scratch

n=1000
m=100000
k1=1
k2=3

gen_test=1

if [[ $gen_test -eq 1 ]]
then
    time ./gen/gen-testcase $n $m $k1 $k2 ./scratch/graph.txt ./scratch/actual-truss.txt
    echo "Generated test data and correct result"

    ./gen/graph-to-bin ./scratch/header.dat ./scratch/graph.gra < ./scratch/graph.txt
fi

threads=1
time mpirun -n $threads ./a3 --startk $k1 --endk $k2 --inputpath "./scratch/graph.gra" --headerpath "./scratch/header.dat" --outputpath "./scratch/_out" > ./scratch/computed-truss.txt

diff --brief ./scratch/computed-truss.txt ./scratch/actual-truss.txt
comp_value=$?

if [[ $comp_value -eq 1 ]]
then
    echo "Testcase differs"
else
    echo "Passed test"
    rm ./scratch/*
fi