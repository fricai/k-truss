#!/bin/bash

n=1000
m=100000
k1=1
k2=20
ranks=6

gen_test=0

if [[ $gen_test -eq 1 ]]
then
    mkdir -p scratch
    rm scratch/*

    time ./gen/gen-testcase $n $m $k1 $k2 ./scratch/graph.txt ./scratch/actual-truss.txt
    echo "Generated test data and correct result"

    ./gen/graph-to-bin ./scratch/header.dat ./scratch/graph.gra < ./scratch/graph.txt
fi

time mpirun -n $ranks ./a3 --startk $k1 --endk $k2 --inputpath "./scratch/graph.gra" --headerpath "./scratch/header.dat" --outputpath "./scratch/_out" > ./scratch/computed-truss.txt

diff --brief ./scratch/computed-truss.txt ./scratch/actual-truss.txt
comp_value=$?

if [[ $comp_value -eq 1 ]]
then
    echo "Testcase differs"
else
    echo "Passed test"
fi
