test_no=$1
threads=$2
time mpirun -n $threads ./a3 --startk=1 --endk=20 --inputpath=data/test$test_no/test-input-$test_no.gra --headerpath=data/test$test_no/test-header-$test_no.dat --verbose=0 --outputpath=op --taskid=$test_no --p=0
