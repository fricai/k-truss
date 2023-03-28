test_no=$1

echo "Test $test_no"
for i in $(seq 1 6); do
    echo "Running $i"
    time mpirun -n $i ./a3 --startk=1 --endk=20 --inputpath=data/test$test_no/test-input-$test_no.gra --headerpath=data/test$test_no/test-header-$test_no.dat --verbose=0 --outputpath=op --taskid=$test_no --p=0;
done
