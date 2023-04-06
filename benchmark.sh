module load suite/intel/parallelStudio/2020

time mpirun ./a3 \
--taskid 1 \
--verbose 0 \
--inputpath data/A2/test10/test-input-10.gra \
--headerpath data/A2/test10/test-header-10.dat \
--outputpath data/A2/test10/a3_output10_task1.txt \
--p 10 \
--startk 1 \
--endk 5 \

time mpirun ./a3 \
--taskid 2 \
--verbose 0 \
--inputpath data/A2/test10/test-input-10.gra \
--headerpath data/A2/test10/test-header-10.dat \
--outputpath data/A2/test10/a3_output10_task2.txt \
--p 10 \
--startk 1 \
--endk 5 \
