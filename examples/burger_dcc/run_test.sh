#! /bin/bash
#include ../../Makefile.inc
DCC=~/svn/new_dcc/src/dcc
CPP=mpicxx
CFLAGS="-ansi -g -O0 -std=c++0x"

$DCC -a f.c 
cp a1_f.c a1_f_dcc.c
#cpp -I. -P -C b1_f_tmp.c >b1_f.c
sed '/#include/g' a1_f.c > a1_f_tmp.c
mv a1_f_tmp.c a1_f.c
sed '/f(nx/g' a1_f.c > a1_f_tmp.c
mv a1_f_tmp.c a1_f.c

# add extern
sed 's/int myid=0/extern int myid/g' a1_f.c > a1_f_tmp.c
mv a1_f_tmp.c a1_f.c
sed 's/int numprocs=0/extern int numprocs/g' a1_f.c > a1_f_tmp.c
mv a1_f_tmp.c a1_f.c

$CPP $CFLAGS -I../../include -O2 -c dcc_alloc.cpp
$CPP $CFLAGS -I../../include -O2 -c dcc_mpi.cpp
$CPP $CFLAGS -I../../include -O2 -L../../src -o main_1.exe main_1.cpp dcc_mpi.o dcc_alloc.o -lm -lAMPI 
mpirun -np 2 ./main_1.exe
diff tlm.ref adm0.out
diff tlm.ref adm1.out


sed '/#include/g' a1_f_dcc.c > a1_f_dcc_tmp.c
mv a1_f_dcc_tmp.c a1_f_dcc.c
$DCC -t -d 2 a1_f_dcc.c 
mv t2_a1_f_dcc.c t2_a1_f.c
# add extern
sed '/f(nx/g' t2_a1_f.c > t2_a1_f_tmp.c
mv t2_a1_f_tmp.c t2_a1_f.c
sed 's/int myid=0/extern int myid/g' t2_a1_f.c > t2_a1_f_tmp.c
mv t2_a1_f_tmp.c t2_a1_f.c
sed 's/int numprocs=0/extern int numprocs/g' t2_a1_f.c > t2_a1_f_tmp.c
mv t2_a1_f_tmp.c t2_a1_f.c
#$CPP -I../../include -O2 -o main_2.exe main_2.cpp 
$CPP $CFLAGS -I../../include -O2 -L../../src -o main_2.exe main_2.cpp dcc_mpi.o dcc_alloc.o -lm -lAMPI 
mpirun -np 2 ./main_2.exe
diff sotlm.ref soadm0.out
diff sotlm.ref soadm1.out
echo "Think about substracting the partial derivatives!"
