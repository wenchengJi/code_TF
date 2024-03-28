
### C program code of a coarse-grained MD simulation 

How to install gcc compiler (version 14.0 tested)
In the command line environment:  brew install gcc

How to run
(1) Compile the code  *.c in the command line environment
gcc -lm polymerMD.c -o a.out

(2) Run the a.out file in the command line environment: ./a.out  ${EB}  ${n}  ${L}  ${k}

EB is binding energy. 
n is the number of beads.
L is the antenna length.
k is the random seed number.

Example
./a.out 11 4 300 1

Output in “Search_n=*.dat”  (An example: Search_n=4_EB=11.00_L=300.0_2.dat)

Three columns:
Time;  Number of beads bound;  distance to the DBD target along the antenna direction

