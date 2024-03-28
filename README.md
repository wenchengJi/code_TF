
### C program code of a coarse-grained MD simulation 
1. System requirements
   1.1 All software dependencies and operating systems (including version numbers):
      (a) gcc compiler version 14.0
      (b) command line environment

   1.2 Versions the software has been tested on:
       gcc compiler version 14.0 tested

   1.3 Any required non-standard hardware：
       No

2. Installation guide
   2.1 Instructions:
       In the command line environment:  brew install gcc.

   2.3 Typical install time on a "normal" desktop computer:
       It takes several hours.


3. Demo
   3.1 Instructions to run on data:
       (1) Compile the code  *.c in the command line environment: gcc -lm polymerMD.c -o a.out
       (2) Run the a.out file in the command line environment: ./a.out  ${EB}  ${n}  ${L}  ${k}
         EB is binding energy. 
         n is the number of beads.
         L is the antenna length.
         k is the random seed number.
         Example: ./a.out 11 4 300 1

   3.2 Expected output: 
       Output is in the “Search_n=*.dat” file (An example "Search_n=4_EB=11.00_L=300.0_2.dat").
       The output file has three columns:
       First column: Time;
       Second column: Number of beads bound;
       Third column:  Distance to the DBD target along the antenna direction

   3.3 Expected run time for demo on a "normal" desktop computer
       It takes from a few minutes to a few hours because of its randomness.


 4. Instructions for use
    How to run the software on your data:
    We do not have any imput data . Our data are output from the numerical simulations.

