# Monte-Carlo-Simulation-for-the-3D-Heisenberg-Model
You will need a C++ compiler to compile it. If you have GNU C++ compiler installed, then you can simply run the makefile by typing "make". It will build an executable named "main.exe".

./main.exe  <input_file_name>  <run_number>  <random_seed>
The input file named "input.txt" specifies the input parameters for the simulation. "run_number" is just a dummy integer to distinguish between multiple runs with the same parameters. "random_seed" is an arbitrary multi digit integer to seed the random number generator.  
I will explain everything in more detail tomorrow during the meeting, but for now, you may run the code by typing the following line:

./main.exe input.txt 1 384318

When you run the code, all output configurations will be saved in a directory named "data_heisenberg____..." within the current working directory. 