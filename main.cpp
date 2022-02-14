#include <iostream>
#include <cstdlib>
#include <vector>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "InputHandler.h"
#include "Vec.h"

#include "HeisenbergModel.h"
#include "R1279Wrap.h"

#include "ParallelTempering.h"

#include <sys/stat.h> 
#include <sys/types.h>

typedef R1279Wrap RndGen;

typedef Vec<3> Vec3;
typedef std::vector<int> int_list;
typedef std::vector<double> double_list;
typedef std::vector<Vec3> vec_list;

int main(int argc, char* argv[])
{
	if(argc<4)
    {
               std::cerr << "Usage: ./main.exe <input_file_name> <run_number> <random_seed>" << std::endl;
               return EXIT_FAILURE;
    }

	std::string inputFileName = argv[1];
	int runNo = atoi(argv[2]);
	long randSeed = atol(argv[3]);

	InputHandler input(inputFileName, runNo, randSeed);

	if (input.read() != 0)
	{
		std::cout << "\nCannot read input file \""<< inputFileName << "\"" << std::endl;
		return EXIT_FAILURE;
	}

	// Define and initialize the model
	HeisenbergModel model(input.L, input.J);

	// Initialize random number generators
	RndGen random(randSeed);

	ParallelTempering pt;

    pt.run(model, input.noTemperatures, input.T_min, input.T_max, input.noSamples, input.equilibrationSweeps, input.sweepsBetweenSamples, input.runNo, random);

	return 0;
}

