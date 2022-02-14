#ifndef INPUT_HANDLER_H
#define INPUT_HANDLER_H

#include <iostream>

#include "INIReader.h"

class InputHandler
{
	public:
        InputHandler(std::string& inputFile, int run_num, long randSeed);
        ~InputHandler(){};

    int read();

    int L;
	double J;

	int noTemperatures;
	double T_min;
	double T_max;
	int noSamples;
    int sweepsBetweenSamples;
    int equilibrationSweeps;

    int runNo;
    long randomSeed;

    private:
        std::string inputFileName;
};

InputHandler::InputHandler(std::string& inputFile, int run_num, long randSeed) : inputFileName(inputFile), runNo(run_num), randomSeed(randSeed)
{
}

int InputHandler::read()
{
    INIReader reader(inputFileName);

    if (reader.ParseError() != 0) {
        return 1;
    }

    L = reader.GetInteger("Model_Parameters", "L", -1);
    J = reader.GetReal("Model_Parameters", "J", -1);

    noTemperatures = reader.GetInteger("Monte_Carlo_Parameters", "no_of_temperatures", -1);
    T_min = reader.GetReal("Monte_Carlo_Parameters", "T_min", -1);
    T_max = reader.GetReal("Monte_Carlo_Parameters", "T_max", -1);
    noSamples = reader.GetInteger("Monte_Carlo_Parameters", "no_of_samples", -1);
    equilibrationSweeps = reader.GetInteger("Monte_Carlo_Parameters", "equilibration_sweeps", -1);
    sweepsBetweenSamples = reader.GetInteger("Monte_Carlo_Parameters", "sweeps_between_samples", -1);

    return 0;
}

#endif