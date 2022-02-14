#ifndef OUTPUT_HANDLER_H
#define OUTPUT_HANDLER_H

#include <iostream>
#include <vector>

#include <sys/stat.h> 
#include <sys/types.h>

#include "Vec.h"
#include "Replica.h"

typedef Vec<3> Vec3;
typedef std::vector<int> int_list;
typedef std::vector<double> double_list;
typedef std::vector<Vec3> vec_list;
typedef std::vector<int_list> int_int_list;

template <class Model>
class OutputHandler
{
	public:
        OutputHandler(Model& mod, double_list& T_vals, int run_num);
        ~OutputHandler(){};

    void save(std::vector<Replica>& replica, int sampleNo);

    private:
        Model& model;
        double_list& T;
        int runNo;

        //std::string mainDirName;
        std::vector<std::string> subDirPaths;
};

template <class Model>
OutputHandler<Model>::OutputHandler(Model& mod, double_list& T_vals, int run_num) : model(mod), T(T_vals), runNo(run_num)
{
    std::ostringstream os_mainDirName;
    os_mainDirName << "data" << model.simTitle() << "runNo_" << runNo;
    std::string mainDirName = os_mainDirName.str();

    // Create main data directory
    mkdir(mainDirName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    // Create sub directories
    subDirPaths.resize(T.size());

    for (int n=0; n<T.size(); n++)
    {
        std::ostringstream os_subDirName;
        os_subDirName << "temp_id_" << n;
        subDirPaths[n] = mainDirName + "/" + os_subDirName.str();
        mkdir(subDirPaths[n].c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }


    // Save the temperature set 
    std::string filePath = mainDirName + "/" + "temp_values.txt";
    std::ofstream  fout(filePath.c_str());

    for (int n=0; n<T.size(); n++)
        fout << std::setw(20) << std::left << std::setprecision(12) << T[n] << std::endl;

    fout.close();
}

template <class Model>
void OutputHandler<Model>::save(std::vector<Replica>& replica, int sampleNo)
{
    for (int n=0; n<T.size(); n++)
    {
        std::ostringstream os_fileName;
        os_fileName << "temp_id_" << n << "_config_" << sampleNo << ".txt";
        std::string filePath =  subDirPaths[n] + "/" + os_fileName.str();

        std::ofstream fout(filePath.c_str());  

        for (int i=0; i<model.noOfSpins(); i++)
        {
		    fout	<< std::setw(20) << std::left << std::setprecision(12) << std::fixed << replica[n].spin[i][0]
			        << std::setw(20) << std::left << std::setprecision(12) << std::fixed << replica[n].spin[i][1]
                    << std::setw(20) << std::left << std::setprecision(12) << std::fixed << replica[n].spin[i][2]
			        << std::endl;
        }

        fout.close();
    }

}

#endif