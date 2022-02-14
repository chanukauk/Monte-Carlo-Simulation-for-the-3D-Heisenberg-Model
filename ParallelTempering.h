#ifndef PARALLEL_TEMPERING_H
#define PARALLEL_TEMPERING_H

#include <iostream>
#include <cmath>
#include <vector>
#include <queue>
#include <numeric>
#include <algorithm>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "Vec.h"
#include "Replica.h"
#include "MCMoves.h"
#include "OutputHandler.h"

typedef Vec<3> Vec3;
typedef std::vector<int> int_list;
typedef std::vector<double> double_list;
typedef std::vector<Vec3> vec_list;
typedef std::vector<int_list> int_int_list;

class Value
{
	public:
        Value();
		Value(int ind, double val) : index(ind), value(val) {}

        int index;
        double value;

        bool operator < (const Value& aVal) const
        {
            return (value < aVal.value);
        }
};



class ParallelTempering
{
	public:
		ParallelTempering();
		~ParallelTempering(){}

		template <class Model, class RndGen>
        void optimizeTemperatureSet(Model& model, int noReplicas, double T_min, double T_max, int noIterations, int runNo, RndGen& random);

        template <class Model, class RndGen>
        void run(Model& model, int noReplicas, double T_min, double T_max, 
                            int noSamples, int equilibrationSweeps, int sweepsBetweenSamples, int runNo, RndGen& random);

	private:

		MCMoves mcmc;
};


ParallelTempering::ParallelTempering() : mcmc()
{
}


template <class Model, class RndGen>
void ParallelTempering::run(Model& model, int noReplicas, double T_min, double T_max, 
                            int noSamples, int equilibrationSweeps, int sweepsBetweenSamples, int runNo, RndGen& random)
{
    double_list T(noReplicas, 1.0);

   
    // Linear temperature schedule
    double dT = (T_max-T_min)/(noReplicas-1);

    for (int i=0; i<noReplicas; i++)
        T[i] = T_min + i*dT;
    

    /*
    // Geometric temperature schedule
    double r = pow( T_max/T_min, 1.0/(noReplicas-1.0) );
    double mult = r;
    double lamda = 1;

    T[0] = T_min;
    T[noReplicas-1] = T_max;

    for (int i=1; i<(noReplicas-1); i++)
    { 
        T[i] = lamda*T_min*mult;

        mult = mult*r;
    }
    */
    

	//Generate replicas
    int noSpins = model.noOfSpins();
	Replica rep(noSpins, 0);

	std::vector<Replica> replica(noReplicas, rep);

    int_list currentReplicaId(noReplicas, 0);

    // Initialize replicas
    Vec3 upSpin;
    upSpin[0] = 0; upSpin[1] = 0; upSpin[2] = 1;

    for (int n=0; n<noReplicas; n++)
    {
        for (int i=0; i<noSpins; i++)
            replica[n].spin[i].randomizeDir(random);
    
        // Calculate the energy of the replica
        double E = model.energy(replica[n].spin);

        replica[n].repId = n;
        replica[n].energy = E;
        
        currentReplicaId[n] = n;
    }

    double_list pairSwapProb(noReplicas-1, 0.0);
    int_list noPairTrialSwaps(noReplicas-1, 0);


    int sampleNo = 0;

    int nextSampleAtIter = equilibrationSweeps;

    int noIterations = equilibrationSweeps + noSamples*sweepsBetweenSamples;


    OutputHandler<Model> output(model, T, runNo);

    std::cout << "Equilibrating system ..." << std::endl;

    for (int iter=0; iter<noIterations; iter++)
    {
        // Perform MC  sweeps
		for (int n=0; n<noReplicas; n++)
            mcmc.compoundSweep(replica[n].spin, model, 1.0/T[n], replica[n].energy, random);
        
        // Perform replica swaps 
        for (int n=0; n<(noReplicas-1); n++)
        {
            double delta = (1.0/T[n+1] - 1.0/T[n])*(replica[n+1].energy - replica[n].energy);

            noPairTrialSwaps[n]++;

            if ( (delta > 0) || ( random() <= exp(delta) )  )
            {
                replica[n].swap(replica[n+1]);
                pairSwapProb[n] += 1.0;
            }
        }          
        
        for (int n=0; n<noReplicas; n++)
            currentReplicaId[ replica[n].repId ] = n;

        if (iter == nextSampleAtIter)
        {
            // write data       
            output.save(replica, sampleNo);
            nextSampleAtIter += sweepsBetweenSamples; 

            sampleNo++;

            std::cout << "Sample " << sampleNo << " of " << noSamples << " taken." << std::endl;
            
            if (sampleNo==noSamples)
            {
                std::cout << "Done!" << std::endl;
                break;
            }
        }

    } // End - for (int iter=0; iter<noIterations; iter++)

}




template <class Model, class RndGen>
void ParallelTempering::optimizeTemperatureSet(Model& model, int noReplicas, double T_min, double T_max, int noIterations, int runNo, RndGen& random)
{
    double minSwapProb = 0.3;

    double_list T(noReplicas, 1.0);

    // Linear temperature schedule
    double dT = (T_max-T_min)/(noReplicas-1);

    for (int i=0; i<noReplicas; i++)
        T[i] = T_min + i*dT;


	//Generate replicas
    int noSpins = model.noOfSpins();
	Replica rep(noSpins, 0);

	std::vector<Replica> replica(noReplicas, rep);

    int_list currentReplicaId(noReplicas, 0);

    // Initialize replicas
    Vec3 upSpin;
    upSpin[0] = 0; upSpin[1] = 0; upSpin[2] = 1;

    for (int n=0; n<noReplicas; n++)
    {
        for (int i=0; i<noSpins; i++)
            replica[n].spin[i].randomizeDir(random);
    
        // Calculate the energy of the replica
        double E = model.energy(replica[n].spin);

        replica[n].repId = n;
        replica[n].energy = E;
        
        currentReplicaId[n] = n;
    }



    for (int iter=0; iter<200; iter++)
    {
        double_list pairSwapProb(noReplicas-1, 0.0);
        int_list noPairTrialSwaps(noReplicas-1, 0);

        // Caculate acceptance ratios for the swaps (swap probabilities)
        for (int ct=0; ct<100; ct++)
        {
            // Perform MC  sweeps
            for (int n=0; n<noReplicas; n++)
		        mcmc.compoundSweep(replica[n].spin, model, 1.0/T[n], replica[n].energy, random);

            // Perform replica swaps 
            for (int n=0; n<(noReplicas-1); n++)
            {
                double delta = (1.0/T[n+1] - 1.0/T[n])*(replica[n+1].energy - replica[n].energy);

                noPairTrialSwaps[n]++;

                if ( (delta > 0) || ( random() <= exp(delta) )  )
                {
                    replica[n].swap(replica[n+1]);
                    pairSwapProb[n] += 1.0;
                }
            }          

        } // End - for (int ct=0; ct<10; ct++)


        for (int n=0; n<noReplicas; n++)
            currentReplicaId[ replica[n].repId ] = n;

        for (int n=0; n<(noReplicas-1); n++)
                pairSwapProb[n] /= noPairTrialSwaps[n];
    

        std::cout << std::endl;

        for (int n=0; n<(noReplicas-1); n++)
            std::cout << std::setw(10) << std::left << std::setprecision(5) << pairSwapProb[n];

        std::cout << std::endl;

        
        int noIndicesToSelect = (int)(0.2*noReplicas);
        
        if (noIndicesToSelect == 0)
            noIndicesToSelect = 1;

 
        // Check if all success probabilities are above the required minimum value
        int condReached = 1;

        for (int n=0; n<(noReplicas-1); n++)
        {
            if (pairSwapProb[n] < minSwapProb)
            {
                condReached = 0;
    
                break;
            }
        }

        if (condReached)
        {
            std::cout << "\nOptimal temperature set reached with " << noReplicas << " temperatures!\n" << std::endl;

            std::ostringstream fileName;
            fileName << "optimized_temp_schedule" << model.simTitle() << "runNo_" << runNo << ".dat";
            std::ofstream  fout(fileName.str().c_str());

            for (int n=0; n<noReplicas; n++)
                fout << std::setw(10) << std::left << std::setprecision(5) << T[n];

            fout.close();

            break;
        }


        std::vector<Value> lowSwapProbs; 

        for (int n=0; n<(noReplicas-1); n++)
        {
            if (pairSwapProb[n] < minSwapProb)
                lowSwapProbs.push_back(Value(n, pairSwapProb[n]));
        }


        std::sort(lowSwapProbs.begin(), lowSwapProbs.end(), std::less<Value>());

        if (lowSwapProbs.size() < noIndicesToSelect)
            noIndicesToSelect = lowSwapProbs.size();

        int_list selectedIndices(noIndicesToSelect);

        for (int n=0; n<noIndicesToSelect; n++)
            selectedIndices[n] = lowSwapProbs[n].index;

        std::sort(selectedIndices.begin(), selectedIndices.end());




        // Add more temperatures into the middle 
        int_list insertIndex;
        double midT;
        double_list insertT;
        int repIndex;

        for (int n=0; n<selectedIndices.size(); n++)
        {
            repIndex = selectedIndices[n];
            insertIndex.push_back(repIndex + 1);

            midT = 0.5*(T[repIndex] + T[repIndex+1]);

            insertT.push_back(midT);
        }

         

        int index;

        for (int i=0; i<insertIndex.size(); i++)
        {
           index =  insertIndex[i] + i;

           T.insert(T.begin() + index, insertT[i]);
           replica.insert(replica.begin() + index, replica[index]);
           
           currentReplicaId.push_back(0);

           noReplicas++;

        }
 


    }

}
// End - Adaptive method for optimizing J set //

#endif

