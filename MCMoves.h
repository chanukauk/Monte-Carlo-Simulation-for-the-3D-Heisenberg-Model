#ifndef MC_MOVES_H
#define MC_MOVES_H

#include <iostream>
#include <cmath>
#include <climits>
#include <vector>
#include <queue>

#include "Vec.h"

typedef Vec<3> Vec3;
typedef std::vector<int> int_list;
typedef std::vector<double> double_list;
typedef std::vector<Vec3> vec_list;
typedef std::vector<int_list> int_int_list;


class MCMoves
{
	public:
		MCMoves();
		~MCMoves(){}

		template <class Model, class RndGen>
		double metropolisSweep(vec_list& spin, Model& model, double beta, double& energy, RndGen& random);

        template <class Model, class RndGen>
        int heatbathStep(vec_list& spin, Model& model, double beta, double& energyDiff, RndGen& random, int i);

        template <class Model, class RndGen>
        void heatbathSweep(vec_list& spin, Model& model, double beta, double& energy, RndGen& random);

        template <class Model>
        void microcanonicalSweep(vec_list& spin, Model& model);

        template <class Model, class RndGen>
        double compoundSweep(vec_list& spin, Model& model, double beta, double& energy, RndGen& random);

		template <class Model, class RndGen>
		void metropolisRun(vec_list& spin, Model& model, double beta, double& energy, int noSweeps, RndGen& random);

	private:
		template <class Model, class RndGen>
		int metropolisStep(vec_list& spin, Model& model, double beta, double& energyDiff, RndGen& random);
};


MCMoves::MCMoves()
{
}

template <class Model, class RndGen>
int MCMoves::metropolisStep(vec_list& spin, Model& model, double beta, double& energyDiff, RndGen& random)
{
	int chosenIndex = (int)(random()*spin.size());

    Vec3 oldSpin = spin[chosenIndex];
	
	double oldEnergyTerm = model.energy(spin, chosenIndex);

    // Randomize the chosen spin
    spin[chosenIndex].randomizeDir(random);

	double newEnergyTerm = model.energy(spin, chosenIndex);

	energyDiff = newEnergyTerm - oldEnergyTerm;

        if ( (energyDiff < 0) || ( random() <= exp(-energyDiff*beta) )  )
        {   
                // Keep the new spin configuration
                return 1;
        }   
        else
        {   
            // Restore the old spin configuration
		    spin[chosenIndex] = oldSpin;
    
            return 0;
        }   
}


template <class Model, class RndGen>
double MCMoves::metropolisSweep(vec_list& spin, Model& model, double beta, double& energy, RndGen& random)
{
	int moveAccepted = 0;
    double energyDiff;
    int acceptedMoves = 0;

	for (int i=0; i<spin.size(); i++)
	{
		moveAccepted = metropolisStep(spin, model, beta, energyDiff, random);

        if (moveAccepted)
        {
            acceptedMoves++;
            energy += energyDiff;
        }
	}

    return (double)acceptedMoves/spin.size();
}



template <class Model, class RndGen>
int MCMoves::heatbathStep(vec_list& spin, Model& model, double beta, double& energyDiff, RndGen& random, int i)
{
    double oldEnergyTerm = model.energy(spin, i);

    Vec3 H = model.localField(spin, i);
    double magH = H.mag();

    double r1 = random();
    double r2 = random();

    double phi = 2.0*PI*r1;
    double cos_phi = cos(phi);
    double sin_phi = sin(phi);

    double cos_theta = log(1.0 + r2*(exp(2.0*beta*magH) - 1.0))/(beta*magH) - 1.0;
    double sin_theta = sqrt(1.0-cos_theta*cos_theta);

    double cos_theta_H = H[2]/magH;
    double sin_theta_H = sqrt(1.0-cos_theta_H*cos_theta_H);
    
    double cos_phi_H = H[0]/(magH*sin_theta_H);
    double sin_phi_H = H[1]/(magH*sin_theta_H);

    double a = cos_theta*sin_theta_H + sin_theta*cos_theta_H*cos_phi;
    double b = sin_theta*sin_phi;

    spin[i][0] = a*cos_phi_H - b*sin_phi_H;
    spin[i][1] = a*sin_phi_H + b*cos_phi_H;
    spin[i][2] = cos_theta*cos_theta_H - sin_theta*sin_theta_H*cos_phi;

    double newEnergyTerm = model.energy(spin, i);

    energyDiff = newEnergyTerm - oldEnergyTerm;

    return 1;
}


template <class Model, class RndGen>
void MCMoves::heatbathSweep(vec_list& spin, Model& model, double beta, double& energy, RndGen& random)
{
    double energyDiff;

	for (int i=0; i<spin.size(); i++)
	{
		heatbathStep(spin, model, beta, energyDiff, random, i);

        energy += energyDiff;
	}
}



template <class Model>
void MCMoves::microcanonicalSweep(vec_list& spin, Model& model)
{
    Vec3 H;
    double magH;

    for (int i=0; i<spin.size(); i++)
    {
        H = model.localField(spin, i);

        magH = H.mag();

        spin[i] = -spin[i] + 2.0*(spin[i]*H)/(magH*magH)*H;
    }
}


template <class Model, class RndGen>
double MCMoves::compoundSweep(vec_list& spin, Model& model, double beta, double& energy, RndGen& random)
{
    for (int i=0; i<10; i++)
       microcanonicalSweep(spin, model);

    heatbathSweep(spin, model, beta, energy, random);

    return 0;
}


template <class Model, class RndGen>
void MCMoves::metropolisRun(vec_list& spin, Model& model, double beta, double& energy, int noSweeps, RndGen& random)
{
	for (int n=0; n<noSweeps; n++)
		metropolisSweep(spin, model, beta, energy, random);
}



#endif
