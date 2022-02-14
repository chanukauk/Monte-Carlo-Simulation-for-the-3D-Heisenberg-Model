#ifndef HEISENBERG_MODEL_H
#define HEISENBERG_MODEL_H

#include <iostream>
#include <vector>

#include "Vec.h"
#include "Weight.h"
#include "Neighbor.h"

typedef Vec<3> Vec3;
typedef std::vector<int> int_list;
typedef std::vector<double> double_list;
typedef std::vector<Vec3> vec_list;
typedef std::vector<Weight> weight_list;
typedef std::vector<Neighbor> neighbor_list;
typedef std::vector<int_list> int_int_list;


class HeisenbergModel
{
	public:
		HeisenbergModel(int Lval, double Jval);
		~HeisenbergModel(){}

        int noOfSpins() { return noSpins; }

        std::string& simTitle() { return title; }

        double energy(vec_list& spin);
        double energy(vec_list& spin, int i);

        Vec3 localField(vec_list& spin, int i);

        double magnetization(vec_list& spin);

	private:
		int getIndex(int i, int j, int k);
        void buildModel();
        void buildNList();

		int L;
		int noSpins;

        double J_val;

        weight_list J;
		std::vector<neighbor_list> nList;

        std::string title;
};


HeisenbergModel::HeisenbergModel(int Lval, double Jval) : L(Lval), J_val(Jval), noSpins(L*L*L)
{
    buildModel();

	std::ostringstream titleOstring;
    titleOstring << "_heisenberg_L_" << L << "_J_" << J_val << "_";
	title = titleOstring.str();
}


void HeisenbergModel::buildModel()
{
	int ipos, ineg, jpos, jneg, kpos, kneg;
	int myIndex, nIndex;

	for (int i=0; i<L; i++)
	{
		for (int j=0; j<L; j++)
		{	
			for (int k=0; k<L; k++)
			{
				myIndex = getIndex(i, j, k); 
				
				ipos = (i+1)%L;
				ineg = (i+L-1)%L;
				jpos = (j+1)%L;
				jneg = (j+L-1)%L;
				kpos = (k+1)%L;
				kneg = (k+L-1)%L;
				
				nIndex = getIndex(ipos, j, k);
				
				if (nIndex > myIndex)
				{
					J.push_back( Weight(myIndex, nIndex, J_val) );
				}

				nIndex = getIndex(ineg, j, k);

				if (nIndex > myIndex)
				{
					J.push_back( Weight(myIndex, nIndex, J_val) );
				}

				nIndex = getIndex(i, jpos, k);

				if (nIndex > myIndex)
				{
					J.push_back( Weight(myIndex, nIndex, J_val) );
				}

				nIndex = getIndex(i, jneg, k);

				if (nIndex > myIndex)
				{
					J.push_back( Weight(myIndex, nIndex, J_val) );
				}

				nIndex = getIndex(i, j, kpos);

				if (nIndex > myIndex)
				{
					J.push_back( Weight(myIndex, nIndex, J_val) );
				}

				nIndex = getIndex(i, j, kneg);

				if (nIndex > myIndex)
				{
					J.push_back( Weight(myIndex, nIndex, J_val) );
				}

			}

		}
	}


    // Build the neighbor list
    buildNList();
}



void HeisenbergModel::buildNList()
{
	nList.resize(noSpins);

	int i, j;

	for (int n=0; n<J.size(); n++)
	{
		i = J[n].index1();
		j = J[n].index2();
		
		Neighbor nb;

		nb.setNeighborIndex(j);
		nb.setWeightIndex(n);

		nList[i].push_back(nb);

		nb.setNeighborIndex(i);
		nb.setWeightIndex(n);

		nList[j].push_back(nb);
	}
}


double HeisenbergModel::energy(vec_list& spin)
{
	double E = 0.0;
	
	int j;
	double J_ij;

	for (int i=0; i<nList.size(); i++)
	{
		for (int nIndex=0; nIndex<nList[i].size(); nIndex++)
		{
			j = nList[i][nIndex].neighborIndex();

			J_ij = J[ nList[i][nIndex].weightIndex() ].value();

			E -= 0.5*J_ij*(spin[i]*spin[j]);
					
		}	
	} // End - for (int i=0; i<nList.size(); i++)

	return E;
}


double HeisenbergModel::energy(vec_list& spin, int i)
{
	double E = 0.0;

	int j;
	double J_ij;

	for (int nIndex=0; nIndex<nList[i].size(); nIndex++)
	{
		j = nList[i][nIndex].neighborIndex();

		J_ij = J[ nList[i][nIndex].weightIndex() ].value();

		E -= J_ij*(spin[i]*spin[j]);				
	}	

	return E;
}


Vec3 HeisenbergModel::localField(vec_list& spin, int i)
{
	Vec3 H;
    H[0] = 0; H[1] = 0; H[2] = 0;

	int j;
	double J_ij;

	for (int nIndex=0; nIndex<nList[i].size(); nIndex++)
	{
		j = nList[i][nIndex].neighborIndex();

		J_ij = J[ nList[i][nIndex].weightIndex() ].value();

		H = H + J_ij*spin[j];				
	}	

	return H;
}



double HeisenbergModel::magnetization(vec_list& spin)
{
	Vec3 M;
    M[0] = 0; M[1] = 0; M[2] = 0;

    for (int i=0; i<noSpins; i++)
        M = M + spin[i];

    return M.mag()/noSpins;
}




int HeisenbergModel::getIndex(int i, int j, int k)
{
	return k + j*L + i*L*L;
}


#endif