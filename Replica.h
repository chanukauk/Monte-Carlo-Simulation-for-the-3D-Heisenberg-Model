#ifndef REPLICA_H
#define REPLICA_H

#include "Vec.h"

typedef Vec<3> Vec3;
typedef std::vector<int> int_list;
typedef std::vector<Vec3> vec_list;

class Replica
{
	public:
		Replica(int noSpins, int myRepId);
		Replica(const Replica& cpRep);

		~Replica() {}

		void replace(Replica& cpRep);

        void swap(Replica& cpRep);

		vec_list spin;

		int repId;

        double energy;
};

Replica::Replica(int noSpins, int myRepId) : repId(myRepId), energy(0.0)
{
	spin.resize(noSpins);	
}

Replica::Replica(const Replica& cpRep)
{
	repId = cpRep.repId;
	energy = cpRep.energy;
	
	spin.resize(cpRep.spin.size());
	
    for (int i=0; i<spin.size(); i++)
        spin[i] = cpRep.spin[i];
}

void Replica::replace(Replica& cpRep)
{
	repId = cpRep.repId;
	energy = cpRep.energy;
	
    for (int i=0; i<spin.size(); i++)
        spin[i] = cpRep.spin[i];
}

void Replica::swap(Replica& cpRep)
{
    std::swap(repId, cpRep.repId);
    std::swap(energy, cpRep.energy);

    spin.swap(cpRep.spin);
}


#endif
