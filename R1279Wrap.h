#ifndef R1279Wrap_H
#define R1279Wrap_H

#include <cmath>

#include "R1279.h"

class R1279Wrap
{
	public:
		R1279Wrap(long randomSeed);
		~R1279Wrap(){}
	
		double operator()();
		double gaussRandom();
		double gaussRandom(double mu, double sigma);
	
		long seed() { return randSeed; }

	private:
		long randSeed;
		R1279 randObj;
};

R1279Wrap::R1279Wrap(long randomSeed) : randSeed(randomSeed), randObj(randomSeed)				
{
}

double R1279Wrap::operator()() 
{
	return randObj();	
}

double R1279Wrap::gaussRandom() 
{
	double x1, x2, w, y1, y2;
 
         do {
                 x1 = 2.0 * randObj() - 1.0;
                 x2 = 2.0 * randObj() - 1.0;
                 w = x1 * x1 + x2 * x2;
         } while ( w >= 1.0 );

         w = sqrt( (-2.0 * log( w ) ) / w );

         y1 = x1 * w;
         //y2 = x2 * w;

	return y1;
}

double R1279Wrap::gaussRandom(double mu, double sigma)
{
	return gaussRandom()*sigma + mu;
}

#endif
