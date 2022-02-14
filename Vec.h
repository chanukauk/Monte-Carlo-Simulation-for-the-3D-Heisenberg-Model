/*
Represents a vector with "DIM" dimensions.
No of dimensions (DIM) should be passed as a template parameter.
*/

#ifndef Vec_H
#define Vec_H

# define PI  3.14159265358979323846

#include <iostream>
#include <vector>
#include <cmath>

template <int DIM>
class Vec
{
	public:
		Vec();
		Vec(const Vec<DIM>& aVec);
		Vec(double val);
		~Vec();
	
		const double& operator[](unsigned int i) const;
		double& operator[](unsigned int i);
		
		Vec<DIM> operator+(const Vec<DIM>& aVec) const;
		Vec<DIM> operator+(double scalar) const;
		
		Vec<DIM> operator-() const;
		Vec<DIM> operator-(const Vec<DIM>& aVec) const;
		Vec<DIM> operator-(double scalar) const;
		
		double operator*(const Vec<DIM>& aVec) const;
		Vec<DIM> operator*(double scalar) const;
		
		Vec<DIM> operator/(double scalar) const;
		
		Vec<DIM>& operator=(const Vec<DIM>& aVec);
		Vec<DIM>& operator=(double scalar);

		double mag() const;
		Vec<DIM> unitVec() const;
		void normalize();
		
		Vec<DIM> cross(const Vec<DIM> &aVec) const;
				
		void rotate(const Vec<DIM>& u, double phi);
		
		template < class RNGenerator >
		void randomize(double mag, RNGenerator &random);
		
		template < class RNGenerator >
		void randomizeDir(RNGenerator &random);

		template < class RNGenerator >
		void randomizeDir(double thetaMax, RNGenerator &random);

		template <int D>
		friend Vec<D> operator+(double scalar, const Vec<D>& aVec);
		
		template <int D>
		friend Vec<D> operator-(double scalar, const Vec<D>& aVec);		
		
		template <int D>
		friend Vec<D> operator*(double scalar, const Vec<D>& aVec);
		
	
	protected:
		double component[DIM];
};



// Constructors

template <int DIM>
Vec<DIM>::Vec()
{
	*this = 0;
}

template <int DIM>
Vec<DIM>::Vec(const Vec<DIM>& aVec)
{
	*this = aVec;
}

template <int DIM>
Vec<DIM>::Vec(double val)
{
	*this = val;
}

// Destructor
template <int DIM>
Vec<DIM>::~Vec()
{
}


// Vector magnitude
template <int DIM>
double Vec<DIM>::mag() const
{
	double result = (*this)*(*this);
	
	return sqrt(result);
}

// Return the unit vector in the vector's direction
template <int DIM>
Vec<DIM> Vec<DIM>::unitVec() const
{
	return (*this)/(*this).mag();
}


// Normalize the vector
template <int DIM>
void Vec<DIM>::normalize()
{
	(*this) = this->unitVec();
}


// Overload operator []

template <int DIM>
const double& Vec<DIM>::operator[](unsigned int i) const
{
	return component[i];
}

template <int DIM>
double& Vec<DIM>::operator[](unsigned int i)
{
	return component[i];
}



// Overload operator +

template <int DIM>
Vec<DIM> Vec<DIM>::operator+(const Vec<DIM>& aVec) const
{
	Vec<DIM> temp;
	
	for (int i=0; i<DIM; i++)
		temp[i] = (*this)[i] + aVec[i];

	return temp;		
}

template <int DIM>
Vec<DIM> Vec<DIM>::operator+(double scalar) const
{
	Vec<DIM> temp;

	for (int i=0; i<DIM; i++)
		temp[i] = (*this)[i] + scalar;

	return temp;
}

template <int DIM>
Vec<DIM> operator+(double scalar, const Vec<DIM>& aVec)
{
	return (aVec + scalar);
}



// Overload operator -

template <int DIM>
Vec<DIM> Vec<DIM>::operator-() const
{
	return -1.0*(*this);
}

template <int DIM>
Vec<DIM> Vec<DIM>::operator-(const Vec<DIM>& aVec) const
{
	return -aVec + (*this);
}

template <int DIM>
Vec<DIM> Vec<DIM>::operator-(double scalar) const
{
	return -scalar + (*this);
}

template <int DIM>
Vec<DIM> operator-(double scalar, const Vec<DIM>& aVec)
{
	return -aVec + scalar;
}



// Overload operator *

template <int DIM>
double Vec<DIM>::operator*(const Vec<DIM>& aVec) const
{
	double result = 0;
	
	for (int i=0; i<DIM; i++)
		result += component[i]*aVec.component[i];
		
	return result;
}

template <int DIM>
Vec<DIM> Vec<DIM>::operator*(double scalar) const
{
	Vec<DIM> result;
	
	for (int i=0; i<DIM; i++)
		result[i] = component[i]*scalar;
		
	return result;
}


template <int DIM>
Vec<DIM> operator*(double scalar, const Vec<DIM>& aVec)
{
	return aVec*scalar;
}


// Overload operator / 
template <int DIM>
Vec<DIM> Vec<DIM>::operator/(double scalar) const
{
	return (1.0/scalar)*(*this);		
}



// Overload operator = 

template <int DIM>
Vec<DIM>& Vec<DIM>::operator=(const Vec<DIM>& aVec)
{
	for (int i=0; i<DIM; i++)
		component[i] = aVec.component[i];

	return *this;
}

template <int DIM>
Vec<DIM>& Vec<DIM>::operator=(double scalar)
{
	for (int i=0; i<DIM; i++)
		component[i] = scalar;
	
	return *this;
}


/* Scales the magnitude of the vector by a factor of "mag", and orient the vector
   in a random direction
   CAUTION: Only works for DIM = 3
*/
template <int DIM>
template < class RNGenerator >
void Vec<DIM>::randomize(double mag, RNGenerator &random)
{ 
	(*this).randomizeDir(random);
	
	(*this) = mag*(*this);
}


template <int DIM>
template < class RNGenerator >
void Vec<DIM>::randomizeDir(RNGenerator &random)
{ 
	double v1, v2, s;
	
	do
	{
		v1 = 2*random()-1;
		v2 = 2*random()-1;
		
		s = v1*v1 + v2*v2;
	} while (s>1);
	
	component[0] = 2*v1*sqrt(1-s);
	component[1] = 2*v2*sqrt(1-s);
	component[2] = 1-2*s; 
}

template <int DIM>
template < class RNGenerator >
void Vec<DIM>::randomizeDir(double thetaMax, RNGenerator &random)
{
	Vec<DIM> u, w;
	Vec<DIM>& s = *this;

	u[0] = 1.0; 
	u[1] = 0.0; 
	u[2] = 0.0;

	double dotProd = u*s;
	
	if ( ( fabs(dotProd-1.0) < 1.0e-20 ) ||  ( fabs(dotProd+1.0) < 1.0e-20 ) )
	{
		w[0] = 0.0;
		w[1] = 1.0;	
		w[2] = 0.0;
	}
	else
		w = u - (u*s)/(s*s)*s;

	w.normalize();

	double phi = 2.0*PI*random();		

	w = cos(phi)*w + s.cross(w)*sin(phi);

	double cosThetaMax = cos(thetaMax);

	//double cosTheta = (1.0 - cosThetaMax)*random() + cosThetaMax;	
	double cosTheta = (cosThetaMax-1.0)*random() + 1.0;
	double sinTheta = sqrt(1-cosTheta*cosTheta);

	//std::cout << "old spin - s[o] = " << (*this)[0] << ", s[1] = " << (*this)[1] << ", s[2] = " << (*this)[2] << std::endl;
	Vec<DIM> old = s;
		
	s = cosTheta*s + sinTheta*w;	

	//std::cout << "new spin - s[o] = " << (*this)[0] << ", s[1] = " << (*this)[1] << ", s[2] = " << (*this)[2] << std::endl;
	//std::cout << "angle =  " << acos(old*s)/PI*180 << std::endl;
}


/* Vector cross product 
   CAUTION: Only works for DIM = 3
*/
template <int DIM>
Vec<DIM> Vec<DIM>::cross(const Vec<DIM> &aVec) const
{
	Vec<DIM> result;

	result[0] = (*this)[1]*aVec[2] - (*this)[2]*aVec[1];
	result[1] = (*this)[2]*aVec[0] - (*this)[0]*aVec[2];
	result[2] = (*this)[0]*aVec[1] - (*this)[1]*aVec[0];	
	
	return result;
}


/* Rotates the vector by an angle of "theta" about the axis represented by 
   the unit vector "u" 
   CAUTION: Only works for DIM = 3
*/
template <int DIM>
void Vec<DIM>::rotate(const Vec<DIM>& u, double phi)
{
    	Vec<DIM>& V = *this;
    	Vec<DIM> aTerm = u*(u*V);
    		
    	V = (V - aTerm)*cos(phi) + u.cross(V)*sin(phi) + aTerm;			
}


#endif
