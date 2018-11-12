/*******************************************************************************
********************************************************************************

	A minimal Point class, which represents a point in a
	Euclidean space of arbitrary dimension.

													John L. Barber, May 2013
													jlbarber@lanl.gov
													505-664-0605

	11/15/2014: Added "Normalize" routine.

	3/1/2016 Added Norm2 routines, and changed the + and - operators so they
	return vectors with a size equal to the LARGEST size of the two input vectors.
	Also, changed the Distance function to a simple "Norm of difference" calculation.

	3/4/2016 Added FillWith.

	2017 Feb 18: Did a little optimization throughout. Largest affect was in "Distance" function,
		which got a speedup of nearly a factor of 100.

	2018 Feb 20: Added formatted output for const Point objects.

********************************************************************************
*******************************************************************************/


#include "Point3D.h"


// Point constructor. Fills r with 3 copies of x.
Point::Point(double x /* = 0 */){

	r[0] = r[1] = r[2] = x;

	return;
}


// Alternate Point constructor. Fills r array with (x, y, z).
Point::Point(double x, double y, double z){

	r[0] = x;
	r[1] = y;
	r[2] = z;

	return;
}


// Alternate Point constructor. Initializes using argument point.
Point::Point(const Point& P){

	r[0] = P.r[0];
	r[1] = P.r[1];
	r[2] = P.r[2];

	return;
}


// Overloaded += operator. Returns reference.
Point& Point::operator+=(const Point& RHS){

	this->r[0] += RHS.r[0];
	this->r[1] += RHS.r[1];
	this->r[2] += RHS.r[2];

	return *this;
}


// Overloaded -= operator. Returns reference.
Point& Point::operator-=(const Point& RHS){

	this->r[0] -= RHS.r[0];
	this->r[1] -= RHS.r[1];
	this->r[2] -= RHS.r[2];

	return *this;
}


// Overloaded *= operator. Returns reference.
Point& Point::operator*=(double x){

	this->r[0] *= x;
	this->r[1] *= x;
	this->r[2] *= x;

	return *this;
}


// Overloaded /= operator. Returns reference.
Point& Point::operator/=(double x){

	this->r[0] /= x;
	this->r[1] /= x;
	this->r[2] /= x;

	return *this;
}


// Overloaded ostream << operator for use with Point objects
ostream& operator<<(ostream &Out, Point &P){

	Out << "(" << P.r[0] << ", " << P.r[1] << ", " << P.r[2] << ")";

	return Out;
}

// Overloaded ostream << operator for use with const Point objects
ostream& operator<<(ostream &Out, const Point &P){

	Out << "(" << P.r[0] << ", " << P.r[1] << ", " << P.r[2] << ")";

	return Out;
}
