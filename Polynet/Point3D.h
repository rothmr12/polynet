/*******************************************************************************
********************************************************************************

	Header file for a minimal Point class, which represents a point in a
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

#ifndef POINT_HEADER
#define POINT_HEADER

#include <iostream>
#include <vector>
#include "math.h"

using namespace std;

class Point{
public:

	double r[3];

	Point(double x = 0);														// Constructor. Fills r array with 3 copies of x
	Point(double x, double y, double z);										// Alternate constructor. Fills r array with (x, y, z).
	Point(const Point& P);														// Alternate constructor. Initializes using argument point.
	Point& operator += (const Point& RHS);										// Overloaded += operator
	Point& operator -= (const Point& RHS);										// Overloaded -= operator
	Point& operator *= (double x);												// Overloaded *= operator
	Point& operator /= (double x);												// Overloaded /= operator
	inline double Norm() const {return sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);}// Euclidean norm of this Point
	inline double Norm2() const {return r[0]*r[0] + r[1]*r[1] + r[2]*r[2];}		// Euclidean norm squared of this Point (to avoid the square root when we can)
	inline void   FillWith(double x = 0){r[0] = r[1] = r[2] = x;}				// Fills the Point's coordinates with 3 copies of x.
	inline void   Normalize(){*this /= this->Norm();}							// Normalizes Point, i.e. sets its magnitude to 1. (Warning: Does NOT check if norm of point is 0)

	friend inline double Norm(const Point& P){return sqrt(P.r[0]*P.r[0] + P.r[1]*P.r[1] + P.r[2]*P.r[2]);}	// Euclidean norm
	friend inline double Norm2(const Point& P){return P.r[0]*P.r[0] + P.r[1]*P.r[1] + P.r[2]*P.r[2];}		// Euclidean norm squared (to avoid the square root when we can)
	friend inline double Dot(const Point& P1, const Point& P2){												// Dot product
		return P1.r[0]*P2.r[0] + P1.r[1]*P2.r[1] + P1.r[2]*P2.r[2];
	}
	friend inline const double Distance(const Point& P1, const Point& P2){ 										// Euclidean distance
		return sqrt((P1.r[0]-P2.r[0])*(P1.r[0]-P2.r[0])+(P1.r[1]-P2.r[1])*(P1.r[1]-P2.r[1])+(P1.r[2]-P2.r[2])*(P1.r[2]-P2.r[2]));
	}

	// Cross product
	friend inline Point Cross3D(const Point& P1, const Point& P2){
		return Point(P1.r[1]*P2.r[2]-P1.r[2]*P2.r[1], P1.r[2]*P2.r[0]-P1.r[0]*P2.r[2], P1.r[0]*P2.r[1]-P1.r[1]*P2.r[0]);
	}

	// 2D Cross product using x and y components to return a scalar 
	friend inline double Cross2D(const Point& P1, const Point& P2){return P1.r[0]*P2.r[1] - P1.r[1]*P2.r[2];}

	// Returns normalized version of the point
	friend inline Point Normalize(const Point& P){return P / P.Norm();}

	friend inline bool  operator==(const Point& P1, const Point& P2){return (P1.r[0]==P2.r[0] && P1.r[1]==P2.r[1] && P1.r[2]==P2.r[2]);}	// Comparison operator
	friend inline bool  operator!=(const Point& P1, const Point& P2){return (P1.r[0]!=P2.r[0] || P1.r[1]!=P2.r[1] || P1.r[2]!=P2.r[2]);}	// "Anti"-comparison operator
	friend inline Point operator+(const Point& P1, const Point& P2){return Point(P1.r[0] + P2.r[0], P1.r[1] + P2.r[1], P1.r[2] + P2.r[2]);}	// Addition
	friend inline Point operator-(const Point& P1, const Point& P2){return Point(P1.r[0] - P2.r[0], P1.r[1] - P2.r[1], P1.r[2] - P2.r[2]);}	// Subtraction
	friend inline Point operator*(const double& d, const Point& P){return Point(d*P.r[0], d*P.r[1], d*P.r[2]);}								// Scalar multiplcation (scalar first)
	friend inline Point operator*(const Point& P, const double& d){return Point(P.r[0]*d, P.r[1]*d, P.r[2]*d);}								// Scalar multiplcation (scalar last)
	friend inline Point operator/(const Point& P, const double& d){	return Point(P.r[0]/d, P.r[1]/d, P.r[2]/d);}							// Scalar division (scalar last)
	friend inline Point operator-(const Point& P){return Point(-P.r[0], -P.r[1], -P.r[2]);}													// Unary minus

	friend ostream& operator<<(ostream &Out, Point &P);			// Formatted output
	friend ostream& operator<<(ostream &Out, const Point &P);	// Formatted output
};


#endif
