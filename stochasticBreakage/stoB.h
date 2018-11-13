#ifndef STOB_H
#define STOB_H

/*
Class to stochastically determine chain breakage 
for polyNet
*/
#include <vector>
#include "stdio.h"

class stoB{
 public:
	stoB(bool=false);
	~stoB();

	// Initialization methods
	void setN(int);

	// Set simulation temperature, and timestep
	void setTempDt(double, double);

	// Sets the rate constant for the i^th bond
	void setRate(int, double, double, double);

	// Set PChainMax
	void setPChainMax(double);

	/* Timestep breakage probabilities
	 Breakage probability of a bond is calculated as an integral:
	               s=t
	               /
	 P_i(t) =      |  k_i * exp( - k_i * s ) ds
	               /
	             s=0 
	
	 This method adds k_i * exp(-k_i * t) dt to P_i
	 This is equivalent to using the rectangle
	 rule for integration.

	 The updated values of P_i are then used to update PChain

	*/
 	void timestep(double);

	bool isBroken();

	// Methods only for testing
	void printAll();
	
 private:
	bool broken;
	bool valid;    // Flag that declares that it is OK to timestep
 	double time;   // Maybe unnecessary
 	double dt;     // seconds
 	double T;      // Kelvin
	double beta;   // 1/(kB T)
 	int    N;      // Number of bonds in chain
	int    Nsteps; // Number of timesteps taken

 	// Parameters for each breakable bond
 	// Breakage rate is calculated as
 	// k  = k0 * exp( -beta * Ea );
 	// Ea = Ea0 + EaSlope * F
 	std::vector<double> k0;      // inverse seconds
 	std::vector<double> Ea0;     // electron-Volts
 	std::vector<double> EaSlope; // eV/Newton
 	std::vector<double> PBond;   // dimensionless
	
 	double PChain;          // dimensionless
	double PChainMax;

	int Nassigned;

	// Variables for debugging/testing
	bool debug;
	FILE *fp;
	int iprint;
	
};


#endif
