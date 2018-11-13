#include "stoB.h"
#include <iostream>
#include "stdio.h"
#include "math.h"

stoB::stoB(bool xdebug){
	// Assign defaults
	debug     = xdebug;
	broken    = false;
	valid     = false;
	time      = 0.0;
	dt        = 0.0;
	T         = 0.0;
	N         = -1;
	Nsteps    = 0;
	PChain    = 0.0;
	PChainMax = -1.0;
	Nassigned = 0;
	iprint    = 0;

	if ( debug )
		fp = fopen("timestepping.data","w");
}

stoB::~stoB(){
	if (debug)
		fclose(fp);
}

void stoB::setN(int i){
	N           = i;
	k0.assign     (N,-1.0);
	Ea0.assign    (N,-1.0);
	EaSlope.assign(N,-1.0);
	PBond.assign  (N, 0.0);
}

void stoB::setTempDt(double xT, double xdt){
	if ( N < 0 )
		throw std::invalid_argument
			("*** stoB: Cannot set T before N");
	T    = xT;
	dt   = xdt;
	beta = 1.0/(8.6173303e-5 * T);
}

void stoB::setRate(int i, double xk0, double xEa0, double xEaSlope){
	if ( T < 0.001 )
		throw std::invalid_argument
			("*** stoB: Cannot set rate before T, dt");
	
	if ( Nassigned == N )
		throw std::invalid_argument
			("*** stoB: Cannot set more rates than bonds");
	for( i=0;i<N;i++){
		k0[i]      = xk0;
		Ea0[i]     = xEa0;
		EaSlope[i] = xEaSlope;
		Nassigned += 1;
	}
}

void stoB::setPChainMax(double xP){
	if ( Nassigned != N )
		throw std::invalid_argument
			("*** stoB: Cannot set PChainMax before assigning rates");

	PChainMax = xP;
	valid     = true;
}

void stoB::timestep(double F){
	if ( !valid )
		throw std::invalid_argument
			("*** stoB: Cannot timestep before assigning rates");
	if ( broken ){
		if(debug)
			std::cout << "Broken chain. No timestep for you!";
		return;
	}

	double ki, Ea, prod;
	
	// Recalculate P_i for all bonds, and PChain
	prod = 1;
	for ( int i=0; i<N; i++ ){
		Ea        = Ea0[i] + EaSlope[i]*F;
		ki        = k0[i] * exp ( - beta * Ea );
		PBond[i] += ki * exp(-ki * time)*dt;
		prod     *= (1-PBond[i]);
	}
	PChain = 1-prod;
	
	
	time += dt;
	if ( PChain >= PChainMax ){
		broken = true;
	}
	
	
	if(debug){
		iprint++;
		if (iprint > 10000){
			fprintf(fp, "%e ",time);
			fprintf(fp, "%e ",PChain);
			for (int i=0; i<N; i++)
				fprintf(fp, "%e ",PBond[i]);
			fprintf(fp,"\n");
			iprint = 0;
		}
	}
}

bool stoB::isBroken(){
	return broken;
}

void stoB::printAll(){
	printf("******** Status of stoB object ********\n");
	printf("   Validity = %s\n", valid?"True":"False");
	printf("    # bonds = %4d\n",N);
	printf("Temperature = %7.3f Kelvin\n",T);
	printf("   Timestep = %7.3f seconds\n",dt);
	printf("  PChainMax = %7.3f\n",PChainMax);
	if (valid){
		printf("\n");
		printf(" Bond   k0             Ea0        EaSlope\n");
		for (int i=0; i<N; i++){
			printf("%5d   %10.5e %10.5f   %12.3e\n",i, k0[i], Ea0[i], EaSlope[i]);
		}
	}
	printf("************** End status *************\n");
}
