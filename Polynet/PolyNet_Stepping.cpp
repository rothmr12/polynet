/*

	PolyNet: A mesoscale model of a polymer network.
	This is based on the original EPnet algorithm of David E. Hanson, with some
	modifications.
													John L. Barber,
													jlbarber@lanl.gov
													505-664-0605

*/

#include "PolyNet.h"


/*

	Performs a single time/strain step (Dt) of evolution of the system.

*/

bool Network::Step(){

	int i;
	double DragForce2 = DragForce*DragForce, ForceCoeff = 0.5*Dt/Gamma;
	double W[3] = {L[0][1]-L[0][0], L[1][1]-L[1][0], L[2][1]-L[2][0]};

	// Perform one half timestep of forcing, based on the chain
	// tensions at the beginning of this time step. Movement of
	// a given node occurs here only if the net force on it exceeds
	// the drag force threshold. Movement is via "dashpot" drag force
	// differential equation. (Ask me for details. It is at once both silly and useful.)

	for(auto & ElementN : NodeList){
		if(FreezeEndsFlag && ElementN.tag == 1) continue;
		if(ElementN.F.Norm2() <= DragForce2) continue;
		
		ElementN.r += ForceCoeff * ElementN.F;
		if(PeriodicFlag) for(i = 0; i < 3; ++i){
			while(ElementN.r.r[i] < L[i][0]) ElementN.r.r[i] += W[i];
			while(ElementN.r.r[i] > L[i][1]) ElementN.r.r[i] -= W[i];
		}
	}


	// Perform one full timestep of the applied strain:
	// Every node is moved perfectly according to the nominal strain field.

	if(Deform(Dt))
		return false;


	// Recalculate the net force on each node. (This also checks
	// chains for rupture, and marks them (changes their tag to 1)
	// if necessary.)

	CalculateForces();


	// Perform one half timestep of forcing, as above, based on the chain
	// tensions at the end of the applied strain step.

	for(auto & ElementN : NodeList){
		if(FreezeEndsFlag && ElementN.tag == 1) continue;
		if(ElementN.F.Norm2() <= DragForce2) continue;
		
		ElementN.r += ForceCoeff * ElementN.F;
		if(PeriodicFlag) for(i = 0; i < 3; ++i){
			while(ElementN.r.r[i] < L[i][0]) ElementN.r.r[i] += W[i];
			while(ElementN.r.r[i] > L[i][1]) ElementN.r.r[i] -= W[i];
		}
	}

	// Recalculate the net force on each node. (This also checks
	// chains for rupture, and marks them if necessary.)

	CalculateForces();


	// Finally, update the time to reflect the timestep just performed

	Time += Dt;

	return true;
}





/*

	Applies one time increment of the relevant deformation to the system.

*/

int Network::Deform(double TimeIncrement){

	double LambdaOld[3] = {Lambda[0], Lambda[1], Lambda[2]};

	UpdateLambdas(Time + TimeIncrement);

	for(auto & ElementN : NodeList){
		ElementN.r.r[0] *= Lambda[0] / LambdaOld[0];
		ElementN.r.r[1] *= Lambda[1] / LambdaOld[1];
		ElementN.r.r[2] *= Lambda[2] / LambdaOld[2];
	}

	for(int i = 0; i < 3; ++i){
		L[i][0] = L0[i][0] * Lambda[i];
		L[i][1] = L0[i][1] * Lambda[i];
	}

	return 0;
}





/*

	Changes the three lambdas that make up the diagonals of the engineering strain tensor
	to what their values should be (for the given deformation) at time t.

*/

void Network::UpdateLambdas(double t){

	if(DeformationType == 0){
		Lambda[2] = 1 + StrainRate*t;
		Lambda[0] = Lambda[1] = 1 / sqrt(Lambda[2]);
	}

	else if(DeformationType == 1){
		Lambda[2] = 1 + StrainRate*t;
		if(Lambda[2] > LambdaMax) Lambda[2] = LambdaMax;
		Lambda[0] = Lambda[1] = 1 / sqrt(Lambda[2]);
	}

	else if(DeformationType == 2){
		Lambda[2] = (LambdaMax+1)/2 - cos(2*StrainRate*t/(LambdaMax-1))*(LambdaMax-1)/2 ;
		Lambda[0] = Lambda[1] = 1 / sqrt(Lambda[2]);
	}

	else{
		cout << "Error in UpdateLambdas: Unknown deformation type." << endl;
		exit(1);
	}

	return;
}





/*

	Cycles over all chains to compute the net vector force on each node/crosslink
	in the network.

*/

void Network::CalculateForces(){

	double rmag, F;
	double W[] = {L[0][1]-L[0][0], L[1][1]-L[1][0], L[2][1]-L[2][0]};
	Point r;

	for(auto & ElementN : NodeList)
		ElementN.F.FillWith(0);

	for(auto & ElementC : ChainList){
		if(ElementC.tag != 0) continue;
		r = NodeList[ElementC.Node1].r - NodeList[ElementC.Node0].r;
		if(PeriodicFlag) Periodify(r, W);

		rmag = r.Norm();
		F = Force(rmag, ElementC.r0, ElementC.N);

		if(F > BreakForce){
			ElementC.tag = 1;
			++NumRupturedChains;
			continue;
		}

		r *= F/rmag;
		NodeList[ElementC.Node0].F += r;
		NodeList[ElementC.Node1].F -= r;
	}

	return;
}





/*

	Returns the magnitude of the tension in a chain consisting of NumberOfMonomers
	monomers, with current end-to-end distance r, and with original end-to-end
	distance r0.

	For tortuosity < 1 uses the empirical polynomial for rubber from Hanson and Martin.

*/

double Network::Force(double r, double r0, int NumberOfMonomers){

	if(NumberOfMonomers <= 0)
		return 0;

	if(r <= r0)
		return 0;

	double L0 = NumberOfMonomers * MonomerLength;

	if(r <= L0){
		double FIa = Keps * Temperature * (r - r0) / L0;
		double FIb = Kt * Temperature * r / L0;
		if(FIa <= FIb)
			return FIa;
		return FIb;
	}

	double s = r/L0 - 1;
//	static double C[5] = {8.6815e-9, -107.29e-9, 1151.7e-9, -3107.9e-9, 2567e-9}; /* For isoprene, Newtons per strain */
//	double FII = Kt*Temperature + s*(C[0] + s*(C[1] + s*(C[2] + s*(C[3] + s*C[4]))));
	double FII = Kt*Temperature + s*(8.6815e-9 + s*(-107.29e-9 + s*(1151.7e-9 + s*(-3107.9e-9 + s*2567e-9)))); // Numbers from David Hanson fit

	return FII;
}
