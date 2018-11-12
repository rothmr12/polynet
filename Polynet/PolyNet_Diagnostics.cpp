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

	Cycles over all chains in the system to calculate the 9 components of
	the virial true stress tensor of the network. This only takes into account
	the force part of the virial, not the velocity part, because this model does
	not have true velocities.

	(mu, nu) component of the true stress:

		\sigma_{\mu\nu} = -\frac{1}{V} \sum_{chains} \vec{F}_{\mu} \vec{r}_{\nu}

	where V is the domain volume.

*/

void Network::ComputeStressTensor(){

	int i, mu, nu;
	double rmag, F;
	double W[3] = {L[0][1] - L[0][0], L[1][1] - L[1][0], L[2][1] - L[2][0]};
	Point r;

	for(mu = 0; mu < 3; ++mu) for(nu = 0; nu < 3; ++nu)
		Stress[mu][nu] = 0;

	for(auto & ElementC : ChainList){
		if(ElementC.tag != 0) continue;
		r = NodeList[ElementC.Node1].r - NodeList[ElementC.Node0].r;
		if(PeriodicFlag) Periodify(r, W);

		rmag = r.Norm();
		if(rmag <= 0) continue;

		F = Force(rmag, ElementC.r0, ElementC.N);
		for(mu = 0; mu < 3; ++mu) for(nu = 0; nu < 3; ++nu)
			Stress[mu][nu] += (F/rmag) * r.r[mu]*r.r[nu];
	}

	for(mu = 0; mu < 3; ++mu) for(nu = 0; nu < 3; ++nu)
		Stress[mu][nu] /= (W[0]*W[1]*W[2]);

	return;
}





/*

	Calculates the average tortuosity of chains in the network.

*/

double Network::CalculateAverageTortuosity(){

	int NumberOfUnrupturedChains = 0;
	double TortuositySum = 0, W[3] = {L[0][1]-L[0][0], L[1][1]-L[1][0], L[2][1]-L[2][0]};
	Point r;

	for(auto & ElementC : ChainList){
		if(ElementC.tag != 0) continue;
		++NumberOfUnrupturedChains;
		r = NodeList[ElementC.Node1].r - NodeList[ElementC.Node0].r;
		if(PeriodicFlag) Periodify(r, W);

		TortuositySum += ElementC.N * MonomerLength / r.Norm();
	}

	return TortuositySum / NumberOfUnrupturedChains;
}





/*

	Calculates the average end-to-end distance of chains in the network.

*/

double Network::CalculateAverageEndToEndDistance(){

	int NumberOfUnrupturedChains = 0;
	double LengthSum = 0, W[3] = {L[0][1]-L[0][0], L[1][1]-L[1][0], L[2][1]-L[2][0]};
	Point r;

	for(auto & ElementC : ChainList){
		if(ElementC.tag != 0) continue;
		++NumberOfUnrupturedChains;
		r = NodeList[ElementC.Node1].r - NodeList[ElementC.Node0].r;
		if(PeriodicFlag) Periodify(r, W);
		LengthSum += r.Norm();
	}

	return LengthSum / NumberOfUnrupturedChains;
}





/*

	Calculates the total mass of the network in in both chains and nodes.

*/

double Network::GetTotalNetworkMass(){

	double TotalMass = 0;
	for(auto & ElementN : NodeList)
		TotalMass += ElementN.m;

	return TotalMass;
}





/*

	Calculates the total number of monomers in the network.

*/

int Network::GetTotalMonomerCount(){

	int TotalMonomers = 0;
	for(auto & ElementC : ChainList)
		TotalMonomers += ElementC.N;

	return TotalMonomers;
}





/*

	Returns a count of the number of ruptured chains in the system.

*/

int Network::CountRupturedChains(){

	NumRupturedChains = 0;
	for(auto & ElementC : ChainList)
		if(ElementC.tag != 0) ++NumRupturedChains;

	return 0;
}





/*

	Returns the length of the domain's diagonal.

*/

double Network::GetDomainDiagonal(){

	return sqrt( (L[0][1]-L[0][0])*(L[0][1]-L[0][0]) + (L[1][1]-L[1][0])*(L[1][1]-L[1][0]) + (L[2][1]-L[2][0])*(L[2][1]-L[2][0]) );
}
