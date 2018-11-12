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

	Loads parameters from the input deck file InputDeckName.

	Also checks loaded parameters, and calculates a few auxiliary parameters.

	Returns true if there is no error in loaded parameters (e.g. bad values),
	otherwise returns false.

*/

bool Network::LoadNetworkParameters(){

	ifstream ParameterStream(InputDeckName.c_str());
	if(!ParameterStream.is_open()){
		cout << "Error in Network::LoadNetworkParameters: Unable to access input deck file " << InputDeckName << endl;
		return false;
	}

	bool ReturnValue = true;
	int TempInt;
	double TempDouble, TempDoubleArray[3];
	string TempString;

	GetParameter(ParameterStream, "Seed", TempInt, Seed_Default);
	if(TempInt <= 0) Seed = long(time(0));
	else             Seed = long(TempInt);

	GetParameter(ParameterStream, "LoadNetworkFlag", LoadNetworkFlag, LoadNetworkFlag_Default);
	if(LoadNetworkFlag)
		GetParameter(ParameterStream, "NetworkFileName", NetworkFileName, NetworkFileName_Default);


	// System geometry
	double Width[3];
	GetParameters(ParameterStream, "InitialWidth", Width, 3, InitialDomainWidth_Default);
	if(Width[0] <= 0 || Width[1] <= 0 || Width[2] <= 0){
		cout << "Error in Network::LoadNetworkParameters: ";
		cout << Width[0] << " x " << Width[1] << " x " << Width[2];
		cout << " is an invalid set of initial domain dimensions." << endl;
		ReturnValue = false;
	}
	else{
		for(int i = 0; i < 3; ++i){
			L0[i][0] = L[i][0] = -Width[i]/2;
			L0[i][1] = L[i][1] = +Width[i]/2;
		}
	}
	GetParameter(ParameterStream, "PeriodicFlag", PeriodicFlag, PeriodicFlag_Default);


	// Material properties
	GetParameter(ParameterStream, "MonomerNumberDensity", MonomerNumberDensity, MonomerNumberDensity_Default);
	if(MonomerNumberDensity <= 0){
		cout << "Error in Network::LoadNetworkParameters: " << MonomerNumberDensity << " is an invalid monomer number density." << endl;
		ReturnValue = false;
	}

	GetParameter(ParameterStream, "NodeNumberDensity", NodeNumberDensity, NodeNumberDensity_Default);
	if(NodeNumberDensity <= 0){
		cout << "Error in Network::LoadNetworkParameters: " << NodeNumberDensity << " is an invalid node/crosslink number density." << endl;
		ReturnValue = false;
	}
	else if(NodeNumberDensity >= MonomerNumberDensity/2){
		cout << "Error in Network::LoadNetworkParameters: NodeNumberDensity (" << NodeNumberDensity << ") must be < MonomerNumberDensity/2 (" << MonomerNumberDensity/2 << ")" << endl;
		ReturnValue = false;
	}
	pxL = 2*NodeNumberDensity / MonomerNumberDensity;

	NumNodes = int( NodeNumberDensity * Width[0] * Width[1] * Width[2] + 0.5 );

	GetParameter(ParameterStream, "MassNode", MassNode, MassNode_Default);
	if(MassNode < 0) MassNode = 0;

	GetParameter(ParameterStream, "MassMonomer", MassMonomer, MassMonomer_Default);
	if(MassMonomer < 0) MassMonomer = 0;

	GetParameter(ParameterStream, "MonomerLength", MonomerLength, MonomerLength_Default);
	if(MonomerLength <= 0){
		cout << "Error in Network::LoadNetworkParameters: " << MonomerLength << " is an invalid monomer length." << endl;
		ReturnValue = false;
	}

	GetParameter(ParameterStream, "KuhnLength", KuhnLength, KuhnLength_Default);
	if(KuhnLength <= 0){
		cout << "Error in Network::LoadNetworkParameters: " << KuhnLength << " is an invalid Kuhn length." << endl;
		ReturnValue = false;
	}

	GetParameter(ParameterStream, "ChainsPerNode", ChainsPerNode, ChainsPerNode_Default);
	if(ChainsPerNode <= 0){
		cout << "Error in Network::LoadNetworkParameters: " << ChainsPerNode << " is an invalid coordination number." << endl;
		ReturnValue = false;
	}

	GetParameter(ParameterStream, "RadiusProbabilityCutoff", TempDouble, RadiusProbabilityCutoff_Default);
	if(TempDouble <= 0 || TempDouble >= 1){
		cout << "Error in Network::LoadNetworkParameters: " << TempDouble << " is an invalid probability cutoff radius." << endl;
		ReturnValue = false;
	}
	else
		MaximumRadius = MaxRadius(TempDouble);


	// Initialize current state
	Time = 0;
	Lambda[0] = Lambda[1] = Lambda[2] = 1;


	// Dynamics
	GetParameter(ParameterStream, "Temperature", Temperature, Temperature_Default);
	if(Temperature <= 0){
		cout << "Error in Network::LoadNetworkParameters: " << Temperature << " is an invalid temperature." << endl;
		ReturnValue = false;
	}

	GetParameter(ParameterStream, "Keps", Keps, Keps_Default);
	if(Keps <= 0){
		cout << "Error in Network::LoadNetworkParameters: " << Keps << " is an invalid value for Keps." << endl;
		ReturnValue = false;
	}

	GetParameter(ParameterStream, "Kt", Kt, Kt_Default);
	if(Kt <= 0){
		cout << "Error in Network::LoadNetworkParameters: " << Kt << " is an invalid value for Kt." << endl;
		ReturnValue = false;
	}

	GetParameter(ParameterStream, "DragForce", DragForce, DragForce_Default);
	if(DragForce <= 0){
		cout << "Error in Network::LoadNetworkParameters: " << DragForce << " is an invalid drag force." << endl;
		ReturnValue = false;
	}

	GetParameter(ParameterStream, "BreakForce", BreakForce, BreakForce_Default);
	if(BreakForce <= 0){
		cout << "Error in Network::LoadNetworkParameters: " << BreakForce << " is an invalid breaking force." << endl;
		ReturnValue = false;
	}

	GetParameter(ParameterStream, "Gamma", Gamma, Gamma_Default);
	if(Gamma <= 0){
		cout << "Error in Network::LoadNetworkParameters: " << Gamma << " is an invalid value for gamma." << endl;
		ReturnValue = false;
	}

	GetParameter(ParameterStream, "DeformationType", DeformationType, DeformationType_Default);
	if(DeformationType < 0 || DeformationType > 2){
		cout << "Error in Network::LoadNetworkParameters: " << DeformationType << " is an unknown deformation type." << endl;
		cout << "0 = Simple Strain, 1 = Stress relaxation, 2 = Oscillation." << endl;
		ReturnValue = false;
	}

	GetParameter(ParameterStream, "StrainRate", StrainRate, StrainRate_Default);

	GetParameter(ParameterStream, "LambdaMax", LambdaMax, LambdaMax_Default);
	if(LambdaMax <= 0){
		cout << "Error in Network::LoadNetworkParameters: " << LambdaMax << " is an invalid max lambda." << endl;
		ReturnValue = false;
	}

	if(DeformationType != 0)
		GetParameter(ParameterStream, "TimeMax", TimeMax, TimeMax_Default);

	if(DeformationType == 0){ // Simple strain
		if(LambdaMax-1 >= 0)
			StrainRate = +fabs(StrainRate);
		else
			StrainRate = -fabs(StrainRate);
		if(StrainRate * (LambdaMax-1) == 0)
			TimeMax = 0;
		else
			TimeMax = (LambdaMax-1)/StrainRate;
	}
	else if(DeformationType == 1){	// Stress relaxation
		if(LambdaMax-1 > 0)
			StrainRate = +fabs(StrainRate);
		else if(LambdaMax-1 < 0)
			StrainRate = -fabs(StrainRate);
		if(StrainRate*(LambdaMax-1) > 0 && TimeMax <= (LambdaMax-1)/StrainRate)
			cout << "Warning in Network::LoadNetworkParameters: TimeMax (" << TimeMax << ") is insufficient to exceed the time (" << (LambdaMax-1)/StrainRate << ") of  max lambda." << endl;
	}

	GetParameter(ParameterStream, "DtCoefficient", TempDouble, DtCoefficient_Default);
	if(TempDouble <= 0){
		cout << "Error in Network::LoadNetworkParameters: " << TempDouble << " is an invalid time step coefficient." << endl;
		ReturnValue = false;
	}
	double AverageDistanceBetweenNodes = pow(NodeNumberDensity, -1./3.);
	double TVisc = AverageDistanceBetweenNodes*Gamma / DragForce;
	double TEps = 1 / fabs(StrainRate);
	Dt = TempDouble * (TVisc <= TEps ? TVisc : TEps);
	NumberOfTimeSteps = int(ceil(TimeMax / Dt));
	Dt = TimeMax / NumberOfTimeSteps;

	GetParameter(ParameterStream, "FreezeEndsFlag", FreezeEndsFlag, FreezeEndsFlag_Default);
	GetParameter(ParameterStream, "FreezeFraction", FreezeFraction, FreezeFraction_Default);
	if(FreezeFraction < 0) FreezeFraction = 0;
	if(FreezeFraction > 0.5) FreezeFraction = 0.5;


	// Output
	GetParameter(ParameterStream, "OutputFolder", OutputFolder, OutputFolder_Default);
	if(OutputFolder[OutputFolder.size() - 1] != '/')
		OutputFolder = OutputFolder + "/";

	GetParameter(ParameterStream, "StrainStressFileFlag", StrainStressFileFlag, StrainStressFileFlag_Default);

	if(StrainStressFileFlag){
		GetParameter(ParameterStream, "StrainStressFile", StrainStressFile, StrainStressFile_Default);
		StrainStressFile = OutputFolder + StrainStressFile;
		StrainStressFileStream.open(StrainStressFile.c_str());
		if(!StrainStressFileStream.is_open()){
			cout << "Error in Network::LoadNetworkParameters: Unable to initialize strain/stress file " << StrainStressFile << endl;
			ReturnValue = false;
		}
	}

	GetParameter(ParameterStream, "DrawNetworkImagesFlag", DrawNetworkImagesFlag, DrawNetworkImagesFlag_Default);

	if(DrawNetworkImagesFlag){
		GetParameter(ParameterStream, "NetworkImageInterval", NetworkImageInterval, NetworkImageInterval_Default);

		GetParameter(ParameterStream, "NumberOfPixels", NumberOfPixels, NumberOfPixels_Default);
		if(NumberOfPixels <= 0){
			cout << "Error in Network::LoadNetworkParameters: " << NumberOfPixels << " is an invalid number of image pixels." << endl;
			ReturnValue = false;
		}

		GetParameter(ParameterStream, "NetworkImageName", NetworkImageName, NetworkImageName_Default);
		NetworkImageName = OutputFolder + NetworkImageName;
		GetParameter(ParameterStream, "ImageFormat", ImageFormat, ImageFormat_Default);

		if(!GetParameters(ParameterStream, "ViewPoint", ViewPoint.r, 3))
			ViewPoint = Point(ViewPointX_Default, ViewPointY_Default, ViewPointZ_Default);
		ViewPoint = GetDomainDiagonal() * ViewPoint;
		GetParameter(ParameterStream, "RescaleViewPointFlag", RescaleViewPointFlag, RescaleViewPointFlag_Default);

		GetParameters(ParameterStream, "ViewDirection", ViewDirection.r, 3, 0);
		if(ViewDirection.Norm2() <= 0){
			if(ViewPoint.Norm2() <= 0)
				ViewDirection = Point(ViewDirectionX_Default, ViewDirectionY_Default, ViewDirectionZ_Default);
			else
				ViewDirection = -ViewPoint;
		}
		ViewDirection.Normalize();

		GetParameters(ParameterStream, "Up", Up.r, 3, 0);
		if(Up.Norm2() <= 0)
			Up = Point(UpX_Default, UpY_Default, UpZ_Default);
		Up.Normalize();

		GetParameters(ParameterStream, "LightDirection", LightDirection.r, 3, 0);
		if(LightDirection.Norm2() <= 0)
			LightDirection = -Up;
		else
			LightDirection.Normalize();

		GetParameters(ParameterStream, "AngularSpan", TempDoubleArray, 2, AngularSpan_Default);
		AngularSpanX = (TempDoubleArray[0] <= 0 ? AngularSpan_Default : TempDoubleArray[0]);
		AngularSpanY = (TempDoubleArray[1] <= 0 ? AngularSpan_Default : TempDoubleArray[1]);

		GetParameter(ParameterStream, "NodeRadius", NodeRadius, NodeRadius_Default);
		if(NodeRadius < 0) NodeRadius = 0;
		NodeRadius *= MonomerLength;

		GetParameter(ParameterStream, "ChainRadius", ChainRadius, ChainRadius_Default);
		if(ChainRadius < 0) ChainRadius = 0;
		ChainRadius *= MonomerLength;

		GetParameter(ParameterStream, "BorderRadius", BorderRadius, BorderRadius_Default);
		if(BorderRadius < 0) BorderRadius = 0;
		BorderRadius *= MonomerLength;

		GetColor(ParameterStream, "BackgroundColor", BackgroundColor, BackgroundColor_Default);
		GetColor(ParameterStream, "Node0Color", Node0Color, Node0Color_Default);
		GetColor(ParameterStream, "Node1Color", Node1Color, Node1Color_Default);
		GetColor(ParameterStream, "LooseChainColor", LooseChainColor, LooseChainColor_Default);
		GetColor(ParameterStream, "TightChainColor", TightChainColor, TightChainColor_Default);
		GetColor(ParameterStream, "BorderColor", BorderColor, BorderColor_Default);
	}

	ParameterStream.close();

	return ReturnValue;
}





/*

	Loads a network of nodes and chains from the specified file, and in the
	specified format.

	This was made to load an output file from Dave Hanson's EPnet code, and in EPnet's format.
	I don't know if this will ever be used again, but here it is.

*/

bool Network::LoadNetwork(string Format /* = "EPNet" */){

	if(Format != "EPNet"){
		cout << "Error in LoadNetwork: " << Format << " is an unknown file format." << endl;
		cout.flush();
		return false;
	}

	int i, j;
	Node TempNode;
	Chain TempChain;

	NodeList.clear();
	ChainList.clear();
	NumNodes = NumChains = 0;

	ifstream NetStream(NetworkFileName.c_str());
	if(!NetStream.is_open()){
		cout << endl << "Error in LoadNetwork: Unable to access the file " << NetworkFileName << endl;
		cout.flush();
		return false;
	}


	/* EPNet */
	PeriodicFlag = true;

	char NetStreamLine[1024];
	char Delimiter = '\n';
	double Lx, Ly, Lz;

	NetStream.clear();
	NetStream.seekg(0, ios::beg);

	NetStream.getline(NetStreamLine, 1024, Delimiter); // Get rid of header line
	NetStream >> NetStreamLine; // Get rid of 's';

	NetStream >> NumNodes;
	NetStream >> NumChains;
	NetStream >> MonomerLength;
	NetStream >> KuhnLength;
	NetStream >> Lambda[1];
	Lambda[1] += 1;
	Lambda[0] = Lambda[2] = 1 / sqrt(Lambda[1]);
	NetStream >> Lx;
	NetStream >> Ly;
	NetStream >> Lz;

	L[0][0] = L0[0][0] = -Lx/2;
	L[0][1] = L0[0][1] = +Lx/2;
	L[1][0] = L0[1][0] = -Ly/2;
	L[1][1] = L0[1][1] = +Ly/2;
	L[2][0] = L0[2][0] = -Lz/2;
	L[2][1] = L0[2][1] = +Lz/2;
	double W[3] = {L[0][1]-L[0][0], L[1][1]-L[1][0], L[2][1]-L[2][0]};

		// Get nodes
	TempNode.m = MassNode;
	TempNode.NumChains = 0;
	TempNode.tag = 0;
	TempNode.r.r[0] = TempNode.r.r[1] = TempNode.r.r[2] = 0;
	TempNode.F.r[0] = TempNode.F.r[1] = TempNode.F.r[2] = 0;

	try{
		NodeList.insert(NodeList.end(), NumNodes, TempNode);
	}
	catch(bad_alloc& ba){
		cout << "Error in Network::LoadNetwork: Unable to allocate memory to hold " << NumNodes << " nodes." << endl;
		return false;
	}

	for(i = 0; i < NumNodes; ++i){
		for(j = 0; j < 3; ++j){
			NetStream >> NodeList[i].r.r[j];
			while(NodeList[i].r.r[j] <  L[j][0]) NodeList[i].r.r[j] += L[j][1] - L[j][0];
			while(NodeList[i].r.r[j] >= L[j][1]) NodeList[i].r.r[j] -= L[j][1] - L[j][0];
		}
	}


		// Get chains
	int Node0, Node1, Status, TotalMonomers = 0;
	double Tortuosity, Tension, Length;
	Point rr;

	TempChain.Node0 = 0;
	TempChain.Node1 = 0;
	TempChain.r0 = 0;
	TempChain.N = 0;
	TempChain.tag = 0;

	MaximumRadius = -DBL_MAX;

	int NumberOfIntactChains = 0;

	for(i = 0; i < NumChains; ++i){
		NetStream >> Node0;
		NetStream >> Node1;
		NetStream >> Tortuosity;
		NetStream >> Tension;
		NetStream >> Status;

		if(Status == 0) continue;

		++NumberOfIntactChains;

		if(Node1 >= Node0){
			TempChain.Node0 = Node0 - 1;
			TempChain.Node1 = Node1 - 1;
		}
		else{
			TempChain.Node0 = Node1 - 1;
			TempChain.Node1 = Node0 - 1;
		}

		rr = NodeList[TempChain.Node1].r - NodeList[TempChain.Node0].r;
		if(PeriodicFlag) Periodify(rr, W);

		TempChain.r0 = rr.Norm();
		if(TempChain.r0 > MaximumRadius) MaximumRadius = TempChain.r0;

		Length = TempChain.r0 * Tortuosity;
		TempChain.N = int(Length/MonomerLength + 0.5);

		TotalMonomers += TempChain.N;

		try{
			ChainList.push_back(TempChain);
		}
		catch(bad_alloc& ba){
			cout << "Error in LoadNetwork: Unable to allocate memory to hold chain " << i+1 << " of " << NumChains << endl;
			cout.flush();
			return false;
		}

		NodeList[TempChain.Node0].NumChains++;
		NodeList[TempChain.Node1].NumChains++;
		NodeList[TempChain.Node0].m += TempChain.N*MassMonomer / 2;
		NodeList[TempChain.Node1].m += TempChain.N*MassMonomer / 2;
	}

	NetStream.close();

	NumChains = NumberOfIntactChains;
	pxL = 2.*NumChains/(1.*TotalMonomers);

	MarkUnconnectedNodes();

	return true;
}





/*

	Removes from the network (i.e. sets their tag equal to 1)
	those nodes which are connected to zero chains.

*/

int Network::MarkUnconnectedNodes(){

	int NumberMarked = 0;

	for(auto & ElementN : NodeList){
		if(ElementN.NumChains > 0) continue;
		ElementN.tag = 1;
		++NumberMarked;
	}

	return NumberMarked;
}





/*

	Generates a random set of NumNodes node positions in the current domain.
	Erases any existing nodes.

*/

bool Network::CreateNodes(){

	clock_t c_begin = clock();
	cout << "Generating nodes... " << endl;
	cout.flush();

	int NumberFrozen = 0;

	// Make dummy Node and fill this Network's NodeList vector with
	// NumNodes copies of it.
	Node TempNode;
	TempNode.r.FillWith(0);
	TempNode.F.FillWith(0);
	TempNode.m = MassNode;
	TempNode.T = Temperature;
	TempNode.NumChains = 0;
	TempNode.tag = 0;

	NodeList.clear();

	try{
		NodeList.insert(NodeList.end(), NumNodes, TempNode);
	}
	catch(bad_alloc& ba){
		cout << "Error in Network::CreateNodes: Unable to allocate memory to hold " << NumNodes << " nodes." << endl;
		return false;
	}

	// Now, run through list of nodes and give each a random position in the box (uniformly)
	double Lx = L[0][1] - L[0][0];
	double Ly = L[0][1] - L[0][0];
	double Lz = L[0][1] - L[0][0];
	for(auto & ElemenN : NodeList){
		ElemenN.r = Point(L[0][0]+rand1(Seed)*Lx, L[1][0]+rand1(Seed)*Ly, L[2][0]+rand1(Seed)*Lz);

		if(FreezeEndsFlag){
			if(ElemenN.r.r[2] <= L[2][0] + FreezeFraction*Lz || ElemenN.r.r[2] >= L[2][1] - FreezeFraction*Lz){
				ElemenN.tag = 1;
				++NumberFrozen;
			}
		}
	}

	cout << "\t" << NumNodes << " nodes created." << endl;
	if(FreezeEndsFlag)
		cout << "\tNumber frozen: " << NumberFrozen << endl;
	cout << "\t" << ((double)(clock() - c_begin)) / CLOCKS_PER_SEC << " seconds" << endl;
	cout.flush();

	return true;
}





/*

	Generates a set of chains between nodes in the network.

*/

bool Network::CreateChains(){

	clock_t c_begin = clock();
	cout << "Generating chains... " << endl;
	cout.flush();

	vector<Chain> TrialChains;
	vector<double> Randoms;
	Chain TempChain;
	int i, j, k, n0, n1, kmax, Ntrial, ChainsNeeded, ChainsSoFar = 0;
	double d01, Weight, Rand, RandMax;
	Point rr;
	r2overNbm = 0;


	// Clear out any old chain structure (if any) to start anew
	ChainList.clear();
	for(auto & ElementN : NodeList)
		ElementN.NumChains = 0;


	// See if the system is large enough to require a cell structure to avoid O(N^2) timing
	double Lx = L[0][1] - L[0][0];
	double Ly = L[1][1] - L[1][0];
	double Lz = L[2][1] - L[2][0];
	double W[3] = {L[0][1]-L[0][0], L[1][1]-L[1][0], L[2][1]-L[2][0]};
	int ncx = int( Lx / MaximumRadius );
	int ncy = int( Ly / MaximumRadius );
	int ncz = int( Lz / MaximumRadius );
	double Dx = Lx / ncx;
	double Dy = Ly / ncy;
	double Dz = Lz / ncz;
	bool CellsNeededFlag = false;
	if((ncx >= 3) && (ncy >= 3) && (ncz >= 3))		// At least 3 cells in each direction
		if((ncx >= 4) || (ncy >= 4) || (ncz >= 4))	// and at least one direction with >= 4 cells.
			CellsNeededFlag = true;


	// System small. Do it the brute force way.
	if(!CellsNeededFlag){
		cout << "\tNo cell structure." << endl;
		cout.flush();

		cout << "\tNode ";
		for(i = 0; i < NumNodes; ++i){
			CountdownWithSlash(i+1, NumNodes, !(i == 0));
			ChainsNeeded = ChainsPerNode - NodeList[i].NumChains;
			if(ChainsNeeded <= 0) continue;

			RandMax = -DBL_MAX;
			Randoms.clear();
			TrialChains.clear();

			TempChain.Node0 = i;
			TempChain.tag = 0;

			for(j = i+1; j < NumNodes; ++j){
				if(NodeList[j].NumChains >= ChainsPerNode) continue;
				rr = NodeList[i].r - NodeList[j].r;
				if(PeriodicFlag) Periodify(rr, W);

				d01 = rr.Norm();
				if(d01 > MaximumRadius) continue;

				Weight = PDF_r(d01);
				if(Weight <= 0) continue;

				TempChain.Node1 = j;
				TempChain.r0 = d01;
				Rand = -log(1 - rand1(Seed)) / Weight;

				if(TrialChains.size() < ChainsNeeded){
					TrialChains.push_back(TempChain);
					Randoms.push_back(Rand);
					if(Rand > RandMax){
						RandMax = Rand;
						kmax = TrialChains.size() - 1;
					}
					continue;
				}

				if(Rand >= RandMax) continue;

				TrialChains[kmax] = TempChain;
				Randoms[kmax] = Rand;
				RandMax = -DBL_MAX;
				for(k = 0; k < TrialChains.size(); ++k){
					if(Randoms[k] > RandMax){
						RandMax = Randoms[k];
						kmax = k;
					}
				}
			}

			for(j = 0; j < TrialChains.size(); ++j){

				do{
					Ntrial = GetNumberOfSteps(TrialChains[j].r0);
				} while(Ntrial*MonomerLength < TrialChains[j].r0);
				TrialChains[j].N = Ntrial;
				++ChainsSoFar;
				try{
					ChainList.push_back(TrialChains[j]);
				}
				catch(bad_alloc& ba){
					cout << "Error in CreateChains: Unable to allocate memory to hold chain " << ChainsSoFar << endl;
					cout.flush();
					return false;
				}

				n0 = TrialChains[j].Node0;
				n1 = TrialChains[j].Node1;
				NodeList[n0].m += 0.5*TrialChains[j].N*MassMonomer;
				NodeList[n1].m += 0.5*TrialChains[j].N*MassMonomer;
				++NodeList[n0].NumChains;
				++NodeList[n1].NumChains;

				r2overNbm += (TrialChains[j].r0*TrialChains[j].r0) / (TrialChains[j].N*MonomerLength);
			}
		}

	}


	// Bigger system. Use cell structure.

	else{

		cout << "\tCell structure: " << ncx << " x " << ncy << " x " << ncz << endl;
		cout.flush();

		int ix, iy, iz, ix1, iy1, iz1, n, nn;

		// Make cell structure: A list of nodes for each cell (CellNode), A trio of cell indices for each cell (CellIndex),
		// and a cell number for each node (NodeCell).

		vector<int> *CellNode = new(nothrow) vector<int>[ncx*ncy*ncz];
		if(!CellNode){
			cout << "Error in Network::CreateChains: Unable to allocate space for CellNode array." << endl;
			return false;
		}

		int *NodeCell = new(nothrow) int[NumNodes];
		if(!NodeCell){
			cout << "Error in Network::CreateChains: Unable to allocate space for NodeCell array." << endl;
			return false;
		}

		for(i = 0; i < NumNodes; ++i){
			ix = int( (NodeList[i].r.r[0] - L[0][0]) / Dx );
			iy = int( (NodeList[i].r.r[1] - L[1][0]) / Dy );
			iz = int( (NodeList[i].r.r[2] - L[2][0]) / Dz );
			CellNode[iz + ncz*(iy + ncy*ix)].push_back(i);
			NodeCell[i] = iz + ncz*(iy + ncy*ix);
		}

		int CellIndex[ncx*ncy*ncz][3];

		for(ix = 0; ix < ncx; ++ix) for(iy = 0; iy < ncy; ++iy) for(iz = 0; iz < ncz; ++iz){
			CellIndex[iz + ncz*(iy + ncy*ix)][0] = ix;
			CellIndex[iz + ncz*(iy + ncy*ix)][1] = iy;
			CellIndex[iz + ncz*(iy + ncy*ix)][2] = iz;
		}

		// Cycle through nodes
		cout << "\tNode ";
		for(i = 0; i < NumNodes; ++i){
			CountdownWithSlash(i+1, NumNodes, !(i == 0));

			ChainsNeeded = ChainsPerNode - NodeList[i].NumChains;
			if(ChainsNeeded <= 0) continue;

			RandMax = -DBL_MAX;
			Randoms.clear();
			TrialChains.clear();

			TempChain.Node0 = i;
			TempChain.tag = 0;

			// Cycle through cells around this node
			for(ix = -1; ix <= +1; ++ix) for(iy = -1; iy <= +1; ++iy) for(iz = -1; iz <= +1; ++iz){
				ix1 = CellIndex[NodeCell[i]][0] + ix;
				iy1 = CellIndex[NodeCell[i]][1] + iy;
				iz1 = CellIndex[NodeCell[i]][2] + iz;
				if(!PeriodicFlag && (ix1 < 0 || ix1 >= ncx || iy1 < 0 || iy1 >= ncy || iz1 < 0 || iz1 >= ncz)) continue;
				if(ix1 < 0) ix1 += ncx;
				else if(ix1 >= ncx) ix1 -= ncx;
				if(iy1 < 0) iy1 += ncy;
				else if(iy1 >= ncy) iy1 -= ncy;
				if(iz1 < 0) iz1 += ncz;
				else if(iz1 >= ncz) iz1 -= ncz;

				// Cycle through nodes in the current neighbor cell
				nn = CellNode[iz1 + ncz*(iy1 + ncy*ix1)].size();
				for(n = 0; n < nn; ++n){
					j = CellNode[iz1 + ncz*(iy1 + ncy*ix1)][n]; 
					if(j <= i) continue;

					if(NodeList[j].NumChains >= ChainsPerNode) continue;
					rr = NodeList[i].r - NodeList[j].r;
					if(PeriodicFlag) Periodify(rr, W);

					d01 = rr.Norm();
					if(d01 > MaximumRadius) continue;

					Weight = PDF_r(d01);
					if(Weight <= 0) continue;

					TempChain.Node1 = j;
					TempChain.r0 = d01;
					Rand = -log(1 - rand1(Seed)) / Weight;

					if(TrialChains.size() < ChainsNeeded){
						TrialChains.push_back(TempChain);
						Randoms.push_back(Rand);
						if(Rand > RandMax){
							RandMax = Rand;
							kmax = TrialChains.size() - 1;
						}
						continue;
					}

					if(Rand >= RandMax) continue;

					TrialChains[kmax] = TempChain;
					Randoms[kmax] = Rand;
					RandMax = -DBL_MAX;
					for(k = 0; k < TrialChains.size(); ++k){
						if(Randoms[k] > RandMax){
							RandMax = Randoms[k];
							kmax = k;
						}
					}
				} // n: Neighboring node loop
			} // ix, iy, iz: Cell loop

			for(j = 0; j < TrialChains.size(); ++j){

				do{
					Ntrial = GetNumberOfSteps(TrialChains[j].r0);
				} while(Ntrial*MonomerLength < TrialChains[j].r0);
				TrialChains[j].N = Ntrial;

				++ChainsSoFar;
				try{
					ChainList.push_back(TrialChains[j]);
				}
				catch(bad_alloc& ba){
					cout << "Error in CreateChains: Unable to allocate memory to hold chain " << ChainsSoFar << endl;
					cout.flush();
					return false;
				}

				n0 = TrialChains[j].Node0;
				n1 = TrialChains[j].Node1;
				NodeList[n0].m += TrialChains[j].N*MassMonomer / 2;
				NodeList[n1].m += TrialChains[j].N*MassMonomer / 2;
				NodeList[n0].NumChains++;
				NodeList[n1].NumChains++;

				r2overNbm += (TrialChains[j].r0*TrialChains[j].r0) / (TrialChains[j].N*MonomerLength);
			}

		} // i: Main node loop

		cout << endl;

		delete [] CellNode;
		delete [] NodeCell;

	}

	NumChains = ChainList.size();
	NumRupturedChains = 0;
	r2overNbm /= NumChains;

	cout << "  " << NumChains << " chains created." << endl;
	cout << "\tAverage value of r^2 / N b_m = " << r2overNbm << endl;
	cout << "\t" << ((double)(clock() - c_begin)) / CLOCKS_PER_SEC << " seconds" << endl;
	cout.flush();

	return true;
}





/*

	Print the parameters that characterize the network.

*/

void Network::PrintNetworkParameters(){

	int* ChainCountFrequencies = new int[ChainsPerNode+1];
	int i, TotalNodeCount = 0;
	for(i = 0; i <= ChainsPerNode; i++)
		ChainCountFrequencies[i] = 0;
	for(i = 0; i < NumNodes; i++){
		ChainCountFrequencies[NodeList[i].NumChains]++;
		TotalNodeCount++;
	}
	cout << "Number of nodes: " << NumNodes << endl;
	cout << "Nominal coordination number: " << ChainsPerNode << endl; 
	for(i = 0; i <= ChainsPerNode; i++)
		cout << "\t" << ChainCountFrequencies[i] << " nodes connected to " << i << " chains" << endl;
	cout << "\t" << "Total: " << TotalNodeCount << endl;
	delete [] ChainCountFrequencies;

	cout << "Number of chains: " << NumChains << endl;
	cout << "< r^2 / N b_m > = " << r2overNbm << " = " << r2overNbm / KuhnLength << " x nominal." << endl;

	cout << "Random number seed: " << Seed << endl;
	if(LoadNetworkFlag)
		cout << "Network loaded from file " << NetworkFileName << endl;
	if(PeriodicFlag)
		cout << "Periodic domain, ";
	else
		cout << "Non-periodic domain, ";
	cout << "[" << L0[0][0] << ", " << L0[0][1] << "] x [" << L0[1][0] << ", " << L0[1][1] << "] x [" << L0[2][0] << ", " << L0[2][1] << "]" << endl;
	if(FreezeEndsFlag)
		cout << "A fraction " << FreezeFraction << " of nodes are frozen at each end of the domain." << endl;

	// Material properties
	cout << "Monomer number density: " << MonomerNumberDensity << endl;
	cout << "Monomer mass: " << MassNode << endl;
	cout << "Monomer length (avg): " << MonomerLength << endl;
	cout << "Kuhn length: " << KuhnLength << endl;
	cout << "Node/Crosslink number density: " << NodeNumberDensity << endl;
	cout << "Node/Crosslink mass: " << MassNode << endl;
	cout << "pxL: " << pxL << endl;
	cout << "Maximum end-to-end distance = " << MaximumRadius << " = " << MaximumRadius/MonomerLength << " x monomer length." << endl;

	// Dynamics
	cout << "Temperature: " << Temperature << endl;
	cout << "Keps: " << Keps << endl;
	cout << "Kt: " << Kt << endl;
	cout << "DragForce: " << DragForce << endl;
	cout << "BreakForce: " << BreakForce << endl;
	cout << "Friction coefficient gamma: " << Gamma << endl;

	if(DeformationType == 0){
		cout << "Simple strain." << endl;
		cout << "Strain rate: " << StrainRate << endl;
		cout << "Ultimate extension (lambda): " << LambdaMax << endl;
	}
	else if(DeformationType == 1){
		cout << "Stress relaxation." << endl;
		cout << "Strain rate: " << StrainRate << endl;
		cout << "Maximum strain of " << LambdaMax - 1 << " obtains at t = " << (LambdaMax - 1)/StrainRate << " = " << (LambdaMax - 1)/(StrainRate*Dt) << " timesteps." << endl;
	}
	else if(DeformationType == 2){
		cout << "Oscillation." << endl;
		cout << "Maximum strain rate: " << StrainRate << endl;
		cout << "Strain oscillates between 0 and " << LambdaMax - 1 << endl;
		double OscillationPeriod = 2*PI / (2*StrainRate/(LambdaMax-1));
		cout << "Oscillation period: " << OscillationPeriod << " = " << OscillationPeriod / Dt << " timesteps." << endl;
	}

	cout << "Timestep: " << Dt << endl;
	cout << "Maximum simulation time: " << TimeMax << " = " << NumberOfTimeSteps << " timesteps" << endl;

	// Output
	if(StrainStressFileFlag || DrawNetworkImagesFlag){
		cout << endl;
		cout << "Output directory: " << OutputFolder << endl;
		if(StrainStressFileFlag)
			cout << "Strain/stress file: " << StrainStressFile << endl;
		if(DrawNetworkImagesFlag){
			if(NetworkImageInterval <= 0)
				cout << "An image " << NetworkImageName << ".0." << ImageFormat << " of the network's initial condition will be created." << endl;
			else
				cout << "An image " << NetworkImageName << ".[timestep]." << ImageFormat << " of the network will be created every " << NetworkImageInterval << " timesteps." << endl;
			cout << "ViewPoint:      " << ViewPoint << endl;
			cout << "ViewDirection:  " << ViewDirection << endl;
			cout << "Up orientation: " << Up << endl;
			cout << "LightDirection: " << LightDirection << endl;
			cout << "Angular dimensions: " << AngularSpanX << " x " << AngularSpanY << endl;

			cout << "Image size: " << NumberOfPixels << " pixels." << endl;
			cout << "Colors (R G B): " << endl;
			cout << "\tBackground:    " << ColorToString(BackgroundColor) << endl;
			cout << "\tNodes (Tag 0): " << ColorToString(Node0Color) << endl;
			cout << "\tNodes (Tag 1): " << ColorToString(Node1Color) << endl;
			cout << "\tLoose Chains:  " << ColorToString(LooseChainColor) << endl;
			cout << "\tTight Chains:  " << ColorToString(TightChainColor) << endl;
			cout << "\tBorder:        " << ColorToString(BorderColor) << endl;
		}
	}

	cout << endl;

	// Calculate some statistical network diagnostics

	double DomainVolume = (L[0][1] - L[0][0]) * (L[1][1] - L[1][0]) * (L[2][1] - L[2][0]);
	double p = pxL;
	double ExpectedDensity = MonomerNumberDensity*MassMonomer + NodeNumberDensity*MassNode;
	double ActualDensity = GetTotalNetworkMass() / DomainVolume;
	double ExpectedAvgEndToEnd = sqrt(8*MonomerLength*KuhnLength/(3*PI)) * p*Li(-0.5, 1-p)/(1-p);
	double ActualAvgEndToEnd = CalculateAverageEndToEndDistance();
	int ExpectedNumberOfMonomers = int(DomainVolume * MonomerNumberDensity + 0.5);
	int ActualNumberOfMonomers = GetTotalMonomerCount();
	double ExpectedMonomersPerChain = 1/p;
	double ActualMonomersPerChain = ActualNumberOfMonomers / (1.*NumChains);
	double ExpectedAvgTortuosity = sqrt(6*MonomerLength/(PI*KuhnLength)) * p*Li(-0.5, 1-p)/(1-p);
	double ActualAvgTortuosity = CalculateAverageTortuosity();

	cout << "Expected system density = " << ExpectedDensity << endl;
	cout << "Actual system density = " << ActualDensity << " = " << ActualDensity / ExpectedDensity << " x Expected" << endl;
	cout << "Expected average end-to-end length = " << ExpectedAvgEndToEnd << endl;
	cout << "Actual average end-to-end length = " << ActualAvgEndToEnd << " = " << ActualAvgEndToEnd / ExpectedAvgEndToEnd << " x Expected" << endl;
	cout << "Expected monomer count = " << ExpectedNumberOfMonomers << endl;
	cout << "Actual monomer count = " << ActualNumberOfMonomers << " = " << ActualNumberOfMonomers / (1.*ExpectedNumberOfMonomers) << " x Expected" << endl;
	cout << "Expected average number of monomers per chain = " << ExpectedMonomersPerChain << endl;
	cout << "Actual average number of monomers per chain = " << ActualMonomersPerChain << " = " << ActualMonomersPerChain / ExpectedMonomersPerChain << " x Expected" << endl;
	cout << "Expected average tortuosity = " << ExpectedAvgTortuosity << endl;
	cout << "Actual average tortuosity = " << ActualAvgTortuosity << " = " << ActualAvgTortuosity / ExpectedAvgTortuosity << " x Expected" << endl;
	cout << endl;

	cout.flush();

	return;
}





/*

	Given the end-to-end length of a random walk, generates a random number of steps.

	(Gaussian limit.)

*/

int Network::GetNumberOfSteps(double r){

/*
	static double scalefactor = 1.5 / (MonomerLength*KuhnLength);
	static double mu = fabs(log(1 - pxL));
	static double lim = 0.001;//HACK

	double a = r*r*scalefactor, w, wtot;
	int N, Npick;

	w = wtot = exp(-(a + mu));
	N = Npick = 1;
	while(true){
		++N;
		w = exp(-(a/N + mu*N)) / pow(N, 1.5);
		wtot += w;
		if(wtot <= 0) continue;
		if(w <= lim*wtot) break;
		if(rand1(Seed)*wtot <= w) Npick = N;
	}

	return Npick;
*/

	static double scalefactor = 1.0 / sqrt(MonomerLength*KuhnLength);
	static double q = 1 - pxL; 
	static double inva = -1./log(q);
	static double Em3_2 = exp(-1.5);
	static double E1    = exp(1.0);

	double x = r * scalefactor, x2 = x*x;
	double Nstar, temp;

	// Exponential, kissing h(N-1)
	double C1 = Em3_2 / (x*x2 * q);
	double A1 = C1 * inva;

	// Gamma, k = 2, kissing h(N)
	double C2 = pow(2.5/(x2*E1), 2.5);
	double A2 = C2 * inva * inva;

	// Gamma, k = 2, kissing h(N-1)
	temp = 3*x2 + 7;
	Nstar = (temp + sqrt(temp*temp - 40)) / 10;
	double C3 = exp(-1.5*x2 / (Nstar-1)) / (q * Nstar * pow(Nstar-1, 1.5));
	double A3 = C3 * inva * inva;

	double Trial, g, h;
	int Trialn;

	do{

		if(A1 <= A2 && A1 <= A3){
			// Sample the exponential distribution
			Trial = -log(1 - rand1(Seed)) * inva;
			g = C1*pow(q, Trial);
		}

		else{
			// Sample the gamma distribution with k = 2
			Trial = -log((1 - rand1(Seed))*(1 - rand1(Seed))) * inva;
			if(A2 <= A3)
				g = C2*Trial*pow(q, Trial);
			else
				g = C3*Trial*pow(q, Trial);
		}

		Trialn = int(Trial);
		h = (Trialn <= 0 ? 0 : pow(q, Trialn)*exp(-1.5*x2/Trialn) / pow(Trialn, 1.5));

	} while(g*rand1(Seed) > h);

	return Trialn;
}





/*

	Probability density function of the end-to-end vector of an n-step random walk.
	Here $r = |\vec{r}|$, and n is the number of steps.

*/

double Network::PDF_r_given_n(double r, int N){

	if(N <= 0) return 0;

	return pow(3/(2*PI*MonomerLength*KuhnLength*N), 1.5) * exp(-3*r*r/(2*MonomerLength*KuhnLength*N));
}





/*

	A close approximation to the pdf of the end to end vector of a random walk in the Gaussian limit.

*/

double Network::PDF_r(double r){

	// If first call, initialize several complex-to-compute variables critical
	// to the small-r behavior of this pdf.

	static int first = 0;
	static double sqrt3pi2 = sqrt(3*PI/2);
	static double sqrtpi = sqrt(PI);
	static double a, xfactor, sqrt6a, returnfactor;
	static double Li_3_2;	/* Polylogarithm function evaluated at 3/2, 1-pxL */
	static double Li_5_2;	/* Polylogarithm function evaluated at 5/2, 1-pxL */
	static double Li_7_2;	/* Polylogarithm function evaluated at 7/2, 1-pxL */

	if(first == 0){
		a = -log(1-pxL);
		xfactor = 1./sqrt(MonomerLength*KuhnLength);
		sqrt6a = sqrt(6.*a);
		returnfactor = 1. / (2*PI*PI * pow(MonomerLength*KuhnLength, 1.5));
		Li_3_2 = Li(1.5, 1-pxL);
		Li_5_2 = Li(2.5, 1-pxL);
		Li_7_2 = Li(3.5, 1-pxL);
		first = 1;
	}


	// Now, determine if r is in the "small regime" or the "large regime"

	if(r < 0)
		return 0;

	double PsiApprox, x = xfactor * r;

	if(x <= 0.6666666666667){
		// If x is small, use Pade approximant
		PsiApprox = (3*pxL/(1-pxL))*sqrt3pi2 * ( Li_3_2 - 6*pow(Li_5_2*x, 2)/(4*Li_5_2 + 3*Li_7_2*x*x) );
	}
	else{
		// Otherwise use an approximation to the full sum of residues
		PsiApprox = (3*PI*pxL / (x*(1-pxL))) * ( exp(-x*sqrt6a) - exp(-x*(sqrtpi*(1+x) + x*x*x)) );
	}

	return returnfactor * PsiApprox;
}





/*

	Calculates the radius R such that the probability
	of a chain of length >= R is f.
	This implicitly is based on the marginal r density given by PDF_r(double r).

*/

double Network::MaxRadius(double f){

	if(f < 0 || f >=1){
		cout << "Error in MaxRadius: f = " << f << " is out of bounds." << endl;
		exit(1);
	}

	double a = fabs(log(1-pxL));

	return -sqrt(MonomerLength*KuhnLength / (6*a)) * ( 1 + LambertWm1(-(1-pxL)*a*f / (exp(1.)*pxL)) );
}
