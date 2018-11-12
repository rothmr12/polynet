/*

	PolyNet: A mesoscale model of a polymer network.
	This is based on the original EPnet algorithm of David E. Hanson, with some
	modifications.

	Main file.

													John L. Barber
													jlbarber@lanl.gov
													505-664-0605

	4/25/2014
	Add cell structure for chain generation.
	Add entire input deck machinery, with default values for each parameter
	in PolyNet_Lib.h

*/


#include "PolyNet.h"

using namespace std;


// Main function

int main(int NumberOfCommandLineArguments, char** CommandLineArguments){

	// Check command line arguments
	if(NumberOfCommandLineArguments <= 1){
		cout << endl;
		cout << "Syntax: " << CommandLineArguments[0] << " [Input Deck Name]" << endl;
		cout << endl;
		return 0;
	}


	// Create Network object and load parameters
	clock_t c_begin = clock(), c_start;
	Network PolymerNetwork;
	PolymerNetwork.InputDeckName = CommandLineArguments[1];
	if(!PolymerNetwork.LoadNetworkParameters())
		exit(1);


	// Load in saved network, if called for...
	if(PolymerNetwork.LoadNetworkFlag){
		c_start = clock();
		cout << "Loading network from file " << PolymerNetwork.NetworkFileName << " ... ";
		cout.flush();

		if(!PolymerNetwork.LoadNetwork("EPNet"))
			exit(1);

		cout << ((double)(clock() - c_start)) / CLOCKS_PER_SEC << " seconds" << endl;
		cout.flush();
	}

	// ...otherwise, create one.
	else{
		// Generate node positions
		if(!PolymerNetwork.CreateNodes())
			exit(1);

		// Generate chains
		if(!PolymerNetwork.CreateChains())
			exit(1);
	}

	// Calculate the net force on each node in the initial state of the
	// newly-created network. (Should be 0 everywhere initially.)
	PolymerNetwork.CalculateForces();

	// Print network parameters, statistics
	cout << endl;
	PolymerNetwork.PrintNetworkParameters();

	// Initialize strain/stress file, if called for
	if(PolymerNetwork.StrainStressFileFlag && !PolymerNetwork.WriteStrainStressData(true))
		exit(1);

	// Draw image of initial network
	string IntString;
	stringstream ss;
	string ImageFileName;
	int StepNumber = 0;

	if(PolymerNetwork.DrawNetworkImagesFlag){
		ImageFileName = PolymerNetwork.NetworkImageName + "." + NumberString(StepNumber, PolymerNetwork.NumberOfTimeSteps);
		if(!PolymerNetwork.WriteNetworkImage(ImageFileName))
			exit(1);
	}

	// Stepping
	double OldDiagonal, NewDiagonal;
	bool MakeImageLastStep = false;

	for(StepNumber = 1; StepNumber <= PolymerNetwork.NumberOfTimeSteps; ++StepNumber){

		if(StepNumber == 1 || MakeImageLastStep){
			c_start = clock();
			cout << "Step: ";
			CountdownWithSlash(StepNumber, PolymerNetwork.NumberOfTimeSteps, false);
			cout.flush();
		}
		else
			CountdownWithSlash(StepNumber, PolymerNetwork.NumberOfTimeSteps, true);

		// One step of network deformation/evolution/whatever
		OldDiagonal = PolymerNetwork.GetDomainDiagonal();
		if(!PolymerNetwork.Step())
			exit(1);
		NewDiagonal = PolymerNetwork.GetDomainDiagonal();
		if(PolymerNetwork.RescaleViewPointFlag)
			PolymerNetwork.ViewPoint = (NewDiagonal / OldDiagonal) * PolymerNetwork.ViewPoint;

		// Initialize strain/stress file, if called for
		if(PolymerNetwork.StrainStressFileFlag && !PolymerNetwork.WriteStrainStressData(false))
			exit(1);


		// Make image of network at this timestep, if called for
		MakeImageLastStep = (PolymerNetwork.DrawNetworkImagesFlag && PolymerNetwork.NetworkImageInterval > 0 && (StepNumber % PolymerNetwork.NetworkImageInterval == 0));
		if(MakeImageLastStep){
			cout << "  " << ((double)(clock() - c_start)) / CLOCKS_PER_SEC << " sec" << endl;
			ImageFileName = PolymerNetwork.NetworkImageName + "." + NumberString(StepNumber, PolymerNetwork.NumberOfTimeSteps);
			if(!PolymerNetwork.WriteNetworkImage(ImageFileName))
				exit(1);
		}
	}
	cout << endl;


	// Timing

	cout << "Total time: " << ((double)(clock() - c_begin)) / CLOCKS_PER_SEC << " sec" << endl;

	return 0;
}
