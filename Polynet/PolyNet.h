/*

	PolyNet: A mesoscale model of a polymer network.
	This is based on the original EPnet algorithm of David E. Hanson, with some
	modifications.

	Header file.

													John L. Barber
													jlbarber@lanl.gov
													505-664-0605

	4/25/2014
	Add cell structure for chain generation.
	Add entire input deck machinery, with default values for each parameter
	in PolyNet_Lib.h

*/

#ifndef POLYNET_HEADER
#define POLYNET_HEADER

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <float.h>
#include <vector>
#include <sstream>
#include <time.h>
#include <ctime>

#include "Point3D.h"
#include "ReadInputDeck.h"

using namespace std;


// Defined constants

#define KB						1.380650e-23		// Boltzmann's constant k_b (J/K)
#define PI						3.14159265358979
#define BACKGROUND_PIXEL_VALUE	-47110815 			// Value to use in an image array to indicate that the given pixel is not a part of the "object", and should be colored as background. 


// Defaults

const long Seed_Default = 0;
const bool LoadNetworkFlag_Default = false;
const string NetworkFileName_Default = "movie2.txt";

const double InitialDomainWidth_Default = 100e-9;
const bool PeriodicFlag_Default = true;

const double MonomerNumberDensity_Default = 8.13347e27;
const double NodeNumberDensity_Default = 4e25;
const double MassNode_Default = 0.0;
const double MassMonomer_Default = 1.13113e-25;
const double MonomerLength_Default = 4.332e-10;
const double KuhnLength_Default = 9.58e-10;
const int ChainsPerNode_Default = 4;
const double RadiusProbabilityCutoff_Default = 0.01;

const double Temperature_Default = 300.;
const double Keps_Default = 2.00e-13;
const double Kt_Default = 5.0e-14;
const double DragForce_Default = 0.01e-9;
const double BreakForce_Default = 6.8e-9;
const double Gamma_Default = 4e0;

const int DeformationType_Default = 0;
const double StrainRate_Default = +1e-4;
const double LambdaMax_Default = +30.0;
const double TimeMax_Default = 1e4;
const double DtCoefficient_Default = 0.01;
const bool FreezeEndsFlag_Default = false;
const double FreezeFraction_Default = 0.05;

const string OutputFolder_Default = "./";

const bool StrainStressFileFlag_Default = true;
const string StrainStressFile_Default = "StrainStress.txt";

const bool DrawNetworkImagesFlag_Default = true;
const int NetworkImageInterval_Default = 100;
const int NumberOfPixels_Default = 1000;
const string NetworkImageName_Default = "NetworkImage";
const string ImageFormat_Default = "png";

const double ViewPointX_Default = 0.35;
const double ViewPointY_Default = 0.65;
const double ViewPointZ_Default = 0.00;
const bool   RescaleViewPointFlag_Default = false;

const double ViewDirectionX_Default = 0; // Should only be relevant if ViewPoint is at origin
const double ViewDirectionY_Default = 0;
const double ViewDirectionZ_Default = 1;

const double UpX_Default = 1;
const double UpY_Default = 0;
const double UpZ_Default = 0;

const double AngularSpan_Default = 1.3;

const double NodeRadius_Default = 0.80;
const double ChainRadius_Default = 0.20;
const double BorderRadius_Default = 0;

const string BackgroundColor_Default = "White";
const string Node0Color_Default = "Black";
const string Node1Color_Default = "Black";
const string LooseChainColor_Default = "Blue";
const string TightChainColor_Default = "Red";
const string BorderColor_Default = "Green";


// Classes

class Node{
public:
	Point r, F;
	double m;
	double T;
	int NumChains;
	int tag;
};

class Chain{
public:
	int Node0, Node1;
	double r0;
	int N;
	int tag;
};

class Network{
public:

	string InputDeckName;

	vector<Node>  NodeList;
	vector<Chain> ChainList;
	int NumNodes;
	int NumChains;

	long Seed;
	bool LoadNetworkFlag;
	string NetworkFileName;
	double L0[3][2];
	bool PeriodicFlag;

	/* Material properties */
	double MonomerNumberDensity;
	double NodeNumberDensity;
	double MassNode;
	double MassMonomer;
	double MonomerLength;
	double KuhnLength;
	int ChainsPerNode;
	double MaximumRadius;
	double pxL;
	double r2overNbm;

	/* Current state */
	double Time;
	double L[3][2];
	double Lambda[3];
	double Stress[3][3];
	int    NumRupturedChains;

	/* Dynamics */

	double Temperature;
	double Keps;
	double Kt;
	double DragForce;
	double BreakForce;
	double Gamma;

	int DeformationType;
	double StrainRate;
	double LambdaMax;
	double TimeMax;
	double Dt;
	int NumberOfTimeSteps;
	bool FreezeEndsFlag;
	double FreezeFraction;

	// Output
	string OutputFolder;
	bool StrainStressFileFlag;
	string StrainStressFile;
	ofstream StrainStressFileStream;
	bool DrawNetworkImagesFlag;
	int NetworkImageInterval;
	int NumberOfPixels;
	string NetworkImageName;
	string ImageFormat;

	Point ViewPoint, ViewDirection, Up, LightDirection;
	Point xhat, yhat, zhat;
	bool RescaleViewPointFlag;
	double AngularSpanX, AngularSpanY, DTan;
	double NodeRadius, ChainRadius, BorderRadius;
	unsigned char BackgroundColor[3], Node0Color[3], Node1Color[3], LooseChainColor[3], TightChainColor[3], BorderColor[3];


	// Setup
	bool   LoadNetworkParameters();
	bool   LoadNetwork(string Format = "EPNet");
	int    MarkUnconnectedNodes();
	bool   CreateNodes();
	bool   CreateChains();
	void   PrintNetworkParameters();
	int    GetNumberOfSteps(double r);
	double PDF_r_given_n(double r, int N);
	double PDF_r(double r);
	double MaxRadius(double f);

	// Stepping
	bool   Step();
	int    Deform(double TimeIncrement);
	void   UpdateLambdas(double t);
	void   CalculateForces();
	double Force(double r, double r0, int NumberOfMonomers);

	// Diagnostics
	void   ComputeStressTensor();
	double CalculateAverageTortuosity();
	double CalculateAverageEndToEndDistance();
	double GetTotalNetworkMass();
	int    GetTotalMonomerCount();
	int    CountRupturedChains();
	double GetDomainDiagonal();

	// Output
	bool   WriteStrainStressData(bool FirstTimeFlag);
	bool   WriteNetworkImage(string ImageFileName);
	bool   RenderSphere(unsigned char* ImageArray, double* Depth, int nx, int ny, const Point& r0, double R, unsigned char* Color, bool PeriodicFlag);
	bool   RenderCylinder(unsigned char* ImageArray, double* Depth, int nx, int ny, const Point& r0, const Point& r1, double R, unsigned char* Color, bool PeriodicFlag);
	void   IlluminationModel(const Point& Normal, unsigned char* Color, unsigned char* FinalColor);
	bool   FindImageTangents(const Point& SurfacePoint, double& TanX, double& TanY);
};





// Miscellaneous functions

int    NumberOfDigits(int i);
void   CountdownWithSlash(int i, int j, bool DeleteFlag);
double Li(double n, double z);
double LambertWm1(double z);
double rand1(long& Seed);
double randn(long& Seed);
void Converter(string FileName, string OldExtension, string NewExtension);
void Periodify(Point& P, double *W);

#endif
