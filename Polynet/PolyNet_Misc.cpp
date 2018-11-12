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

	Returns the number of characters (digits and minus sign, if present) in
	the integer i.

*/

int NumberOfDigits(int i){

	int NumberOfDigits = 0;

	if(i < 0 || i == -0){
		++NumberOfDigits;
		i = -i;
	}

	while(i != 0){
		i = i/10;
		++NumberOfDigits;
	}

	return NumberOfDigits;
}





/*

	For use in accomplishing running countdowns (or ups) on a single line.
	Prints [i] / [j]. If DeleteFlag = true, it first deletes a number of spaces
	equal to the [i] / [j] string it is about to print.

	Typical usage:

	for(i = 0; i < n; i++){
		CountdownWithSlash(i+1, n, !(i == 0));
		...
	}

*/

void CountdownWithSlash(int i, int j, bool DeleteFlag){

	int jDigits = NumberOfDigits(j);
	if(DeleteFlag){
		for(int k = 0; k < 2*jDigits + 3; ++k)
			cout << "\b";
	}
	cout << setw(jDigits) << i;
	cout << " / " << j;
	cout.flush();

	return;
}





/*

	Generates values of the polylogarithm function,

		Li_n(z) = \sum_{k=1}^{\infty} \frac{z^k}{k^n}

	for real z < 1 but also close to 1. This function is incomplete,
	and only works for certain special values of n which are integer
	or half-integer. The general evaluation is difficult.

*/

double Li(double n, double z){

	if(z <= 0 || z >= 1){
		cout << "Error in Li: z = " << z << " is out of bounds" << endl;
		flush(cout);
		exit(1);
	}

	if(n == 1)
		return -log(1-z);

	if(n == 0)
		return z / (1-z);

	if(n == -1)
		return z / pow(1-z, 2);

	if(n == -2)
		return z*(1+z) / pow(1-z, 3);

	if(n == -3)
		return z*(1+z*(4+z)) / pow(1-z, 4);

	if(n == -4)
		return z*(1+z*(11+z*(11+z))) / pow(1-z, 5);

	double zcut;
	if     (n == -0.5) zcut = 0.45;
	else if(n ==  0.5) zcut = 0.265;
	else if(n ==  1.5) zcut = 0.370;
	else if(n ==  2.5) zcut = 0.555;
	else if(n ==  3.5) zcut = 0.743;
	else{
		cout << "Error in Li: n = " << n << " is not yet defined" << endl;
		flush(cout);
		exit(1);
	}

	if(z <= zcut)
		return 0.5*z*( 2 + z*(pow(2,-n) + z/pow(3,n) + 1/(pow(2,n)-z)) );

	static double Zeta_m3_2 = -0.025485202; /* Riemann zeta fn at -3/2 */
	static double Zeta_m1_2 = -0.207886225; /* Riemann zeta fn at -1/2 */
	static double Zeta_1_2  = -1.460354509; /* Riemann zeta fn at +1/2 */
	static double Zeta_3_2  =  2.612375349; /* Riemann zeta fn at +3/2 */
	static double Zeta_5_2  =  1.341487257; /* Riemann zeta fn at +5/2 */
	static double Zeta_7_2  =  1.126733867; /* Riemann zeta fn at +7/2 */

	if(n == -0.5){
		double p = 1-z;
		return sqrt(PI)*(32-p*(24+p))/(64*pow(p,1.5)) + Zeta_m1_2;
	}

	double a = -log(z);

	if(n == +0.5)
		return sqrt(PI/a) + Zeta_1_2 - a*(Zeta_m1_2 - Zeta_m3_2*a/2);

	if(n == +1.5)
		return -2*sqrt(PI*a) + Zeta_3_2 - a*(Zeta_1_2 - Zeta_m1_2*a/2);

	if(n == +2.5)
		return Zeta_5_2 + (a/6)*(8*sqrt(PI*a) + 3*Zeta_1_2*a - 6*Zeta_3_2);

	if(n == +3.5)
		return Zeta_7_2 + a*(-Zeta_5_2 + a*(Zeta_3_2/2 - 8*sqrt(PI*a)/15));

	return 0;
}





/*

	Returns the -1 branch of the Lambert W(z) function.

	This function only works for -1/e <= z < 0, and is based on the
	scheme outlined in figure 2 of
	IEEE TRANSACTIONS ON SIGNAL PROCESSING, VOL. 50, NO. 9, pp. 2160-64, SEPTEMBER 2002

*/

double LambertWm1(double z){

	double W = 0;

	if(-1/exp(1.) <= z && z < -0.333){
		double p = -sqrt(2*(exp(1.)*z + 1));
		W = -1 + p*(1 + p*(-1./3. + p*(11./72. + p*(-43./540. + p*(769./17280. - 221*p/8505.)))));
	}

	else if(-0.333 <= z && z <= -0.033){
		double Numerator = -8.0960 + z*(391.0025 + z*(-47.4252 + z*(-4877.6330 - 5532.7760*z)));
		double Denominator = 1 + z*(-82.9423 + z*(433.8688 + 1515.3060*z));
		W = Numerator / Denominator;
	}

	else if(-0.033 < z && z < 0){
		double L1 = log(-z);
		double L2 = log(-log(-z));
		W  = L1 - L2 + L2/L1;
		W += (-2+L2)*L2 / (2*L1*L1);
		W += (6 + L2*(-9 + 2*L2))*L2 / (6*pow(L1,3));
		W += (-12 + L2*(36 + L2*(-22 + 3*L2)))*L2 / (12*pow(L1,4));
		W += (60 + L2*(-300 + L2*(350 + L2*(-125 + 12*L2))))*L2 / (60*pow(L1,5));
	}

	return W;
}





/*

	Generates a uniform random deviate in [0,1)

*/

double rand1(long& Seed){

	static double a = 16807.0;
	static double m = 2147483647.0;

	Seed = (long)( fmod(a*Seed, m) );

	return Seed/m;
}





/*

	Generates a standard normal deviate, mean 0, standard deviation 1

*/

double randn(long& Seed){

	static int    rand_avail = 0;
	static double next_rand;

	double c, u1, u2;

	if(rand_avail == 0){
		do{
			u1 = 2.*rand1(Seed) - 1.;
			u2 = 2.*rand1(Seed) - 1.;
			c = u1*u1 + u2*u2;
		} while(c >= 1. || c == 0.);
		c = sqrt( -2.*log(c)/c );
		next_rand  = u1*c;
		rand_avail = 1;
		return u2*c;
	}
	else{
		rand_avail = 0;
		return next_rand;
	}
}





/*

	Uses whatever image conversion command is defined in the DRACO_Misc.h header as ImageConvertCommand
	to convert an original file [FileName].[OldExtension] to [FileName].[NewExtension], and then
	deletes the original file. (Deletion is done using FileRemoveCommand, also defined in the header.)

	Note that the "." should not be included in the file extension strings.

*/

const string ImageConvertCommand = "convert";	// Image conversion command available on this terminal, for converting ppm to gif. (ImageMagick = "convert")
const string FileRemoveCommand = "rm";			// Terminal command for deleting files. (Unix: "rm")

void Converter(string FileName, string OldExtension, string NewExtension){

	string OldFileName = FileName + "." + OldExtension;
	string NewFileName = FileName + "." + NewExtension;
	string SystemString = ImageConvertCommand + " " + OldFileName + " " + NewFileName;
	system(SystemString.c_str());
	SystemString = FileRemoveCommand + " " + OldFileName;
	system(SystemString.c_str());

	return;
}



/*

	Alters the Point P by assuming periodicity and restricting P.r[i] to [-W[i]/2, +W[i]/2] for
	i = 0, 1, 2

*/

void Periodify(Point& P, double *W){

	if     (P.r[0] > +0.5*W[0]) P.r[0] -= W[0];
	else if(P.r[0] < -0.5*W[0]) P.r[0] += W[0];
	if     (P.r[1] > +0.5*W[1]) P.r[1] -= W[1];
	else if(P.r[1] < -0.5*W[1]) P.r[1] += W[1];
	if     (P.r[2] > +0.5*W[2]) P.r[2] -= W[2];
	else if(P.r[2] < -0.5*W[2]) P.r[2] += W[2];

	return;
}
