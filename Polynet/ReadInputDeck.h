/*

	Header file for functions associated with reading parameter values (strings, ints, doubles,
	bools, or unsigned chars) from an input deck file.

	All of these functions assume a text-based input deck syntax with lines like:

	[ParameterName] <Any amount of whitespace> [One or more values separated by whitespace] <Any amount of whitespace> [# optional comment]

	Here [ParameterName] is a name or label for a given parameter. "#" is the comment character,
	so everything after a "#" on a line is ignored. The order and spacing of the lines in an
	input deck is irrelevant, these functions will search through the entire deck for a given
	parameter label. The first instance (from the top) in which a given parameter label is
	encountered is the one which will be used. Parameter labels are case-insensitive. Boolean
	values are represented in the input deck as numbers: 0 = false, not-0 = true.

	GetColor is a special function for reading in the value(s) of a parameter that represents
	a color, either as 3 RGB values 0-255 or as one of a number of specially-defined strings
	("Red", "Green", "Pink", etc) which are then mapped to RGB. 

	John L. Barber
	jlbarber@lanl.gov
	
*/

#ifndef READINPUTDECK_HEADER
#define READINPUTDECK_HEADER

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;	// Yeah, yeah, I know. But who has time to type "std::" over and over?

bool GetParameters (ifstream& ParameterStream, string ParameterName, string *ParameterValues, int n, string DefaultValue);
bool GetParameters2(ifstream& ParameterStream, string ParameterName, string *ParameterValues, int n, string DefaultValue);
bool GetParameters (ifstream& ParameterStream, string ParameterName, string *ParameterValues, int n);
bool GetParameter  (ifstream& ParameterStream, string ParameterName, string &ParameterValue, string DefaultValue);
bool GetParameter2 (ifstream& ParameterStream, string ParameterName, string &ParameterValue, string DefaultValue);
bool GetParameter  (ifstream& ParameterStream, string ParameterName, string &ParameterValue);

bool GetParameters(ifstream& ParameterStream, string ParameterName, int *ParameterValues, int n, int DefaultValue);
bool GetParameters2(ifstream& ParameterStream, string ParameterName, int *ParameterValues, int n, int DefaultValue);
bool GetParameters(ifstream& ParameterStream, string ParameterName, int *ParameterValues, int n);
bool GetParameter (ifstream& ParameterStream, string ParameterName, int &ParameterValue, int DefaultValue);
bool GetParameter2(ifstream& ParameterStream, string ParameterName, int &ParameterValue, int DefaultValue);
bool GetParameter (ifstream& ParameterStream, string ParameterName, int &ParameterValue);

bool GetParameters (ifstream& ParameterStream, string ParameterName, double *ParameterValues, int n, double DefaultValue);
bool GetParameters2(ifstream& ParameterStream, string ParameterName, double *ParameterValues, int n, double DefaultValue);
bool GetParameters(ifstream& ParameterStream, string ParameterName, double *ParameterValues, int n);
bool GetParameter (ifstream& ParameterStream, string ParameterName, double &ParameterValue, double DefaultValue);
bool GetParameter2(ifstream& ParameterStream, string ParameterName, double &ParameterValue, double DefaultValue);
bool GetParameter (ifstream& ParameterStream, string ParameterName, double &ParameterValue);

bool GetParameters (ifstream& ParameterStream, string ParameterName, bool *ParameterValues, int n, bool DefaultValue);
bool GetParameters2(ifstream& ParameterStream, string ParameterName, bool *ParameterValues, int n, bool DefaultValue);
bool GetParameters(ifstream& ParameterStream, string ParameterName, bool *ParameterValues, int n);
bool GetParameter (ifstream& ParameterStream, string ParameterName, bool &ParameterValue, bool DefaultValue);
bool GetParameter2(ifstream& ParameterStream, string ParameterName, bool &ParameterValue, bool DefaultValue);
bool GetParameter (ifstream& ParameterStream, string ParameterName, bool &ParameterValue);

bool GetParameters(ifstream& ParameterStream, string ParameterName, unsigned char *ParameterValues, int n, unsigned char DefaultValue);
bool GetParameters2(ifstream& ParameterStream, string ParameterName, unsigned char *ParameterValues, int n, unsigned char DefaultValue);
bool GetParameters(ifstream& ParameterStream, string ParameterName, unsigned char *ParameterValues, int n);
bool GetParameter (ifstream& ParameterStream, string ParameterName, unsigned char &ParameterValue, unsigned char DefaultValue);
bool GetParameter2(ifstream& ParameterStream, string ParameterName, unsigned char &ParameterValue, unsigned char DefaultValue);
bool GetParameter (ifstream& ParameterStream, string ParameterName, unsigned char &ParameterValue);

bool GetColor(ifstream& ParameterStream, string ParameterName, unsigned char *ParameterValues, string DefaultColor = "NoDefault");
string ColorToString(unsigned char *RGB);

bool CapitalCompare(const string& s1, const string& s2);
string NumberString(int n, int m = 0);
string NumberString(double x);
int HowManyDigits(int n);

#endif
