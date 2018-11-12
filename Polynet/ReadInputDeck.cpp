/*

	Main file for functions associated with reading parameter values (strings, ints, doubles,
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

#include "ReadInputDeck.h"





/*

	For reading a sequence of n strings from the input deck file connected to the stream
	ParameterStream. The n strings found will be placed in the array ParameterValues. If
	there are < n values given in the input deck file, or if the label ParameterName is
	missing from the input deck altogether, then the remaining entries in the array
	ParameterValues will be filled by DefaultValue.

*/

bool GetParameters(ifstream& ParameterStream, string ParameterName, string *ParameterValues, int n, string DefaultValue){

	if(n <= 0)
		return true;

	char ParameterLine[1024];
	int  i;
	string Substring;
	istringstream iss;

	ParameterStream.clear();
	ParameterStream.seekg(0, ios::beg);

	// Cycle through each line in the file
	while(ParameterStream.getline(ParameterLine, 1024)){
		i = 0;	// Keeps track of which substring of this line we're on
		iss.str(ParameterLine);

		// Check each substring of this line
		while(iss >> Substring){

			// Check if the rest of this line is commented out
			if(Substring[0] == '#')
				break;

			// If the first substring is anything other than our
			// ParameterName, then this isn't the right line
			if(i == 0 && !CapitalCompare(Substring, ParameterName))
				break;

			// If we've gotten past the first substring, then this substring
			// is one of the values we're looking for 
			if(i > 0)
				ParameterValues[i-1] = Substring;

			// If, at this point, we've gotten past one label plus n values,
			// so we're done
			if(++i >= n+1)
				return true;
		}

		// If i >= 2, then we've found the label ParameterName followed
		// by at least one value, so fill the remainder of the ParameterValues array
		// with the default. If we've found the label followed by no values, we don't
		// count it, and keep looking
		if(i >= 2){
			for(--i; i < n; ++i)
				ParameterValues[i] = DefaultValue;
			return true;
		}

		iss.clear();
	}

	// If we get to this point, we've gone through every line in the file without
	// finding the label ParameterName followed by one or more value, so fill the
	// whole array with the default.
	for(i = 0; i < n; ++i)
		ParameterValues[i] = DefaultValue;

	return true;
}





/*

	For reading a sequence of n strings from the input deck file connected to the stream
	ParameterStream. The n strings found will be placed in the array ParameterValues. If
	there are < n values given in the input deck file, then the remaining entries in the array
	ParameterValues will be filled by DefaultValue, and true will be returned.

	If the label ParameterName is missing from the input deck altogether, false will be
	returned.

*/

bool GetParameters2(ifstream& ParameterStream, string ParameterName, string *ParameterValues, int n, string DefaultValue){

	if(n <= 0)
		return true;

	char ParameterLine[1024];
	int  i;
	string Substring;
	istringstream iss;

	ParameterStream.clear();
	ParameterStream.seekg(0, ios::beg);

	// Cycle through each line in the file
	while(ParameterStream.getline(ParameterLine, 1024)){
		i = 0;	// Keeps track of which substring of this line we're on
		iss.str(ParameterLine);

		// Check each substring of this line
		while(iss >> Substring){

			// Check if the rest of this line is commented out
			if(Substring[0] == '#')
				break;

			// If the first substring is anything other than our
			// ParameterName, then this isn't the right line
			if(i == 0 && !CapitalCompare(Substring, ParameterName))
				break;

			// If we've gotten past the first substring, then this substring
			// is one of the values we're looking for 
			if(i > 0)
				ParameterValues[i-1] = Substring;

			// If, at this point, we've gotten past one label plus n values,
			// so we're done
			if(++i >= n+1)
				return true;
		}

		// If i >= 2, then we've found the label ParameterName followed
		// by at least one value, so fill the remainder of the ParameterValues array
		// with the default. If we've found the label followed by no values, we don't
		// count it, and keep looking
		if(i >= 2){
			for(--i; i < n; ++i)
				ParameterValues[i] = DefaultValue;
			return true;
		}

		iss.clear();
	}

	return false;
}





/*

	For reading a sequence of n strings from the input deck file connected to the stream
	ParameterStream. The n strings found will be placed in the array ParameterValues. If
	there are < n values given in the input deck file, or if the label ParameterName is
	missing from the input deck altogether, then false is returned.

*/

bool GetParameters(ifstream& ParameterStream, string ParameterName, string *ParameterValues, int n){

	if(n <= 0)
		return true;

	char ParameterLine[1024];
	int  i;
	string Substring;
	istringstream iss;

	ParameterStream.clear();
	ParameterStream.seekg(0, ios::beg);

	// Cycle through each line in the file
	while(ParameterStream.getline(ParameterLine, 1024)){
		i = 0;	// Keeps track of which substring of this line we're on
		iss.str(ParameterLine);

		// Check each substring of this line
		while(iss >> Substring){

			// Check if the rest of this line is commented out
			if(Substring[0] == '#')
				break;

			// If the first substring is anything other than our
			// ParameterName, then this isn't the right line
			if(i == 0 && !CapitalCompare(Substring, ParameterName))
				break;

			// If we've gotten past the first substring, then this substring
			// is one of the values we're looking for 
			if(i > 0)
				ParameterValues[i-1] = Substring;

			// If, at this point, we've gotten past one label plus n values,
			// so we're done
			if(++i >= n+1)
				return true;
		}

		// If i >= 2, then we've found the label ParameterName followed
		// by at least one value. However, if i < n+1, then it wasn't followed
		// by enough values, and we return false. If we've found the label
		// followed by no values (i <= 1), we don't count it, and keep looking.
		if(i >= 2 && i < n+1)
			return false;

		iss.clear();
	}

	// If we've reached this point, we've searched every line and failed to find
	// the label ParameterName followed by a suffient number of values
	return false;
}





/*

	For reading a single string into ParameterValue from the input deck file connected
	to the stream ParameterStream. If there is no value given in the input deck file, or
	if the label ParameterName is missing altogether, then ParameterValue will be filled
	by DefaultValue.

*/

bool GetParameter (ifstream& ParameterStream, string ParameterName, string &ParameterValue, string DefaultValue){

	return GetParameters(ParameterStream, ParameterName, &ParameterValue, 1, DefaultValue);
}





/*

	For reading a single string into ParameterValue from the input deck file connected
	to the stream ParameterStream. If there is no value given in the input deck file then
	ParameterValue will be filled by DefaultValue, and true will be returned.

	If the label ParameterName is missing altogether, then false will be returned.

*/

bool GetParameter2(ifstream& ParameterStream, string ParameterName, string &ParameterValue, string DefaultValue){

	return GetParameters2(ParameterStream, ParameterName, &ParameterValue, 1, DefaultValue);
}





/*

	For reading a single string into ParameterValue from the input deck file connected
	to the stream ParameterStream. If there is no value given in the input deck file, or
	if the label ParameterName is missing altogether, then false is returned.

*/

bool GetParameter (ifstream& ParameterStream, string ParameterName, string &ParameterValue){

	return GetParameters(ParameterStream, ParameterName, &ParameterValue, 1);
}





/*

	For reading a sequence of n ints from the input deck file connected to the stream
	ParameterStream. The n ints found will be placed in the array ParameterValues. If
	there are < n values given in the input deck file, or if the label ParameterName is
	missing from the input deck altogether, then the remaining entries in the array
	ParameterValues will be filled by DefaultValue.

*/

bool GetParameters(ifstream& ParameterStream, string ParameterName, int *ParameterValues, int n, int DefaultValue){

	// Make an auxiliary array of strings so we can use the string version of GetParameters
	// to do the heavy lifting, and read in the values (initially) as strings
	string *Params = new string[n];
	GetParameters(ParameterStream, ParameterName, Params, n, NumberString(DefaultValue, 0));

	// Now we need to try and convert each of the loaded-in strings into ints
	stringstream ss;
	for(int i = 0; i < n; ++i){
		ss.clear();
		ss.str(Params[i]);

		// If >> fails here, one of the loaded-in strings cannot be
		// translated into an int. Use the default instead.
		if(!(ss >> ParameterValues[i]))
			ParameterValues[i] = DefaultValue;
	}

	delete [] Params;

	return true;
}





/*

	For reading a sequence of n ints from the input deck file connected to the stream
	ParameterStream. The n ints found will be placed in the array ParameterValues. If
	there are < n values given in the input deck file, then the remaining entries in the array
	ParameterValues will be filled by DefaultValue, and true will be returned.

	If the label ParameterName is missing from the input deck altogether, then false will be returned.

*/

bool GetParameters2(ifstream& ParameterStream, string ParameterName, int *ParameterValues, int n, int DefaultValue){

	// Make an auxiliary array of strings so we can use the string version of GetParameters
	// to do the heavy lifting, and read in the values (initially) as strings
	string *Params = new string[n];
	if(!GetParameters2(ParameterStream, ParameterName, Params, n, NumberString(DefaultValue, 0)))
		return false;

	// Now we need to try and convert each of the loaded-in strings into ints
	stringstream ss;
	for(int i = 0; i < n; ++i){
		ss.clear();
		ss.str(Params[i]);

		// If >> fails here, one of the loaded-in strings cannot be
		// translated into an int. Use the default instead.
		if(!(ss >> ParameterValues[i]))
			ParameterValues[i] = DefaultValue;
	}

	delete [] Params;

	return true;
}





/*

	For reading a sequence of n ints from the input deck file connected to the stream
	ParameterStream. The n ints found will be placed in the array ParameterValues. If
	there are < n values given in the input deck file, or if the label ParameterName is
	missing from the input deck altogether, then false will be returned.

*/

bool GetParameters(ifstream& ParameterStream, string ParameterName, int *ParameterValues, int n){

	// Make an auxiliary array of strings so we can use the string version of GetParameters
	// to do the heavy lifting, and read in the values (initially) as strings
	string *Params = new string[n];
	if(!GetParameters(ParameterStream, ParameterName, Params, n))
		return false;

	// Now we need to try and convert each of the loaded-in strings into ints
	stringstream ss;
	for(int i = 0; i < n; ++i){
		ss.clear();
		ss.str(Params[i]);

		// If >> fails here, one of the loaded-in strings cannot be
		// translated into an int
		if(!(ss >> ParameterValues[i]))
			return false;
	}

	delete [] Params;

	return true;
}





/*

	For reading a single int into ParameterValue from the input deck file connected
	to the stream ParameterStream. If there is no value given in the input deck file, or
	if the label ParameterName is missing altogether, then ParameterValue will be filled
	by DefaultValue.

*/

bool GetParameter (ifstream& ParameterStream, string ParameterName, int &ParameterValue, int DefaultValue){

	return GetParameters(ParameterStream, ParameterName, &ParameterValue, 1, DefaultValue);
}





/*

	For reading a single int into ParameterValue from the input deck file connected
	to the stream ParameterStream. If there is no value given in the input deck file
	then ParameterValue will be filled by DefaultValue, and true will be returned.

	If the label ParameterName is missing altogether, then false will be returned.

*/

bool GetParameter2(ifstream& ParameterStream, string ParameterName, int &ParameterValue, int DefaultValue){

	return GetParameters2(ParameterStream, ParameterName, &ParameterValue, 1, DefaultValue);
}





/*

	For reading a single int into ParameterValue from the input deck file connected
	to the stream ParameterStream. If there is no value given in the input deck file, or
	if the label ParameterName is missing altogether, then false will be returned.

*/

bool GetParameter (ifstream& ParameterStream, string ParameterName, int &ParameterValue){

	return GetParameters(ParameterStream, ParameterName, &ParameterValue, 1);
}





/*

	For reading a sequence of n doubles from the input deck file connected to the stream
	ParameterStream. The n doubles found will be placed in the array ParameterValues. If
	there are < n values given in the input deck file, or if the label ParameterName is
	missing from the input deck altogether, then the remaining entries in the array
	ParameterValues will be filled by DefaultValue.

*/

bool GetParameters(ifstream& ParameterStream, string ParameterName, double *ParameterValues, int n, double DefaultValue){

	// Make an auxiliary array of strings so we can use the string version of GetParameters
	// to do the heavy lifting, and read in the values (initially) as strings
	string *Params = new string[n];
	GetParameters(ParameterStream, ParameterName, Params, n, NumberString(DefaultValue));

	// Now we need to try and convert each of the loaded-in strings into doubles
	stringstream ss;
	for(int i = 0; i < n; ++i){
		ss.clear();
		ss.str(Params[i]);

		// If >> fails here, one of the loaded-in strings cannot be
		// translated into a double. Use the default instead.
		if(!(ss >> ParameterValues[i]))
			ParameterValues[i] = DefaultValue;
	}

	delete [] Params;

	return true;
}





/*

	For reading a sequence of n doubles from the input deck file connected to the stream
	ParameterStream. The n doubles found will be placed in the array ParameterValues. If
	there are < n values given in the input deck file then the remaining entries in the array
	ParameterValues will be filled by DefaultValue and true will be returned.

	If the label ParameterName is missing from the input deck altogether, then false will be returned.

*/

bool GetParameters2(ifstream& ParameterStream, string ParameterName, double *ParameterValues, int n, double DefaultValue){

	// Make an auxiliary array of strings so we can use the string version of GetParameters
	// to do the heavy lifting, and read in the values (initially) as strings
	string *Params = new string[n];
	if(!GetParameters2(ParameterStream, ParameterName, Params, n, NumberString(DefaultValue)))
		return false;

	// Now we need to try and convert each of the loaded-in strings into doubles
	stringstream ss;
	for(int i = 0; i < n; ++i){
		ss.clear();
		ss.str(Params[i]);

		// If >> fails here, one of the loaded-in strings cannot be
		// translated into a double. Use the default instead.
		if(!(ss >> ParameterValues[i]))
			ParameterValues[i] = DefaultValue;
	}

	delete [] Params;

	return true;
}





/*

	For reading a sequence of n doubles from the input deck file connected to the stream
	ParameterStream. The n doubles found will be placed in the array ParameterValues. If
	there are < n values given in the input deck file, or if the label ParameterName is
	missing from the input deck altogether, then false will be returned.

*/

bool GetParameters(ifstream& ParameterStream, string ParameterName, double *ParameterValues, int n){

	// Make an auxiliary array of strings so we can use the string version of GetParameters
	// to do the heavy lifting, and read in the values (initially) as strings
	string *Params = new string[n];
	if(!GetParameters(ParameterStream, ParameterName, Params, n))
		return false;

	// Now we need to try and convert each of the loaded-in strings into doubles
	stringstream ss;
	for(int i = 0; i < n; ++i){
		ss.clear();
		ss.str(Params[i]);

		// If >> fails here, one of the loaded-in strings cannot be
		// translated into a double
		if(!(ss >> ParameterValues[i]))
			return false;
	}

	delete [] Params;

	return true;
}





/*

	For reading a single double into ParameterValue from the input deck file connected
	to the stream ParameterStream. If there is no value given in the input deck file, or
	if the label ParameterName is missing altogether, then ParameterValue will be filled
	by DefaultValue.

*/

bool GetParameter (ifstream& ParameterStream, string ParameterName, double &ParameterValue, double DefaultValue){

	return GetParameters(ParameterStream, ParameterName, &ParameterValue, 1, DefaultValue);
}





/*

	For reading a single double into ParameterValue from the input deck file connected
	to the stream ParameterStream. If there is no value given in the input deck file
	then ParameterValue will be filled by DefaultValue, and true will be returned.

	If the label ParameterName is missing altogether, then false will be returned.

*/

bool GetParameter2(ifstream& ParameterStream, string ParameterName, double &ParameterValue, double DefaultValue){

	return GetParameters2(ParameterStream, ParameterName, &ParameterValue, 1, DefaultValue);
}





/*

	For reading a single double into ParameterValue from the input deck file connected
	to the stream ParameterStream. If there is no value given in the input deck file, or
	if the label ParameterName is missing altogether, then false will be returned.

*/

bool GetParameter (ifstream& ParameterStream, string ParameterName, double &ParameterValue){

	return GetParameters(ParameterStream, ParameterName, &ParameterValue, 1);
}





/*

	For reading a sequence of n bools from the input deck file connected to the stream
	ParameterStream. The n bools found will be placed in the array ParameterValues. If
	there are < n values given in the input deck file, or if the label ParameterName is
	missing from the input deck altogether, then the remaining entries in the array
	ParameterValues will be filled by DefaultValue.

*/

bool GetParameters(ifstream& ParameterStream, string ParameterName, bool *ParameterValues, int n, bool DefaultValue){

	// Make an auxiliary array of strings so we can use the string version of GetParameters
	// to do the heavy lifting, and read in the values (initially) as strings
	string *Params = new string[n];
	GetParameters(ParameterStream, ParameterName, Params, n, (DefaultValue ? "1" : "0"));

	// Now we need to try and convert each of the loaded-in strings into bools
	stringstream ss;
	double Temp;
	for(int i = 0; i < n; ++i){
		ss.clear();
		ss.str(Params[i]);

		// If >> fails here, one of the loaded-in strings cannot be
		// translated into a double, which would then be made into a bool.
		// Use the default instead.
		if(!(ss >> Temp))
			ParameterValues[i] = DefaultValue;
		else
			ParameterValues[i] = (Temp != 0);
	}

	delete [] Params;

	return true;
}





/*

	For reading a sequence of n bools from the input deck file connected to the stream
	ParameterStream. The n bools found will be placed in the array ParameterValues. If
	there are < n values given in the input deck file, then the remaining entries in the array
	ParameterValues will be filled by DefaultValue and true will be returned.

	If the label ParameterName is missing from the input deck altogether, then false will be returned.

*/

bool GetParameters2(ifstream& ParameterStream, string ParameterName, bool *ParameterValues, int n, bool DefaultValue){

	// Make an auxiliary array of strings so we can use the string version of GetParameters
	// to do the heavy lifting, and read in the values (initially) as strings
	string *Params = new string[n];
	if(!GetParameters2(ParameterStream, ParameterName, Params, n, (DefaultValue ? "1" : "0")))
		return false;

	// Now we need to try and convert each of the loaded-in strings into bools
	stringstream ss;
	double Temp;
	for(int i = 0; i < n; ++i){
		ss.clear();
		ss.str(Params[i]);

		// If >> fails here, one of the loaded-in strings cannot be
		// translated into a double, which would then be made into a bool.
		// Use the default instead.
		if(!(ss >> Temp))
			ParameterValues[i] = DefaultValue;
		else
			ParameterValues[i] = (Temp != 0);
	}

	delete [] Params;

	return true;
}





/*

	For reading a sequence of n bools from the input deck file connected to the stream
	ParameterStream. The n bools found will be placed in the array ParameterValues. If
	there are < n values given in the input deck file, or if the label ParameterName is
	missing from the input deck altogether, then false will be returned.

*/

bool GetParameters(ifstream& ParameterStream, string ParameterName, bool *ParameterValues, int n){

	// Make an auxiliary array of strings so we can use the string version of GetParameters
	// to do the heavy lifting, and read in the values (initially) as strings
	string *Params = new string[n];
	if(!GetParameters(ParameterStream, ParameterName, Params, n))
		return false;

	// Now we need to try and convert each of the loaded-in strings into bools
	stringstream ss;
	double Temp;
	for(int i = 0; i < n; ++i){
		ss.clear();
		ss.str(Params[i]);

		// If >> fails here, one of the loaded-in strings cannot be
		// translated into a double (which will then be made into a bool)
		if(!(ss >> Temp))
			return false;
		ParameterValues[i] = (Temp != 0);
	}

	delete [] Params;

	return true;
}





/*

	For reading a single bool into ParameterValue from the input deck file connected
	to the stream ParameterStream. If there is no value given in the input deck file, or
	if the label ParameterName is missing altogether, then ParameterValue will be filled
	by DefaultValue.

*/

bool GetParameter (ifstream& ParameterStream, string ParameterName, bool &ParameterValue, bool DefaultValue){

	return GetParameters(ParameterStream, ParameterName, &ParameterValue, 1, DefaultValue);
}





/*

	For reading a single bool into ParameterValue from the input deck file connected
	to the stream ParameterStream. If there is no value given in the input deck file
	then ParameterValue will be filled by DefaultValue and true will be returned.

	If the label ParameterName is missing altogether then false will be returned.

*/

bool GetParameter2(ifstream& ParameterStream, string ParameterName, bool &ParameterValue, bool DefaultValue){

	return GetParameters2(ParameterStream, ParameterName, &ParameterValue, 1, DefaultValue);
}





/*

	For reading a single bool into ParameterValue from the input deck file connected
	to the stream ParameterStream. If there is no value given in the input deck file, or
	if the label ParameterName is missing altogether, then false will be returned.

*/

bool GetParameter (ifstream& ParameterStream, string ParameterName, bool &ParameterValue){

	return GetParameters(ParameterStream, ParameterName, &ParameterValue, 1);
}


/*

	For reading a sequence of n unsigned chars from the input deck file connected to the stream
	ParameterStream. The n unsigned chars found will be placed in the array ParameterValues. If
	there are < n values given in the input deck file, or if the label ParameterName is
	missing from the input deck altogether, then the remaining entries in the array
	ParameterValues will be filled by DefaultValue.

*/

bool GetParameters(ifstream& ParameterStream, string ParameterName, unsigned char *ParameterValues, int n, unsigned char DefaultValue){

	// Make an auxiliary array of ints so we can use the int version of GetParameters
	// to do the heavy lifting, and read in the values (initially) as ints
	int *Params = new int[n];
	GetParameters(ParameterStream, ParameterName, Params, n, int(DefaultValue));

	// Convert each of the loaded-in ints into unsigned chars
	for(int i = 0; i < n; ++i){
		if(Params[i] <= 0)
			ParameterValues[i] = 0;
		else if(Params[i] >= 255)
			ParameterValues[i] = 255;
		else
			ParameterValues[i] = (unsigned char)(Params[i]);
	}

	delete [] Params;

	return true;
}





/*

	For reading a sequence of n unsigned chars from the input deck file connected to the stream
	ParameterStream. The n unsigned chars found will be placed in the array ParameterValues. If
	there are < n values given in the input deck file, then the remaining entries in the array
	ParameterValues will be filled by DefaultValue and true will be returned.

	If the label ParameterName is missing from the input deck altogether then false will be returned.

*/

bool GetParameters2(ifstream& ParameterStream, string ParameterName, unsigned char *ParameterValues, int n, unsigned char DefaultValue){

	// Make an auxiliary array of ints so we can use the int version of GetParameters
	// to do the heavy lifting, and read in the values (initially) as ints
	int *Params = new int[n];
	if(!GetParameters2(ParameterStream, ParameterName, Params, n, int(DefaultValue)))
		return false;

	// Convert each of the loaded-in ints into unsigned chars
	for(int i = 0; i < n; ++i){
		if(Params[i] <= 0)
			ParameterValues[i] = 0;
		else if(Params[i] >= 255)
			ParameterValues[i] = 255;
		else
			ParameterValues[i] = (unsigned char)(Params[i]);
	}

	delete [] Params;

	return true;
}





/*

	For reading a sequence of n unsigned chars from the input deck file connected to the stream
	ParameterStream. The n unsigned chars found will be placed in the array ParameterValues. If
	there are < n values given in the input deck file, or if the label ParameterName is
	missing from the input deck altogether, then false will be returned.

*/

bool GetParameters(ifstream& ParameterStream, string ParameterName, unsigned char *ParameterValues, int n){

	// Make an auxiliary array of ints so we can use the int version of GetParameters
	// to do the heavy lifting, and read in the values (initially) as ints
	int *Params = new int[n];
	if(!GetParameters(ParameterStream, ParameterName, Params, n))
		return false;

	// Convert each of the loaded-in ints into unsigned chars
	for(int i = 0; i < n; ++i){
		if(Params[i] <= 0)
			ParameterValues[i] = 0;
		else if(Params[i] >= 255)
			ParameterValues[i] = 255;
		else
			ParameterValues[i] = (unsigned char)(Params[i]);
	}

	delete [] Params;

	return true;
}





/*

	For reading a single unsigned char into ParameterValue from the input deck file connected
	to the stream ParameterStream. If there is no value given in the input deck file, or
	if the label ParameterName is missing altogether, then ParameterValue will be filled
	by DefaultValue.

*/

bool GetParameter (ifstream& ParameterStream, string ParameterName, unsigned char &ParameterValue, unsigned char DefaultValue){

	return GetParameters(ParameterStream, ParameterName, &ParameterValue, 1, DefaultValue);
}





/*

	For reading a single unsigned char into ParameterValue from the input deck file connected
	to the stream ParameterStream. If there is no value given in the input deck file then
	ParameterValue will be filled by DefaultValue and true will be returned.

	If the label ParameterName is missing altogether, then false will be returned.

*/

bool GetParameter2(ifstream& ParameterStream, string ParameterName, unsigned char &ParameterValue, unsigned char DefaultValue){

	return GetParameters2(ParameterStream, ParameterName, &ParameterValue, 1, DefaultValue);
}





/*

	For reading a single unsigned char into ParameterValue from the input deck file connected
	to the stream ParameterStream. If there is no value given in the input deck file, or
	if the label ParameterName is missing altogether, then false will be returned.

*/

bool GetParameter (ifstream& ParameterStream, string ParameterName, unsigned char &ParameterValue){

	return GetParameters(ParameterStream, ParameterName, &ParameterValue, 1);
}





/*

	Special function for loading in a color, labelled by ParameterName, from the input deck file connected
	to the stream ParameterStream. Colors may be specified in the input deck either by a trio of integers in
	the range 0 to 255, or via a string representing one of the color names defined below. The color name
	will be translated into 3 unsigned chars and placed in the array ParameterValues.

*/

const int NumberOfColorsDefined = 30; // Change this if you define any more colors below
static string ColorNames[NumberOfColorsDefined]={
	"White", "Gray", "Grey", "Black", "Red", "Lime", "Blue", "Yellow", "Cyan", "Aqua",
	"Magenta", "Fuchsia", "Silver", "Maroon", "Olive", "Green", "Purple", "Teal",
	"Navy", "Brown", "Coral", "Salmon", "Orange", "Gold", "Turquoise", "Sky", "Indigo",
	"Violet", "Pink", "Cerulean"
};
static unsigned char ColorRGBValues[3*NumberOfColorsDefined] = {
	255, 255, 255, // 0,  White
	128, 128, 128, // 1,  Gray
	128, 128, 128, // 2,  Grey
	  0,   0,   0, // 3,  Black
	255,   0,   0, // 4,  Red
	  0, 255,   0, // 5,  Lime
	  0,   0, 255, // 6,  Blue
	255, 255,   0, // 7,  Yellow
	  0, 255, 255, // 8,  Cyan
	  0, 255, 255, // 9,  Aqua
	255,   0, 255, // 10, Magenta
	255,   0, 255, // 11, Fuchsia
	192, 192, 192, // 12, Silver
	128,   0,   0, // 13, Maroon
	128, 128,   0, // 14, Olive
	  0, 128,   0, // 15, Green
	128,   0, 128, // 16, Purple
	  0, 128, 128, // 17, Teal
	  0,   0, 128, // 18, Navy
	165,  42,  42, // 19, Brown
	255, 127,  80, // 20, Coral
	250, 128, 114, // 21, Salmon
	255, 165,   0, // 22, Orange
	255, 215,   0, // 23, Gold
	 64, 224, 208, // 24, Turquoise
	135, 206, 235, // 25, Sky
	 75,   0, 130, // 26, Indigo
	238, 130, 238, // 27, Violet
	255, 192, 203, // 28, Pink
	 29, 172, 214, // 29, Cerulean
};
bool GetColor(ifstream& ParameterStream, string ParameterName, unsigned char *ParameterValues, string DefaultColor /* = "NoDefault" */){

	// Check if the color is specified via 3 numbers
	if(GetParameters(ParameterStream, ParameterName, ParameterValues, 3))
		return true;

	// Check if the color is specified via a string which will (hopefully) correspond
	// to one of the colors defined below
	string ColorString;
	if(!GetParameter(ParameterStream, ParameterName, ColorString)){
		if(DefaultColor == "NoDefault")
			return false;
		ColorString = DefaultColor;
	}

	for(int i = 0; i < NumberOfColorsDefined; ++i){
		if(CapitalCompare(ColorString, ColorNames[i])){
			ParameterValues[0] = ColorRGBValues[0 + 3*i];
			ParameterValues[1] = ColorRGBValues[1 + 3*i];
			ParameterValues[2] = ColorRGBValues[2 + 3*i];
			return true;
		}
	}

	return false;
}





/*

	Given an array containing unsigned char RGB values in the range 0 to 255, returns a
	string describing the resulting color. If that color is one of the colors defined
	in the global arrays above, returns the corresponding color string. If not, returns
	a string containing the 3 RGB values separated by a space.

*/
string ColorToString(unsigned char *RGB){

	// First check if this is one of the defined colors
	for(int i = 0; i < NumberOfColorsDefined; ++i){
		if(RGB[0] != ColorRGBValues[0 + 3*i] || RGB[1] != ColorRGBValues[1 + 3*i] || RGB[2] != ColorRGBValues[2 + 3*i])
			continue;
		return ColorNames[i];
	}

	// If not, the return string is a triplet of 3-digit RGB values
	return NumberString(int(RGB[0]), 999) + " " + NumberString(int(RGB[1]), 999) + " " + NumberString(int(RGB[2]), 999); 
}

/*

	Determines if the capitalized versions (i.e. all a...z changed to A...Z, with other
	characters unaltered) of the two input strings are identical. If so returns true,
	otherwise false. This is for comparing two strings in a non-case-sensitive manner.

*/

bool CapitalCompare(const string& s1, const string& s2){

	unsigned int n = s1.size();
	if(n != s2.size()) return false;

	char c1, c2;
	for(unsigned int i = 0; i < n; ++i){
		if(s1[i] == s2[i]) continue;
		c1 = s1[i];
		c2 = s2[i];
		if(c1 >= 97 && c1 <= 122) c1 -= 32;
		if(c2 >= 97 && c2 <= 122) c2 -= 32;
		if(c1 != c2) return false;
	}

	return true;
}



 


/*

	Returns the string version of the integer n, padded with leading zeroes
	to bring it up to the length of m. If n is longer than - or equal in length - to m,
	no leading zeroes are added.

*/

string NumberString(int n, int m /* = 0 */){

	stringstream ss;
	ss << n;
	string ReturnString = ss.str();

	int nDigits = HowManyDigits(n);
	int mDigits = HowManyDigits(m);

	if(nDigits >= mDigits)
		return ReturnString;

	for(int i = 0; i < mDigits - nDigits; ++i)
		ReturnString = "0" + ReturnString;

	return ReturnString;
}





/*

	Returns the string version of the double x.

*/

string NumberString(double x){

	stringstream ss;
	ss << x;

	return ss.str();
}





/*

	Returns the number of characters (digits and minus sign, if present) in the integer n. 

*/

int HowManyDigits(int n){

	stringstream ss;
	ss << n;

	return ss.str().size();
}
