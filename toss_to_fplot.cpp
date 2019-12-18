//========================================================================
// Name        : toss_to_fplot.cpp
// Author      : Michael Knörzer
// Version     : 1.0 (2019-03-30)
// Copyright   : Copyright (c) 2019
// Description : Transforms lines from TOSS format (wvl+log gf) into
//               WRPLOT idents to be used in an f over lambda plot
//========================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <tuple>
#include <vector>
#include <algorithm>
using namespace std;

int main(int argc, char* argv[])
{
	ifstream in;
	string line;
	// custom struct vector: <wavelength, f-value, loggf>
	std::vector<std::tuple<double, double, string>> values;
	double scale = 1.0;
	bool asUnit = false;

	if(argc < 2)
	{
		cout << "Transforms lines in TOSS format (wvl+log gf) into" << endl;
		cout << "WRPLOT idents to use in a f over lambda plot" << endl << "------------------------------------------------" << endl;
		cout << "Usage: toss_to_fplot <filename> <scalefactor=1.0> <u=false>" << endl;
		return(0);
	}
	else if(argc >= 3)
	{
		// get scale factor
		stringstream ss(argv[2]);
		ss >> scale;
		cout << "** scale factor: " << scale << endl;

		// get true/false for U/cm
		if(argc >= 4)
		{
			stringstream ss2(argv[3]);
			ss2 >> std::boolalpha >> asUnit;
		}
		cout << "** output units (cm/U): " << (asUnit ? "U" : "cm") << endl;
	}

	// open file and read line by line
	cout << "** attempting to open file: " << argv[1] << endl;
	cout << fixed << setprecision(4);
	in.open(argv[1]);
	if(in.is_open())
	{
		while(getline(in,line))
		{
			double wvl,j_low,j_up,loggf,gA,f;
			int g_low;
			string dump;
			stringstream ss(line);

			ss >> wvl >> dump >> dump >> j_low >> dump >> dump >> j_up >> loggf >> gA;
			ss.clear();
			g_low = (2*j_low)+1;

			f = pow(10,loggf) / g_low;

			// check if f-value and gA deviate
			double f2 = gA * 1.49919E-16 * wvl * wvl / g_low;
			double ratio = f/f2;
			double diff = abs(1-ratio);
			if(diff > 0.5)
			{
				// error: print out some useful information
				cout << "** deviating f-value/gA found:\n";
				cout << "*** " << wvl << " gA:" << gA << " f:" << f << " f2:" << f2 << " ratio:" << ratio << " diff:" << diff << endl;
				cout << "*** jlow:" << j_low << " glow:" << g_low << endl;
			}
			else
				values.push_back(std::make_tuple(wvl, f, to_string(loggf)));
		}

		// sort by first value, i.e., wavelength
		std::sort(values.begin(),values.end());
		for(const auto &v:values)
		{
			// ID length by units or cm
			cout << "\\IDLENG " << std::get<1>(v) * scale << (asUnit ? "U" : "") << endl;
			cout << "\\IDENT  " << std::get<0>(v) << "    " << std::get<2>(v) << endl;
		}
	}
	return 0;
}