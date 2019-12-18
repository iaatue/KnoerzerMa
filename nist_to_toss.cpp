//========================================================================
// Name        : nist_to_toss.cpp
// Author      : Michael Knörzer
// Version     : 1.0 (2018-06-13)
// Copyright   : Copyright (c) 2018
// Description : Reads levels and transitions from NIST
//             : and converts to TOSS-readable format
//             : C++11 !
//========================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <math.h>
using namespace std;

// struct for level and transition
struct level
{
	string name;
	string config;
	string term;
	double energy;
	double J;
	string parity;
};

struct transition
{
	double wvl;
	level low;
	level up;
	double log_gf;
	double gA;
};

int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		cout << "Usage: nist_to_toss <tmad-file>" << endl;
		return 0;
	}

	// buffers for input, in/out stream, line buffer
	vector<level> vec_levels;
	vector<transition> vec_trans;
	ifstream in;
	ofstream out;
	string line;

	// open file and read line by line
	cout << "attempting to open file: " << argv[1] << endl;
	in.open(argv[1]);
	if(in.is_open())
	{
		// read in transitions
		while(getline(in,line))
		{
			// loop through items
			int bars = 0;
			int energies = 0;
			level l_low,l_up;
			transition t;
			string tmp;
			stringstream ss(line);

			while(!ss.eof())
			{
				ss >> tmp;

				// skip ---------
				if(tmp.size() > 10 && "-----" == tmp.substr(1,5))
					break;

				if ("|" == tmp)
				{
					bars++;
					continue;
				}

				std::size_t found;
				// get data
				switch(bars)
				{
					// wavelength
					case 0:
						try
						{
				    		t.wvl = std::stof(tmp);
						}
						catch(std::invalid_argument const &e)
						{
				    		// bad line, skip
				    		cout << "bad line (b=0): " << line << endl;
				    		bars=99;
						}
				    	break;

				    // gA
					case 5:
						try
						{
				    		t.gA = std::stof(tmp);
						}
						catch(std::invalid_argument const &e)
						{
				    		// bad line, skip
				    		cout << "bad line (b=5): " << line << endl;
				    		bars=99;
						}
				    	break;

					// log(gf)
					case 6:
						try
						{
							t.log_gf = std::stof(tmp);
						}
						catch(std::invalid_argument const &e)
						{
							// bad line, skip
							cout << "bad line (b=6): " << line << endl;
							bars=99;
						}
						break;

					// energies
					case 8:
						// case 0: try to get first energy
						// case 1: try to get 2nd energy
						try
						{
							double d = std::stof(tmp);
							//cout << "energy: " << d << endl;
							if(0 == energies)
								l_low.energy = d;
							else if (1 == energies)
								l_up.energy = d;
							else
							{
								// should not happen
								cout << "strange error (b=8): " << line << endl;
								bars=99;
							}
							energies++;
						}
						catch(std::invalid_argument const &e)
						{
							// bad line, skip
						}
						break;

					// 9-11: lower level
					case 9:
						tmp.erase(std::remove(tmp.begin(), tmp.end(), '?'), tmp.end());
						l_low.config = tmp;
						break;
					case 10:
						l_low.term = tmp;
						found = tmp.find('*');
						if (found!=std::string::npos)
							l_low.parity = "o";
						else
							l_low.parity = "e";
						break;
					case 11:
						try
						{
							found = tmp.find('/');
							if (found!=std::string::npos)
							{
								// found a fraction, convert
								l_low.J = std::stof(tmp.substr(0,found)) / std::stof(tmp.substr(found+1));
							}
							else
								l_low.J = std::stof(tmp);

							// make name
							l_low.name = l_low.config + "_" + l_low.term;
						}
						catch(std::invalid_argument const &e)
						{
							// bad line, skip
							cout << "bad J (b=11): " << line << endl;
							bars=99;
						}
						break;

					// 12-14: upper level
					case 12:
						tmp.erase(std::remove(tmp.begin(), tmp.end(), '?'), tmp.end());
						l_up.config = tmp;
						break;
					case 13:
						l_up.term = tmp;
						found = tmp.find('*');
						if (found!=std::string::npos)
							l_up.parity = "o";
						else
							l_up.parity = "e";
						break;
					case 14:
						try
						{
							std::size_t found = tmp.find('/');
							if (found!=std::string::npos)
							{
								// found a fraction, convert
								l_up.J = std::stof(tmp.substr(0,found)) / std::stof(tmp.substr(found+1));
							}
							else
								l_up.J = std::stof(tmp);

							// make name
							l_up.name = l_up.config + "_" + l_up.term;

							// finish
							bars = 50;
						}
						catch(std::invalid_argument const &e)
						{
							// bad line, skip
							cout << "bad J (b=11): " << line << endl;
							bars=99;
						}
						break;

					// otherwise skip
					default:
				    	break;
				}

				// finished
				if(50 == bars)
				{
					// reverse if necessary
					if(l_low.energy > l_up.energy)
					{
						// reverse it
						level ltmp = l_low;
						l_low = l_up;
						l_up = ltmp;
						cout << "Info: levels reversed, check gA/gf for consistency!" << endl;
					}
					// store levels
					vec_levels.push_back(l_low);
					vec_levels.push_back(l_up);

					// store transition
					t.low = l_low;
					t.up = l_up;
					// check
					//cout << "transition finished, " << t.wvl;
					//cout << ", E low: " << t.low.energy << ", E up: " << t.up.energy << endl;
					vec_trans.push_back(t);
					break;
				}

				// on error goto next line
				if(99 == bars)
					break;
			}// end: while(!ss.eof())
		}// end: while(getline(in,line))

		// info
		cout << vec_trans.size() << " transitions found !" << endl;

		// open output file
		out.open((string(argv[1])+"_out_toss").c_str());

		// special line for toss
		out << endl << "  Wavelength         Lower Level         Upper Level   log gf        gA       CF" << endl << endl;

		for(const transition &t : vec_trans)
		{
			out << setw(12) << fixed << setprecision(3) << t.wvl << " "
					<< setw(10) << setprecision(1) << t.low.energy << " (" << t.low.parity << ") " << setw(4) << t.low.J << " "
					<< setw(10) << setprecision(1) << t.up.energy << " (" << t.up.parity << ") " << setw(4) << t.up.J
					<< "  " << setprecision(3) << setw(7) << t.log_gf << " " << scientific << setw(5) << setprecision(3) << t.gA
					<< "    0.000" << endl;
		}

		// sort/unique levels
		cout << endl << "levels: " << vec_levels.size() << endl;
		std::sort(vec_levels.begin(),vec_levels.end(),[](const level& lhs, const level& rhs){return lhs.energy < rhs.energy;});
		vec_levels.erase(std::unique(vec_levels.begin(), vec_levels.end(), [](const level& lhs, const level& rhs){return lhs.energy == rhs.energy;}), vec_levels.end());
		cout << "levels after unique: " << vec_levels.size() << endl;
		cout << endl;

		// output levels
		for(const level &l : vec_levels)
		{
			cout << fixed << setw(9) << setprecision(2);
			cout << l.energy << ": " << l.config << " " << l.term << " (" << l.parity << ") " << setprecision(1) << l.J << endl;
		}
	}
	else
	{
		cout << "ERROR: couldn't open file: " << argv[2] << endl;
	}
	in.close();

	return 0;
}