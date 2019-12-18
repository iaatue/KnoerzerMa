//========================================================================
// Name        : adamant_to_toss.cpp
// Author      : Michael Knörzer
// Version     : 1.0 (2018-06-08)
// Copyright   : Copyright (c) 2018
// Description : Reads levels and transitions from adamant
//             : and converts to TOSS-readable format
//========================================================================
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <math.h>

using std::string;
using std::cout;
using std::endl;


// struct for level and transition
struct level
{
	level(): id(0), config(""), energy(0.0), J(0.0), parity("?"){}
	level(int _id, string _config, double _energy, double _J, string _parity):
		id(_id), config(_config), energy(_energy), J(_J), parity(_parity){}

	int id;
	string config;
	double energy;
	double J;
	string parity;
};

struct transition
{
	transition(double _wvl, level _low, level _up, double _loggf, double _gA): wvl(_wvl), low(_low), up(_up), loggf(_loggf), gA(_gA){}

	double wvl;
	level low;
	level up;
	double loggf;
	double gA;
};

// program start
int main(int argc, char* argv[])
{
	if(argc < 3)
	{
		cout << "Usage: adamant_to_toss <level-file> <line-file>" << std::endl;
		return 0;
	}

	// buffers for input, in/out stream, line buffer
	std::vector<level> vec_levels;
	std::vector<transition> vec_lines;
	std::ifstream in;
	std::stringstream out;
	string line;

	// read level file
	cout << "** attempting to open file: " << argv[1] << endl;
	in.open(argv[1]);
	if(in.is_open())
	{
		// read in levels
		while(getline(in,line))
		{
			std::stringstream s(line);
			int id;
			string P,conf;
			double E,J;
			s >> id >> E >> J >> P >> conf >> conf;

			// add to vector
			vec_levels.push_back(level{id,conf,E,J,P});
		}
		in.close();
	}
	else
	{
		std::cout << "** ERROR: couldn't open file: " << argv[1] << std::endl;
		return -1;
	}

	// read line file
	std::cout << "** attempting to open file: " << argv[2] << std::endl;
	in.open(argv[2]);
	if(in.is_open())
	{
		out << std::endl << "  Wavelength         Lower Level         Upper Level   log gf        gA" << std::endl << std::endl;
		// read in transitions
		while(getline(in,line))
		{
			std::stringstream s(line);
			int id_low, id_up;
			string s1;
			double wvl, gf, A;
			s >> id_low >> s1 >> id_up >> s1 >> s1 >> wvl >> A >> gf;

			level l_low,l_up;
			int found = 0;

			// check / assign level references
			for(const level &elem : vec_levels)
			{
				if(found == 2)
					break;

				if(id_low == elem.id || id_up == elem.id)
				{
					found++;
					if(1 == found)
						l_low = elem;
					else
					{
						if(l_low.energy > elem.energy)
						{
							l_up = l_low;
							l_low = elem;
						}
						else
							l_up = elem;
					}
				}
			}

			if(found != 2)
			{
				cout << "** Error: couldn't find corresponding levels to " << endl << line << endl;
				return 1;
			}

			// add to vector
			vec_lines.push_back(transition{wvl, l_low, l_up, log10(gf), A * (2*l_up.J + 1)});
		}
		in.close();

		// sort lines
		std::sort(vec_lines.begin(),vec_lines.end(),[](const transition& lhs, const transition& rhs){return lhs.wvl < rhs.wvl;});
		// prepare / output
		for(const transition &t : vec_lines)
		{
			out << std::setw(12) << std::fixed << std::setprecision(3) << t.wvl << " "
					<< std::setw(10) << std::setprecision(1) << t.low.energy << " (" << t.low.parity << ") " << std::setw(4) << t.low.J << " "
					<< std::setw(10) << std::setprecision(1) << t.up.energy << " (" << t.up.parity << ") " << std::setw(4) << t.up.J
					<< "  " << std::setprecision(3) << std::setw(7) << t.loggf << " " << std::scientific << std::setw(5) << std::setprecision(3) << t.gA << endl;
		}
		cout << out.str();
	}
	else
	{
		cout << "** ERROR: couldn't open file: " << argv[2] << endl;
		return -1;
	}

	return 0;
}