//========================================================================
// Name        : tmad_to_grotrian.cpp
// Author      : Michael Knörzer
// Version     : 1.0 (2018-12-03)
// Copyright   : Copyright (c) 2018
// Description : Reads in levels + lines in TMAP format and
//             : creates a Grotrian diagram
//             : C++11 !
//========================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <math.h>
using namespace std;

// struct for levels
struct level
{
	int mult;
	int n;		// config
	int p;		// config/term, should be the same
	int l;		// term
	double energy;
	double J;

	string name;
	string conf;
	string L;
	string term;
	string parity;

	// compare to any A10-string (name)
	bool operator==(const string& n)
	{
	   return (name == n);
	}
};
// different sorting
struct
{
	bool operator()(const level& i, const level& j)
	{
		if(i.l == j.l)
			return i.energy < j.energy;
		return i.l < j.l;
	}
}sort_by_L;
struct
{
	bool operator()(const level& i, const level& j)
	{
		return i.energy < j.energy;
	}
}sort_by_E;

// struct to keep track of top labels
struct mul_lp
{
	int mult;
	int l;
	int p; // parity

	// sorty by l, parity
	bool operator()(const mul_lp& i, const mul_lp& j)
	{
		if(i.l == j.l)
			return i.p < j.p;
		return i.l < j.l;
	}
}sort_by_lp;
// compare 2 items of same multiplicity
bool operator==(const mul_lp& lhs, const mul_lp& rhs )
{
   return (lhs.l == rhs.l) && (lhs.p == rhs.p);
}

// collection of levels for one type of multiplicity
struct levels_mult
{
	int mult;
	vector<mul_lp> multis;
	vector<level> levels;

	bool operator==(const int& m)
	{
	   return (mult == m);
	}
	bool operator()(const levels_mult& i, const levels_mult& j)
	{
		return i.mult < j.mult;
	}
}sort_by_mult;

auto sum_labels = [](int lhs, const levels_mult& rhs)
{
	return lhs + rhs.multis.size();
};


// struct for lines
struct transition
{
	level low,up;
	double wvl;
	double gf;
	double gA;
	string name;

	// sort by wavelength
	bool operator()(const transition& i, const transition& j)
	{
		if(i.wvl == j.wvl)
			return i.gf < j.gf;
		return i.wvl < j.wvl;
	}
}sort_by_wvl;
struct
{
	bool operator()(const transition& i, const transition& j)
	{
		if(i.gf == j.gf)
			return i.wvl < j.wvl;
		return i.gf < j.gf;
	}
}sort_by_gf;

// functions to convert angular momentum letters to numbers
// and vice versa
int det_L(char c);
string get_L(int l);

// enumerate different states
enum state {SEARCH_ATOM, READ_ATOM, SEARCH_CONTENT, READ_LEVELS, READ_RBB};

int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		cout << endl << "Usage: tmad_to_grotrian <TMAD file> <options>" << endl;
		cout << endl << "Options: e=<number>, n=<number>, l=<number>, c=<Term><parity>" << endl;
		cout << "Exclude levels/configurations from the diagram which have" << endl;
		cout << "energy >= e, principal quantum number >= n, angular momentum qn >= l" << endl;
		cout << "or which have a certain configuration i.e. 3Po or 4Se" << endl;
		return 0;
	}

	// default values
	double skip_e = 9.9e+30;
	int skip_n = 26;
	int skip_l = 23;
	vector<string> skip_conf;
	double ionlimit;

	// get all options, start with arg #3
	for(int i = 3; i<argc; i++)
	{
		string s(argv[i]);
		if (s.substr(0,2) == "e=")
		{
			stringstream ss(s.substr(2));
			ss >> skip_e;
		}
		else if (s.substr(0,2) == "n=")
		{
			stringstream ss(s.substr(2));
			ss >> skip_n;
		}
		else if (s.substr(0,2) == "l=")
		{
			stringstream ss(s.substr(2));
			ss >> skip_l;
		}
		else if (s.substr(0,2) == "c=")
		{
			skip_conf.push_back(s.substr(2));
		}
	}

	// buffers for input, in/out stream, line buffer
	vector<level> vec_levels;
	vector<transition> vec_lines;
	ifstream in;
	string line;
	// different multiplicities
	vector<levels_mult> all_multiplets;

	state s = SEARCH_ATOM;
	string atom;
	int alen;
	int charge;
	int total;
	double unit;
	const double c = 2.99792458e10;	// cms/s 3*10^10

	// open file and read line by line
	cout << "** attempting to open file: " << argv[1] << endl;
	in.open(argv[1]);
	if(in.is_open())
	{
		// read in transitions
		while(getline(in,line))
		{
			// skip comments
			if(line.substr(0,1) == ".")
				continue;

			stringstream ss(line);
			transition tr;
			level le;
			double g;
			double eHz;
			bool skip = false;

			// read in file according to status s
			switch(s)
			{
			case SEARCH_ATOM:
				if(line.substr(0,4) == "ATOM")
					s = READ_ATOM;
				break;

			case READ_ATOM:
				ss >> atom;
				ss >> charge;
				if(atom.size() == 2)
				{
					// if the element has 2 letters already, the total
					// element code should always be 3 characters long
					alen = 3;
				}
				else
				{
					ss.str("");
					ss.clear();
					ss << atom << (charge+1);
					atom = ss.str();
					alen = atom.size();
				}
				s = SEARCH_CONTENT;
				break;

			case SEARCH_CONTENT:
				if(line == "L")
					s = READ_LEVELS;
				else if(line == "LTE")
					s = READ_LEVELS;
				else if(line == "RBB")
					s = READ_RBB;
				break;

			case READ_LEVELS:
				if(line == "0")
				{
					s = SEARCH_CONTENT;
					continue;
				}
				le.name = line.substr(0,10);
				// first 7 characters are atom + config
				// 8-10 is the term
				ss.str(line.substr(alen,7-alen));
				ss.clear();
				ss >> le.conf;
				std::transform(le.conf.begin(), le.conf.end(), le.conf.begin(), ::tolower);

				// get n from conf
				ss.str(le.conf.substr(0,2));
				ss.clear();
				ss >> le.n;

				ss.str(line.substr(7,3));
				ss.clear();
				ss >> le.term;
				// get parity and correct term if necessary
				if(le.term.size() == 2)
				{
					// even, term size is fine already
					le.p = 0;
					le.parity = "e";
				}
				else if(le.term.size() == 3)
				{
					// must be odd, remove parity from term string
					le.p = 1;
					le.parity = "o";
					le.term = le.term.substr(0,2);
				}
				else
				{
					// error
					cout << "** error with level parity:" << endl << "** " << line << endl;
					continue;
				}
				// convert multiplicity to int, convert L to int
				ss.str(le.term.substr(0,1));
				ss.clear();
				ss >> le.mult;
				if(le.mult < 1 || le.mult > 9)
				{
					cout << "** Error with multiplicity:" << endl << line << endl;
					continue;
				}
				le.L = le.term.substr(1,1);
				le.l = det_L(le.L.c_str()[0]);
				if(le.l < 0)
				{
					cout << "** Error with total angular momentum L:" << endl << line << endl;
					continue;
				}

				// get the remainder of the line
				ss.str(line.substr(20));
				ss.clear();
				ss >> eHz;
				// check if we have determined the ionization limit yet
				// ground state level should be the first so determine it from there
				ss >> g;
				le.J = (g-1)/2;
				if(vec_levels.size() == 0)
				{
					ionlimit = eHz / c;
				}
				// convert energy back to cm^-1
				le.energy = ionlimit - (eHz / c);

				// check if we should skip this e,p, or l
				if(le.energy >= skip_e || le.n >= skip_n || le.l >= skip_l)
					continue;

				// check if should skip this term
				ss.str(le.term);
				ss.clear();
				ss << le.parity;
				for(const auto &x:skip_conf)
				{
					if(x == ss.str())
					{
						skip = true;
						break;
					}
				}
				if(skip)
				{
					continue;
				}

				// all good -> add to vector
				vec_levels.push_back(le);
				break;

			case READ_RBB:
				if(line == "0")
				{
					s = SEARCH_CONTENT;
					break;
				}
				// check if we find both levels
				auto it = find(vec_levels.begin(),vec_levels.end(),line.substr(0,10));
				if(it == vec_levels.end())
					continue;
				else
					tr.low = *it;
				it = find(vec_levels.begin(),vec_levels.end(),line.substr(10,10));
				if(it == vec_levels.end())
					continue;
				else
					tr.up = *it;

				tr.name = line;
				tr.wvl = pow(10.0,8.0) / (tr.up.energy - tr.low.energy);
				ss.str(line.substr(20));
				ss.clear();

				// should be like " 3 3 f_ik"
				ss >> tr.gf;
				ss >> tr.gf;
				ss >> tr.gf;
				// gf = g_low * f_ik
				tr.gf = (tr.low.J * 2 + 1) * tr.gf;
				tr.gA = tr.gf / 1.49919E-16 / tr.wvl / tr.wvl;

				vec_lines.push_back(tr);
				break;
			}
			// end switch
		}
		// end getline

		// check if we have found any levels
		if(vec_levels.size() < 1)
		{
			cout << "** found no levels **" << endl;
			return 0;
		}

		// determine different terms and sort
		sort(vec_levels.begin(),vec_levels.end(),sort_by_E);
		for(const auto &i:vec_levels)
		{
			// check if we had this multiplicity before
			auto it = find (all_multiplets.begin(), all_multiplets.end(), i.mult);
			if (it == all_multiplets.end())
			{
				// does not exist yet
				vector<level> v;
				v.push_back(i);
				vector<mul_lp> m;
				all_multiplets.push_back({i.mult,m,v});
			}
			else
			{
				// add to correct collection
				it->levels.push_back(i);
			}
		}
		sort(all_multiplets.begin(), all_multiplets.end(), sort_by_mult);

		// # of total columns
		// one for each mlp combination
		// one for each separator between groups (== all_multis size - 1)
		// + 1 extra for space at left/right border --> all_multis.size + size of all mlp_vectors
		total = all_multiplets.size();
		for(auto &i:all_multiplets)
		{
			// add for top labels
			sort(i.levels.begin(),i.levels.end(),sort_by_L);
			for(const auto &j:i.levels)
			{
				mul_lp mlp;
				mlp.mult = j.mult;
				mlp.l = j.l;
				mlp.p = j.p;
				// check if we had this mult/L/P before
				auto it = find(i.multis.begin(),i.multis.end(),mlp);
				if(it == i.multis.end())
					i.multis.push_back(mlp);
			}
			sort(i.multis.begin(),i.multis.end(),sort_by_lp);
			total += i.multis.size();
		}

		// get x-positions / width of one column
		unit = 100.0 / total;

		// energy low/high
		double low,high;
		low = vec_levels.front().energy;
		high = vec_levels.back().energy;
		double yoffset = (ionlimit * 0.02);


		// make plot
		cout << fixed << setprecision(2);
		cout << endl << "PAPERFORMAT A3Q" << endl;
		cout << "MULTIPLOT START" << endl;
		cout << "** y min/max: " << low << "/" << high << endl;
		cout << "** y offset: " <<  yoffset << endl << endl;

		cout << "PLOT: labels" << endl;
		cout << "\\OFS 2.0 2.0" << endl;
		cout << "\\INBOX" << endl;
		cout << "\\PEN 1" << endl;
		cout << "\\FONT=HELVET" << endl;
		cout << "\\LETTERSIZE=0.25" << endl;
		cout << "\\NOCOPYRIGHT" << endl;
		cout << "\\LUN 50.0 " << (ionlimit+2*yoffset)/1000*1.03 << " -2.9 0.0 0.30 Grotrian diagram of " << argv[1] << endl;
		cout << "HEADER :\\CENTER\\" << endl;
		cout << "X-ACHSE:\\CENTER\\" << endl;
		cout << "Y-ACHSE:\\CENTER\\ energy / 1000 cm&H-1&M" << endl;
		cout << "    MASSTAB       MINIMUM       MAXIMUM    TEILUNGEN     BESCHRIFT.    DARUNTER" << endl;
		cout << "X: 38.00CM              0.0         100.0         10.0          10            0.0 NOLAB NOTICK-BOTH" << endl;
		cout << "Y: 25.70CM            " << (-yoffset)/1000 << "        " << (ionlimit+2*yoffset)/1000 << "         ";
		cout << (ionlimit < 1.0e+6 ? 10 : (ionlimit < 8.0e+6 ? 50 : (ionlimit < 16.0e+6 ? 100 : 200))) << "           ";
		cout << (ionlimit < 1.0e+6 ? 100 : (ionlimit < 8.0e+6 ? 500 : (ionlimit < 16.0e+6 ? 1000 : 2000))) << "            0.0" << endl;
		cout << "N=  ?  PLOTSYMBOL 9 SYMBOLSIZE 0.1 PEN 1 XYTABLE SELECT 1 2 COLOR=1" << endl;
		cout << "FINISH" << endl;
		cout << "END" << endl << endl;

		cout << "PLOT: Grotrian Diagram of TMAD File: " << argv[1] << endl;
		cout << "\\OFS 2.0 2.0" << endl;
		cout << "\\INBOX" << endl;
		cout << "\\PEN 1" << endl;
		cout << "\\FONT=HELVET" << endl;
		cout << "\\LETTERSIZE=0.25" << endl;
		cout << "\\NOCOPYRIGHT" << endl;
		cout << "HEADER :\\CENTER\\" << endl;
		cout << "X-ACHSE:\\CENTER\\" << endl;
		cout << "Y-ACHSE:\\CENTER\\" << endl;
		cout << "    MASSTAB       MINIMUM       MAXIMUM    TEILUNGEN     BESCHRIFT.    DARUNTER" << endl;
		cout << "X: 38.00CM              0.0         100.0         10.0          10            0.0 NOTICK-BOTH" << endl;
		cout << "Y: 25.70CM         " << (-yoffset) << "      " << (ionlimit+2*yoffset) << "      10000        100000            0.0 NOTICK-BOTH" << endl;
		cout << "N=  ?  PLOTSYMBOL 9 SYMBOLSIZE 0.1 PEN 1 XYTABLE SELECT 1 2 COLOR=1" << endl;
		cout << "0 " << ionlimit << endl;
		cout << "100 " << ionlimit << endl;
		cout << "FINISH" << endl;
		cout << "** ionization limit: " << ionlimit << endl << endl;

		stringstream ss_levels,ss_labels,ss_top,ss_seps;
		ss_levels << fixed << setprecision(2);
		ss_labels << fixed << setprecision(2);
		ss_top << fixed << setprecision(2);
		ss_seps << fixed << setprecision(1);

		int before = 0;
		for(const auto &i:all_multiplets)
		{
			// separators
			// 0.5 units space left side + #units before + #units current
			int width = i.multis.size();
			double xpos = unit*(0.5 + before + width + 0.5);
			if(xpos < 100)
				ss_seps << "\\LINUN " << xpos << " YMIN " << xpos << " YMAX 0.0 0.0 SIZE=0.1 SYMBOL=9" << endl;
			// label, centered in that area
			ss_seps << "\\LUN " << (unit*before + (unit*width*0.5) + unit*0.5) << " " << ionlimit + yoffset*0.4 << " -0.2 0.0 0.20 S=" << ((i.mult - 1.0)*0.5) << endl;

			// top labels
			cout << fixed << setprecision(3);
			ss_labels << fixed << setprecision(3);
			int top_offset = 0;
			for(const auto &j:i.multis)
			{
				// todo: percentage instead of 0.4?
				double xlabelpos = unit*(before + top_offset + 0.5 + 0.4);
				ss_top << "\\LUN " << xlabelpos << " YMAX 0.000 0.080 0.2 " << "&H" << j.mult << "&M" << get_L(j.l) << (j.p == 0 ? "" : "&Ho&M") << endl;
				top_offset++;
			}
			// levels
			for(const auto &j:i.levels)
			{
				// find position
				mul_lp lp;
				lp.l = j.l;
				lp.p = j.p;
				auto it = find(i.multis.begin(),i.multis.end(),lp);
				int pos = it - i.multis.begin() + before;
				double xlevelpos = unit*(pos+0.5+0.5);
				// actual level - offset to the left
				ss_levels << "\\LINUN " << (xlevelpos-unit*0.3) << " " << j.energy << " " << (xlevelpos) << " " << j.energy << " 0.0 0.0" << endl;
				// label next to level
				ss_labels << "\\LUN " << (xlevelpos+unit*0.1) << " " << j.energy << " -0.0 -0.05 0.10 " << j.conf << endl;
			}
			// increase offset
			before+= width + 1;
		}

		// check if we have found any lines
		if(vec_lines.size() < 1)
		{
			cout << "** found no lines **" << endl;
		}
		else
		{
			stringstream sslines;
			sslines << fixed << setprecision(2);
			for(const auto &i:vec_lines)
			{
				// connecting line between levels
				// find position
				auto getxpos = [&](level lev)
				{
					mul_lp lp;
					lp.l = lev.l;
					lp.p = lev.p;
					int ibefore = 0;
					int res = 0;
					for(const auto &j:all_multiplets)
					{
						if(lev.mult == j.mult)
						{
							auto it = find(j.multis.begin(),j.multis.end(),lp);
							res = (it - j.multis.begin() + ibefore);
							break;
						}
						else
						{
							// increase by size of that area + separator
							ibefore += j.multis.size() + 1;
						}
					}
					return res;
				};

				// width of a level = 0.3 units, position = xpos + 0.5 + 0.5
				// --> 0.15 to the left from the right end of level line
				double lowpos = unit*(getxpos(i.low)+0.85);
				double highpos = unit*(getxpos(i.up)+0.85);
				sslines << "\\LINUN " << (lowpos) << " " << i.low.energy <<  " " << (highpos) << " " << (i.up.energy) << " 0.0 0.0" << endl;
			}
			cout << "** connecting lines: **" << endl;
			cout << "\\DEFINECOLOR 9 0.6 0.6 0.6" << endl;
			cout << "\\PEN=1" << endl;
			cout << "\\COLOR=9" << endl;
			cout << sslines.str();
			cout << "\\COLOR=1" << endl;
			cout << "** total # lines: " << vec_lines.size() << " " << endl;
			cout << "** end connecting lines **" << endl << endl;
		}

		cout << "** start levels **" << endl;
		cout << "\\PEN=1" << endl;
		cout << "\\COLOR=1" << endl;
		cout << ss_levels.str();
		cout << "** total # levels: " << vec_levels.size() << " " << endl;
		cout << "** end levels **" << endl << endl;

		cout << "** start inside labels **" << endl;
		cout << "\\COLOR=3" << endl;
		cout << ss_labels.str();
		cout << "\\COLOR=1" << endl;
		cout << "** total # inside labels: " << vec_levels.size() << " " << endl;
		cout << "** end inside labels **" << endl << endl;

		cout << "** start top labels **" << endl;
		cout << "\\PEN=5" << endl;
		cout << "\\COLOR=1" << endl;
		cout << ss_top.str();
		cout << "\\PEN=1" << endl;
		cout << "** total # top labels: " << std::accumulate(all_multiplets.begin(), all_multiplets.end(),0,sum_labels) << " " << endl;
		cout << "** end top labels **" << endl << endl;

		cout << "** start separators ** " << endl;
		cout << ss_seps.str();
		cout << "** end separators ** " << endl << endl;

		cout << "END" << endl << "MULTIPLOT END" << endl << endl;
		in.close();
	}
	else
	{
		cout << "Could not open file: " << argv[1] << endl;
		return 0;
	}

	// end
	return 0;
}

int det_L(char c)
{
	int l = -1;
	switch(c)
	{
	case 's':
	case 'S':
		l = 0;
		break;
	case 'p':
	case 'P':
		l = 1;
		break;
	case 'd':
	case 'D':
		l = 2;
		break;
	case 'f':
	case 'F':
		l = 3;
		break;
	case 'g':
	case 'G':
		l = 4;
		break;
	case 'h':
	case 'H':
		l = 5;
		break;
	case 'i':
	case 'I':
		l = 6;
		break;
	case 'k':
	case 'K':
		l = 7;
		break;
	case 'l':
	case 'L':
		l = 8;
		break;
	case 'm':
	case 'M':
		l = 9;
		break;
	case 'n':
	case 'N':
		l = 10;
		break;
	case 'o':
	case 'O':
		l = 11;
		break;
	case 'q':
	case 'Q':
		l = 12;
		break;
	case 'r':
	case 'R':
		l = 13;
		break;
	case 't':
	case 'T':
		l = 14;
		break;
	case 'u':
	case 'U':
		l = 15;
		break;
	case 'v':
	case 'V':
		l = 16;
		break;
	case 'w':
	case 'W':
		l = 17;
		break;
	case 'x':
	case 'X':
		l = 18;
		break;
	case 'y':
	case 'Y':
		l = 19;
		break;
	case 'z':
	case 'Z':
		l = 20;
		break;
	default:
		break;
	}
	return l;
}

string get_L(int l)
{
	switch(l)
	{
	case 0:
		return "S";
	case 1:
		return "P";
	case 2:
		return "D";
	case 3:
		return "F";
	case 4:
		return "G";
	case 5:
		return "H";
	case 6:
		return "I";
	case 7:
		return "K";
	case 8:
		return "L";
	case 9:
		return "M";
	case 10:
		return "N";
	case 11:
		return "O";
	case 12:
		return "Q";
	case 13:
		return "R";
	case 14:
		return "T";
	case 15:
		return "U";
	case 16:
		return "V";
	case 17:
		return "W";
	case 18:
		return "X";
	case 19:
		return "Y";
	case 20:
		return "Z";

	default:
		cout << "** problem converting l: " << l << endl;
		return "?";
	}
}