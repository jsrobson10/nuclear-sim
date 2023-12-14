
#pragma once

#include <unordered_map>
#include <random>

#include "atom.hpp"

namespace sim
{

typedef std::unordered_map<atom, long> atom_map;

class fuel_rod
{
private:

	long mass = 0;
	atom_map all;
	
	void update_neutrons(atom_map& am, std::mt19937& gen);
	void update_decays(atom_map& am, std::mt19937& gen, double secs);
	
public:
	
	long calculate_mass();
	void add_mass(atom a, long c);
	void update(std::mt19937& gen, double secs);
	void display(int top);
	void display();
};

}

