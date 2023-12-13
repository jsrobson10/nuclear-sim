
#pragma once

#include <unordered_map>
#include <random>

#include "atom.hpp"

namespace sim
{

class fuel_rod
{
private:

	long mass = 0;
	std::unordered_map<atom, long> all;
	
	void update_neutrons(std::mt19937& gen);
	
public:
	
	long calculate_mass();
	void add_mass(atom a, long c);
	void update(std::mt19937& gen, double secs);
	void display(int top);
	void display();
};

}

