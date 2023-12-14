
#pragma once

#include <map>
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
	std::map<int, atom_map> all;
	
	void update_neutrons(std::map<int, atom_map>& ac, std::mt19937& gen);
	void update_decays(atom_map& cluster, atom_map& am, std::mt19937& gen, double secs);
	
public:
	
	long calculate_mass();
	void add_mass(atom a, long c);
	void add_mass(int id, atom a, long c);
	void update(std::mt19937& gen, double secs);
	void display_cluster(int id, int top);
	void display_cluster(int id);
	void display(int top);
	void display();
};

}

