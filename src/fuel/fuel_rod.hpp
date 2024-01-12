
#pragma once

#include <map>
#include <unordered_map>
#include <random>

#include "atom.hpp"

namespace sim::fuel
{

class fuel_rod
{
private:

	struct section
	{
		long mass;
		atom_map cluster;
	};

	long mass = 0;
	long energy = 0;
	std::map<int, section> all;
	
	void update_charges(std::map<int, atom_map>& ac);
	void update_neutrons(std::map<int, atom_map>& ac, std::mt19937& gen);
	void update_decays(atom_map& cluster, atom_map& am, std::mt19937& gen, double secs);
	
public:
	
	long extract_energy();
	long calculate_mass() const;
	long get_mass() const;
	long get_mass(int id) const;
	long get_mass(atom a) const;
	long get_mass(int id, atom a) const;
	long remove_mass(double amount, atom_map& dst);
	long remove_mass(int id, double amount, atom_map& dst);
	void add_mass(atom a, long c);
	void add_mass(int id, atom a, long c);
	long add_mass(double amount, atom_map& src);
	long add_mass(int id, double amount, atom_map& src);
	void update(std::mt19937& gen, double secs);
	void display_cluster(int id, int top) const;
	void display_cluster(int id) const;
	void display(int top) const;
	void display() const;
};

}

