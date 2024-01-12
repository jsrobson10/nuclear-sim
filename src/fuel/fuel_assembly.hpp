
#pragma once

#include <ostream>

#include "fuel_rod.hpp"
#include "atom.hpp"

namespace sim::fuel
{


class fuel_assembly
{
	fuel_rod fr;
	atom_map neutron_absorber;
		
	double temperature;
	double reactivity;

public:

	fuel_assembly(double t, long b);
	void add_mass(atom a, long c);
	void update(std::mt19937& gen, double secs);
	void set_reactivity(double r);
	double get_temperature();
	double extract_temperature(double outer, double k);

	friend std::ostream& operator <<(std::ostream& o, const fuel_assembly& fa)
	{
		o << "\nFuel Assembly:\n\n";
		o << "reactivity: " << fa.reactivity << "\n";
		o << "rod mass: " << fa.fr.get_mass(1) << "\n";
		o << "temperature: " << fa.temperature << "\n";
		o << "neutrons: " << fa.fr.get_mass(0, atom(0, 1)) << "\n";
		o << "fuel: " << fa.fr.get_mass(0, atom::from_symbol("U", 235)) << "\n";
		o << "mass: " << fa.fr.get_mass(0) << "\n";

		return o;
	}

};

}

