
#include "fuel_assembly.hpp"

#include <iostream>

using namespace sim::fuel;

fuel_assembly::fuel_assembly(double t, long b)
{
	temperature = t;
	reactivity = 0;

	fr.add_mass(1, atom::from_symbol("B", 10), b);
}

void fuel_assembly::add_mass(atom a, long c)
{
	fr.add_mass(0, a, c);
}

void fuel_assembly::update(std::mt19937& gen, double secs)
{
	fr.update(gen, secs);
	
	temperature += fr.extract_energy();
}

void fuel_assembly::set_reactivity(double r)
{
	if(r > reactivity)
	{
		double a = 1 - (1 - r) / (1 - reactivity);
	
		fr.remove_mass(1, a, neutron_absorber);
	}

	else if(r < reactivity)
	{
		double a = 1 - r / reactivity;

		fr.add_mass(1, a, neutron_absorber);
	}

	reactivity = r;
}

double fuel_assembly::get_temperature()
{
	return temperature;
}

double fuel_assembly::extract_temperature(double outer, double k)
{
	double Q = k * (temperature - outer);

	temperature -= Q;

	return Q;
}

