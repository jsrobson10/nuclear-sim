
#include "decay.hpp"
#include "half_life.hpp"

#include <climits>
#include <iostream>

using namespace sim::fuel;

static const atom FISSILE_TYPES[] = {
	atom::from_mass(92, 233),
	atom::from_mass(92, 235),
	atom::from_mass(94, 239),
	atom::from_mass(94, 241)
};

/*std::array<atoms, 5> decay::calc_for_types(std::mt19937& gen, atom a)
{
	std::array<atoms, 5> atoms_arr;
	std::array<double 5> weights;

	double weight_at = 0;

	for(int i = 0; i < decay::TYPES.size(); i++)
	{
		weight_at += half_life::get(a - decay::TYPES[i]);
		weights[i] = weight_at;

		std::cout << 
	}

	return atoms_arr;
}*/

atom decay::get_mode(atom a)
{
	atom stable_type{0, 0};
	double stable_hl = LONG_MIN;

	for(atom type : decay::TYPES)
	{
		double hl = half_life::get(a - type);
		
		if(hl > stable_hl)
		{
			stable_hl = hl;
			stable_type = type;
		}
	}


	return stable_type;
}

bool decay::is_fissile(atom a)
{
	for(atom type : FISSILE_TYPES)
	{
		if(type == a)
		{
			return true;
		}
	}

	return false;
}

