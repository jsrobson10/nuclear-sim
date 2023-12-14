
#include "fuel_rod.hpp"
#include <iostream>

int main()
{
	std::random_device rd;
	std::mt19937 gen(rd());

	sim::fuel_rod fr;
	fr.add_mass(sim::atom::from_symbol("U", 235), 300000);
	fr.add_mass(sim::atom::from_symbol("U", 238), 1000000);
	fr.add_mass({0, 1, true}, 1000);

	long mass_1 = fr.calculate_mass();

	for(int i = 0; i < 10000; i++)
	{
		fr.update(gen, 1);
		fr.display(40);
	}

	long mass_2 = fr.calculate_mass();

	// verify whether the conservation of mass remains true in this simulation (it should be)
	std::cout << "initial mass is " << mass_1 << std::endl;
	std::cout << "final mass is: " << mass_2 << std::endl;
	std::cout << "conservation of mass is " << (mass_1 == mass_2 ? "" : "NOT ") << "maintained.\n";
}

