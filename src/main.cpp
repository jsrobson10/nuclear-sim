
#include "fuel/fuel_rod.hpp"
#include <iostream>

int main()
{
	std::random_device rd;
	std::mt19937 gen(rd());

	sim::fuel::fuel_rod fr;
	fr.add_mass(0, sim::fuel::atom::from_symbol("U", 235), 500000);
	fr.add_mass(1, sim::fuel::atom::from_symbol("H", 1), 20000000);
	fr.add_mass(1, sim::fuel::atom::from_symbol("O", 16), 10000000);
	fr.add_mass(2, sim::fuel::atom::from_symbol("B", 10), 10000000);
	fr.add_mass({0, 1}, 1000);

	long mass_1 = fr.calculate_mass();
	
	for(int i = 0; i < 10000; i++)
	{
		fr.update(gen, 1);
		fr.display(20);
	}

	long mass_2 = fr.calculate_mass();

	// verify whether the conservation of mass remains true in this simulation (it should be)
	std::cout << "initial mass is " << mass_1 << std::endl;
	std::cout << "final mass is: " << mass_2 << std::endl;
	std::cout << "conservation of mass is " << (mass_1 == mass_2 ? "maintained." : "VIOLATED.") << std::endl;
}

