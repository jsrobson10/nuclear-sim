
#include "fuel_rod.hpp"
#include <iostream>

int main()
{
	std::random_device rd;
	std::mt19937 gen(rd());

	sim::fuel_rod fr;
	fr.add_mass(sim::atom::from_symbol("U", 235), 1000000);

	for(int i = 0; i < 100000; i++)
	{
		fr.update(gen, 1000000000);
		fr.display(40);
	}

	// verify whether the conservation of mass remains true in this simulation (it should be)
	std::cout << "final mass is: " << fr.calculate_mass() << std::endl;
}

