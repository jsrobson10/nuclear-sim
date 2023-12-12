
#include "fuel_rod.hpp"
#include <iostream>

int main()
{
	std::random_device rd;
	std::mt19937 gen(rd());

	sim::fuel_rod fr;
	fr.add_mass(sim::atom::from_mass(88, 207), 1000000000000);

	for(int i = 0; i < 10000; i++)
	{
		fr.display();
		fr.update(gen, 1);
	}
}

