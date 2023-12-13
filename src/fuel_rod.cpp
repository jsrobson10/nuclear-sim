
#include "fuel_rod.hpp"
#include "half_life.hpp"
#include "decay.hpp"

#include <iostream>
#include <random>
#include <cmath>

using namespace sim;

void fuel_rod::add_mass(atom a, long c)
{
	all[a] += c;
}

void fuel_rod::update(std::mt19937& gen, double secs)
{
	for(auto& [a, c] : all)
	{
		if(c <= 0 || a.mass() == 0) continue;

		long double hl = secs / std::pow(10, half_life::get(a));
		long double p = std::pow(0.5, hl);

		std::binomial_distribution<long> dist(c, 1 - p);
		long loop = dist(gen);

		if(loop == 0) continue;

		if(decay::is_fissile(a))
		{
			
		}

		else
		{
			atom n = decay::get_mode(a);

			add_mass(n, loop);
			add_mass(a - n, loop);
			add_mass({0, 0}, loop);
		}
		
		c -= loop;
	}
}

void fuel_rod::display()
{
	std::cout << "Sample contents:" << std::endl;

	for(auto [a, c] : all)
	{
		if(c <= 0) continue;

		std::cout << "  " << a << ": " << c << std::endl;
	}

	std::cout << std::endl;
}

