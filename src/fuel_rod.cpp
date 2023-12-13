
#include "fuel_rod.hpp"
#include "half_life.hpp"
#include "decay.hpp"

#include <map>
#include <iostream>
#include <random>
#include <cmath>
#include <climits>

using namespace sim;

void fuel_rod::add_mass(atom a, long c)
{
	all[a] += c;
}

static void update_decay(fuel_rod& p, atom a, long loop)
{
	atom n = decay::get_mode(a);
	atom o = {0, 0};
	
	if(a.excited && n == atom{-1, 1} && a.neutrons() > 2)
	{
		o = atom{0, 1, true};
		p.add_mass(o, loop);
		
//		std::cout << "decay " << a << " via " << n << " into " << o << " and " << (a - n - o) << std::endl;
	}

	else
	{
//		std::cout << "decay " << a << " via " << n << " into " << (a - n - o) << std::endl;
	}

	p.add_mass(n, loop);
	p.add_mass(a - n - o, loop);
	p.add_mass({0, 0}, loop);
}

static atom do_single_fission(std::mt19937& gen, atom a, int mass_1)
{
	std::uniform_real_distribution<double> dist;
	
	int p = a.protons();
	int n = a.neutrons();
	
	for(int i = 0; i < mass_1; i++)
	{
		double diff = (double)p / (double)(p + n);

		if(dist(gen) > diff)
		{
			n -= 1;
		}

		else
		{
			p -= 1;
		}
	}

	return {a.protons() - p, a.neutrons() - n};
}

static void update_fissile(fuel_rod& p, std::mt19937& gen, atom a, long loop)
{
	std::binomial_distribution<int> dist(a.mass(), (235.0 / 2.0 - 25.0) / 235.0);
	std::uniform_int_distribution<int> neutron_dist(1, 5);
	
	for(long i = 0; i < loop; i++)
	{
		int mass_1 = dist(gen);
		int n = neutron_dist(gen);

		atom a_n = a - atom{0, n};
		atom r1 = do_single_fission(gen, a_n, mass_1);

		while(half_life::get(r1) < -7 || half_life::get(a - r1) < -7)
		{
			r1 = do_single_fission(gen, a_n, mass_1);
		}
		
		atom r2 = a_n - r1;

		r1.excited = true;
		r2.excited = true;

		p.add_mass(r1, 1);
		p.add_mass(a_n - r1, 1);
		p.add_mass({0, 1, true}, n);

//		std::cout << "fissioning " << n << " into " << r1 << ", " << (a - r1) << ", and " << atom{0, 1, true} << " x " << n << std::endl;
	}
	
	p.add_mass({0, 0}, loop);
}

void fuel_rod::update(std::mt19937& gen, double secs)
{
	auto next = all.begin();
	
	for(auto it = all.begin(); it != all.end(); it = next)
	{
		next = it;
		next++;
		auto& [a, c] = *it;

		if(c <= 0)
		{
			all.erase(it);
			continue;
		}

		if(a.mass() == 0) continue;

		long double hl = secs / std::pow(10, half_life::get(a));
		long double p = std::pow(0.5, hl);

		std::binomial_distribution<long> dist(c, 1 - p);
		long loop = dist(gen);

		if(loop == 0) continue;

		if(decay::is_fissile(a))
		{
			update_fissile(*this, gen, a, loop);
		}

		else
		{
			update_decay(*this, a, loop);
		}
		
		c -= loop;
	}
}

long fuel_rod::calculate_mass()
{
	long mass = 0;

	for(auto [a, c] : all)
	{
		mass += a.mass() * c;
	}

	return mass;
}

void fuel_rod::display(int top)
{
	std::multimap<long, atom, std::greater<int>> all_sorted;
	
	for(auto [a, c] : all)
	{
		if(c <= 0) continue;
		all_sorted.insert({c, a});
	}
	
	std::cout << std::endl << "Sample contents:" << std::endl;
	
	int at = 0;

	for(auto [c, a] : all_sorted)
	{
		if(++at > top) break;
		std::cout << "  " << a << ": " << c << std::endl;
	}

	if(at != all_sorted.size())
	{
		std::cout << "... " << (all_sorted.size() - at) << " hidden" << std::endl;
	}

	std::cout << std::endl;
}

void fuel_rod::display()
{
	display(INT_MAX);
}

