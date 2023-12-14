
#include "fuel_rod.hpp"
#include "half_life.hpp"
#include "decay.hpp"

#include <set>
#include <map>
#include <iostream>
#include <random>
#include <cmath>
#include <climits>

using namespace sim;


void fuel_rod::add_mass(atom a, long c)
{
	all[a] += c;
	mass += c * a.mass();
}

static void update_decay(atom_map& am, atom a, long loop)
{
	atom n = decay::get_mode(a);
	atom o = {0, 0};
	
	if(a.excited && n == atom{-1, 1} && a.neutrons() > 2)
	{
		o = atom{0, 1, true};
		am[o] += loop;
		
//		std::cout << "decay " << a << " via " << n << " into " << o << " and " << (a - n - o) << std::endl;
	}

	else
	{
//		std::cout << "decay " << a << " via " << n << " into " << (a - n - o) << std::endl;
	}

	am[n] += loop;
	am[a - n - o] += loop;
	am[atom(0, 0)] += loop;
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

static void update_fissile(atom_map& am, std::mt19937& gen, atom a, long loop)
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

		am[r1] += 1;
		am[r2] += 1;
		am[atom(0, 1, true)] += n;

//		std::cout << "fissioning " << n << " into " << r1 << ", " << (a - r1) << ", and " << atom{0, 1, true} << " x " << n << std::endl;
	}
	
	am[atom(0, 0)] += loop;
}

void fuel_rod::update_neutrons(atom_map& am, std::mt19937& gen)
{
	// get all the neutrons
	std::set<long> neutrons;
	std::uniform_int_distribution<long> dist_n(0, mass - 1);
	long& n_count = all[atom(0, 1, true)];

	for(long i = 0; i < n_count; i++)
	{
		neutrons.insert(dist_n(gen));
	}

	auto it_n = neutrons.begin();
	long mass_at = 0;
	long mass_n = 0;

	for(auto& [a, c] : all)
	{
		mass_at += a.mass() * c;

		if(c <= 0 || a.mass() == 0 || a.protons() == 0)
		{
			continue;
		}

		long n = 0;

		while(it_n != neutrons.end() && mass_n < mass_at)
		{
			mass_n = *it_n;
			it_n++;
			n++;
		}

		if(n == 0)
		{
			continue;
		}

		if(n > c)
		{
			n = c;
		}
		
		atom a_n = a + atom{0, 1};
		a_n.excited = true;
		
		if(decay::is_fissile(a))
		{
			update_fissile(am, gen, a_n, n);
		}

		else
		{
			am[a_n] += n;
		}

		n_count -= n;
		c -= n;
	}
}

void fuel_rod::update_decays(atom_map& am, std::mt19937& gen, double secs)
{
	auto next = all.begin();

	for(auto it = all.begin(); it != all.end(); it = next)
	{
		auto& [a, c] = *it;
		next = it;
		next++;

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
			update_fissile(am, gen, a, loop);
		}

		else
		{
			update_decay(am, a, loop);
		}
		
		c -= loop;
	}
}

void fuel_rod::update(std::mt19937& gen, double secs)
{
	atom_map am;
	update_neutrons(am, gen);
	update_decays(am, gen, secs);

	for(auto [a, n] : am)
	{
		all[a] += n;
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

