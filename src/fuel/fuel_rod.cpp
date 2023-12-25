
#include "fuel_rod.hpp"
#include "half_life.hpp"
#include "decay.hpp"

#include <set>
#include <map>
#include <iostream>
#include <random>
#include <cmath>
#include <climits>

using namespace sim::fuel;

void fuel_rod::add_mass(int id, atom a, long c)
{
	all[id][a] += c;
	mass += c * a.mass();
}

void fuel_rod::add_mass(atom a, long c)
{
	add_mass(0, a, c);
}

static void update_decay(atom_map& am, atom a, long loop)
{
	atom n = decay::get_mode(a);
	
	if(a.excited && n == atom{-1, 1} && a.n > 2)
	{
		atom o {0, 1};
		am[o] += loop;
		am[a - n - o] += loop;
	}

	else
	{
		am[a - n] += loop;
	}

	am[n] += loop;
	am[atom(0, 0)] += loop;
}

static atom do_single_fission(std::mt19937& gen, atom a, int mass_1)
{
	std::uniform_real_distribution<double> dist;
	
	int p = a.p;
	int n = a.n;
	
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

	return {a.p - p, a.n - n};
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

		am[!r1] += 1;
		am[!r2] += 1;
		am[atom(0, 1)] += n;

//		std::cout << "fissioning " << n << " into " << r1 << ", " << (a - r1) << ", and " << atom{0, 1, true} << " x " << n << std::endl;
	}
	
	am[atom(0, 0)] += loop;
}

void fuel_rod::update_charges(std::map<int, atom_map>& ac)
{
	for(auto& [id, cluster] : all)
	{
		long& count_p = cluster[atom(1, -1)];
		long& count_e = cluster[atom(-1, 1)];
		long count = std::min(count_p, count_e);
		
		ac[id][atom(0, 0)] += count;

		count_p -= count;
		count_e -= count;
	}
}

void fuel_rod::update_neutrons(std::map<int, atom_map>& ac, std::mt19937& gen)
{
	// get all the neutrons
	std::set<long> neutrons;
	std::uniform_int_distribution<long> dist_n(0, mass - 1);
	long& n_count = all[0][atom(0, 1)];

	for(long i = 0; i < n_count; i++)
	{
		neutrons.insert(dist_n(gen));
	}

	auto it_n = neutrons.begin();
	long mass_at = 0;
	long mass_n = 0;

	for(auto& [id, cluster] : all)
	{
		atom_map& am = ac[id];

		for(auto& [a, c] : cluster)
		{
			mass_at += a.mass() * c;

			if(c <= 0 || a.mass() == 0 || a.p == 0)
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
			
			atom a_n = !(a + atom{0, 1});
			
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
}

void fuel_rod::update_decays(atom_map& cluster, atom_map& am, std::mt19937& gen, double secs)
{
	auto next = cluster.begin();

	for(auto it = next; it != cluster.end(); it = next)
	{
		auto& [a, c] = *it;
		next = it;
		next++;

		if(c <= 0)
		{
			cluster.erase(it);
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
	std::map<int, atom_map> ac;
	
	update_charges(ac);
	update_neutrons(ac, gen);

	for(auto& [id, cluster] : all)
	{
		update_decays(cluster, ac[id], gen, secs);
	}

	atom_map& cluster_0 = all[0];

	for(auto& [id, am] : ac)
	{
		atom_map& cluster = all[id];

		for(auto [a, n] : am)
		{
			if(a.p == 0 && a.n == 1)
			{
				cluster_0[a] += n;
			}

			else
			{
				cluster[a] += n;
			}
		}
	}
}

long fuel_rod::calculate_mass()
{
	long mass = 0;

	for(auto& [id, cluster] : all)
	{
		for(auto [a, c] : cluster)
		{
			mass += a.mass() * c;
		}
	}

	return mass;
}

void fuel_rod::display_cluster(int id, int top)
{
	std::multimap<long, atom, std::greater<int>> sorted;
	atom_map& cluster = all[id];
	int at = 0;
	
	for(auto [a, c] : cluster)
	{
		if(c <= 0) continue;
		sorted.insert({c, a});
	}

	std::cout << "\nsample contents (cluster " << id << "):\n";

	for(auto [c, a] : sorted)
	{
		if(++at > top) break;
		std::cout << "  " << a << ": " << c << std::endl;
	}

	if(at != sorted.size())
	{
		std::cout << "... " << (sorted.size() - at) << " hidden" << std::endl;
	}

	std::cout << std::endl;
}

void fuel_rod::display_cluster(int id)
{
	display_cluster(id, INT_MAX);
}

void fuel_rod::display(int top)
{
	int cluster_count = 0;

	for(auto& [id, cluster] : all)
	{
		if(cluster.size() > 0)
		{
			cluster_count++;
		}
	}

	if(cluster_count > 1)
	{
		std::cout << "\nshowing all clusters:\n";
	}

	for(auto& [id, cluster] : all)
	{
		if(cluster.size() > 0)
		{
			display_cluster(id, top);
		}
	}
}

void fuel_rod::display()
{
	display(INT_MAX);
}

