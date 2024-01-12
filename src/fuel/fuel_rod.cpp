
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

long fuel_rod::extract_energy()
{
	long e = energy;

	energy = 0;

	return e;
}

void fuel_rod::add_mass(int id, atom a, long c)
{
	long m = c * a.mass();
	section& s = all[id];

	s.cluster[a] += c;
	s.mass += m;
	mass += m;
}

void fuel_rod::add_mass(atom a, long c)
{
	add_mass(0, a, c);
}

long fuel_rod::add_mass(double amount, atom_map& src)
{
	return add_mass(0, amount, src);
}

long fuel_rod::add_mass(int id, double amount, atom_map& src)
{
	section& s = all[id];
	long count = 0;

	for(auto& [a, c] : src)
	{
		long n = c * amount;

		c -= n;
		count += n;
		s.mass += n;
		mass += n;

		s.cluster[a] += n;
	}

	return count;
}

long fuel_rod::get_mass() const
{
	return mass;
}

long fuel_rod::get_mass(int id) const
{
	auto it = all.find(id);

	if(it == all.end()) return 0;
	
	return it->second.mass;
}

long fuel_rod::get_mass(atom a) const
{
	long size = 0;

	for(auto& [id, s] : all)
	{
		size += s.mass;
	}

	return size;
}

long fuel_rod::get_mass(int id, atom a) const
{
	auto it_a = all.find(id);

	if(it_a == all.end()) return 0;
	
	auto it_s = it_a->second.cluster.find(a);

	if(it_s == it_a->second.cluster.end()) return 0;

	return it_s->second;
}

long fuel_rod::remove_mass(double amount, atom_map& dst)
{
	return remove_mass(0, amount, dst);
}

long fuel_rod::remove_mass(int id, double amount, atom_map& dst)
{
	section& s = all[id];
	long count = 0;

	for(auto& [a, c] : s.cluster)
	{
		long n = c * amount;

		c -= n;
		count += n;
		s.mass -= n;
		mass -= n;
		
		dst[a] += n;
	}

	return count;
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
	for(auto& [id, s] : all)
	{
		long& count_p = s.cluster[atom(1, -1)];
		long& count_e = s.cluster[atom(-1, 1)];
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

	section& section0 = all[0];
	long& n_count = section0.cluster[atom(0, 1)];

	for(long i = 0; i < n_count; i++)
	{
		neutrons.insert(dist_n(gen));
	}

	auto it_n = neutrons.begin();
	long mass_at = 0;
	long mass_n = 0;

	for(auto& [id, s] : all)
	{
		atom_map& am = ac[id];

		for(auto& [a, c] : s.cluster)
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

			s.mass -= n;
			section0.mass += n;

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

	for(auto& [id, s] : all)
	{
		update_decays(s.cluster, ac[id], gen, secs);
	}

	section& section0 = all[0];

	for(auto& [id, am] : ac)
	{
		section& s = all[id];
		energy += am[atom(0, 0)];

		for(auto [a, n] : am)
		{
			if(a.p == 0 && a.n == 1)
			{
				section0.cluster[a] += n;
				section0.mass += n;
				s.mass -= n;
			}

			else
			{
				s.cluster[a] += n;
			}
		}
	}
}

long fuel_rod::calculate_mass() const
{
	long mass = 0;

	for(auto& [id, s] : all)
	{
		for(auto [a, c] : s.cluster)
		{
			mass += a.mass() * c;
		}
	}

	return mass;
}

void fuel_rod::display_cluster(int id, int top) const
{
	std::multimap<long, atom, std::greater<int>> sorted;
	auto it_s = all.find(id);
	
	std::cout << "\nsample contents (cluster " << id << "):\n";

	if(it_s == all.end())
	{
		std::cout << "empty\n";
		return;
	}

	const atom_map& cluster = it_s->second.cluster;
	int at = 0;
	
	for(auto [a, c] : cluster)
	{
		if(c <= 0) continue;
		sorted.insert({c, a});
	}


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

void fuel_rod::display_cluster(int id) const
{
	display_cluster(id, INT_MAX);
}

void fuel_rod::display(int top) const
{
	int cluster_count = 0;

	for(auto& [id, s] : all)
	{
		if(s.cluster.size() > 0)
		{
			cluster_count++;
		}
	}

	if(cluster_count > 1)
	{
		std::cout << "\nshowing all clusters:\n";
	}

	for(auto& [id, s] : all)
	{
		if(s.cluster.size() > 0)
		{
			display_cluster(id, top);
		}
	}
}

void fuel_rod::display() const
{
	display(INT_MAX);
}

