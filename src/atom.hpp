
#pragma once

#include <ostream>
#include <unordered_map>
#include <cstdint>

namespace sim
{
		
class atom
{
private:
	
	short p, n;

public:
	
	bool excited;

	constexpr atom(int p, int n, bool e) : p(p), n(n), excited(e) {
	}

	constexpr atom(int p, int n) : p(p), n(n), excited(false) {
	}

	constexpr atom() : p(0), n(0), excited(false) {
	}

	constexpr int mass() const {
		return p + n;
	}

	constexpr int protons() const {
		return p;
	}

	constexpr int neutrons() const {
		return n;
	}

	static constexpr atom from_mass(int p, int m) {
		return atom(p, m - p);
	}

	static atom from_symbol(const char* sym, int mass);

	constexpr bool friend operator==(const atom& a, const atom& b) {
		return a.p == b.p && a.n == b.n && a.excited == b.excited;
	}
	
	constexpr atom friend operator+(const atom& a, const atom& b) {
		return {a.p + b.p, a.n + b.n};
	}

	constexpr atom friend operator-(const atom& a, const atom& b) {
		return {a.p - b.p, a.n - b.n};
	}

	friend std::ostream& operator<<(std::ostream& o, const atom& b)
	{
		if(b.excited)
		{
			o << "*";
		}

		o << b.get_symbol();

		if(b.mass() == 0)
		{
			if(b.protons() > 0) o << "+";
			else if(b.protons() < 0) o << "-";
		}

		else if(b.mass() > 1)
		{
			o << "-" << b.mass();
		}

		return o;
	}

	constexpr atom& operator+=(const atom& b) {
		p += b.p;
		n += b.n;
		return *this;
	}
	
	constexpr atom& operator-=(const atom& b) {
		p -= b.p;
		n -= b.n;
		return *this;
	}

	constexpr atom& operator=(const atom& b) {
		p = b.p;
		n = b.n;
		excited = b.excited;
		return *this;
	}
	
	const char* get_symbol() const;

	static const int U = 92;
	static const int Pu = 94;
};

struct atoms : atom
{
	long c;

	atoms(int p, int n, long c) : atom(p, n), c(c) {
	}
	
	atoms(atom a, long c) : atom(a), c(c) {
	}

	atoms() : atom(), c(0) {
	}
};

}

template<> struct std::hash<sim::atom>
{
	size_t operator() (sim::atom const& a) const noexcept
	{
		size_t h1 = std::hash<size_t>{}(a.protons());
		size_t h2 = std::hash<size_t>{}(h1 ^ a.neutrons());
		return h2 ^ a.excited;
	}
};

