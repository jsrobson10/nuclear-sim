
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

	constexpr atom(int p, int n) : p(p), n(n) {
	}

	constexpr atom() : p(0), n(0) {
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

	constexpr bool friend operator==(const atom& a, const atom& b) {
		return a.p == b.p && a.n == b.n;
	}
	
	constexpr atom friend operator+(const atom& a, const atom& b) {
		return {a.p + b.p, a.n + b.n};
	}

	constexpr atom friend operator-(const atom& a, const atom& b) {
		return {a.p - b.p, a.n - b.n};
	}

	friend std::ostream& operator<<(std::ostream& o, const atom& b) {
		o << b.get_symbol() << "-" << b.mass();
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
		size_t h1 = std::hash<int>{}(a.protons());
		size_t h2 = std::hash<int>{}(a.neutrons());
		return h1 ^ h2;
	}
};

