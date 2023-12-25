#pragma once

#include <random>
#include <array>

#include "atom.hpp"

namespace sim::fuel::decay
{

atom get_mode(atom a);
bool is_fissile(atom a);

const std::array<atom, 3> TYPES{{
	{2, 2},
	{1, -1},
	{-1, 1}
}};

};

