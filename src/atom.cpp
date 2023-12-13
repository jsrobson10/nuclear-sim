
#include "atom.hpp"

#include <cstring>

using namespace sim;

static const int SYMBOLS_LEN = 118;
static const char* SYMBOLS[SYMBOLS_LEN]
{
	"H",
	"He",
	"Li",
	"Be",
	"B",
	"C",
	"N",
	"O",
	"F",
	"Ne",
	"Na",
	"Mg",
	"Al",
	"Si",
	"P",
	"S",
	"Cl",
	"Ar",
	"K",
	"Ca",
	"Sc",
	"Ti",
	"V",
	"Cr",
	"Mn",
	"Fe",
	"Co",
	"Ni",
	"Cu",
	"Zn",
	"Ga",
	"Ge",
	"As",
	"Se",
	"Br",
	"Kr",
	"Rb",
	"Sr",
	"Y",
	"Zr",
	"Nb",
	"Mo",
	"Tc",
	"Ru",
	"Rh",
	"Pd",
	"Ag",
	"Cd",
	"In",
	"Sn",
	"Sb",
	"Te",
	"I",
	"Xe",
	"Cs",
	"Ba",
	"La",
	"Ce",
	"Pr",
	"Nd",
	"Pm",
	"Sm",
	"Eu",
	"Gd",
	"Tb",
	"Dy",
	"Ho",
	"Er",
	"Tm",
	"Yb",
	"Lu",
	"Hf",
	"Ta",
	"W",
	"Re",
	"Os",
	"Ir",
	"Pt",
	"Au",
	"Hg",
	"Tl",
	"Pb",
	"Bi",
	"Po",
	"At",
	"Rn",
	"Fr",
	"Ra",
	"Ac",
	"Th",
	"Pa",
	"U",
	"Np",
	"Pu",
	"Am",
	"Cm",
	"Bk",
	"Cf",
	"Es",
	"Fm",
	"Md",
	"No",
	"Lr",
	"Rf",
	"Db",
	"Sg",
	"Bh",
	"Hs",
	"Mt",
	"Ds",
	"Rg",
	"Cn",
	"Nh",
	"Fl",
	"Mc",
	"Lv",
	"Ts",
	"Og"
};

constexpr int* fill_lookup_page(char a)
{
	int* table = new int[27];

	for(int i = 0; i < SYMBOLS_LEN; i++)
	{
		const char* sym = SYMBOLS[i];

		if(sym[0] != a)
		{
			continue;
		}

		int loc = 0;

		if(sym[1] != '\0')
		{
			loc = (sym[1] & 255) - 'a' + 1;
		}
		
		table[loc] = i + 1;
	}

	return table;
}

constexpr int** create_lookup()
{
	int** lookup = new int*[26];

	for(int i = 0; i < 26; i++)
	{
		lookup[i] = fill_lookup_page('A' + i);
	}

	return lookup;
}

static const int *const *const SYMBOLS_R = create_lookup();

const char* atom::get_symbol() const {
	if(p == 0 && n == 0) return "γ";
	else if(p + n == 0) return "β";
	else if(p == 1 && n == 0) return "P";
	else if(p == 0 && n == 1) return "N";
	else if(p > 0 && p < 119) return SYMBOLS[p - 1];
	else return "Unknown";
}


atom atom::from_symbol(const char* sym, int mass)
{
	int len = strlen(sym);

	if((len != 1 && len != 2) || sym[0] < 'A' || sym[0] > 'Z')
	{
		return {0, 0};
	}

	int loc1 = (sym[0] - 'A') & 255;
	int loc2 = 0;

	if(sym[1] != '\0')
	{
		if(sym[1] < 'a' || sym[1] > 'z')
		{
			return {0, 0};
		}

		loc2 = (sym[1] - 'a' + 1) & 255;
	}

	return from_mass(SYMBOLS_R[loc1][loc2], mass);
}

