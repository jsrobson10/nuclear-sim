
#include "fuel/fuel_assembly.hpp"

#include <cmath>
#include <sstream>
#include <random>

#include <unistd.h>
#include <curses.h>

int main()
{
	std::mt19937 gen(std::random_device{}());

	initscr();
	cbreak();
	noecho();
	keypad(stdscr, TRUE);
	nodelay(stdscr, TRUE);

	sim::fuel::fuel_assembly fa(25, 100000000);
	
	fa.add_mass(sim::fuel::atom::from_symbol("U", 235), 500000);
	fa.add_mass({0, 1}, 1000);

	double r = 0;
	double temp_out;
	double temp_last = 0;
	double temp = 0;
	bool mode_a = false;

	fa.set_reactivity(r);

	for(;;)
	{
		move(0, 0);

		fa.update(gen, 1);
		temp_out = fa.extract_temperature(25, 0.01);
		temp = fa.get_temperature();

		fa.set_reactivity(r);

		double m = 0;
		std::stringstream ss;

		ss << fa;
		ss << "power output: " << temp_out << std::endl;
		ss << "auto mode: " << mode_a << std::endl;

		addstr(ss.str().c_str());
		refresh();

		if(mode_a)
		{
			m = (temp - temp_last > 0) ? -1 : 1;
		}

		temp_last = temp;

		int c = getch();

		switch(c)
		{
		case KEY_DOWN:
			m = -1;
			break;
		case KEY_UP:
			m = 1;
			break;
		case 'a':
			mode_a = !mode_a;
			break;
		}

		r += m * 0.01;

		usleep(10000);
	}

	endwin();
}

