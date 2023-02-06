// Author - Oliver Marx ojm40@cam.ac.uk

#include <iostream>
#include <cstddef>
#include <vector>
#include <array>

#ifndef FVCELL_HEADER
#define FVCELL_HEADER

class fvCell{
	public:
		fvCell(const std::array<double,3>&);
		fvCell(const std::array<double,3>&, bool);
		fvCell();
		fvCell  operator+(const fvCell&) const;
		fvCell  operator-(const fvCell&) const;
		//fvCell  operator*(const double) const;
		void    operator=(const fvCell&);
		double& operator[](std::size_t);
		double  operator[](std::size_t) const;
		bool  operator==(const fvCell&) const;
		fvCell  toPrim(const double&) const;
		fvCell  toCons(const double&);

	private:
		//Private member functions
		double calc_e(const double&) const;
		std::array<double,3> valC;

		// Private member variables
		bool m_isConservative;
};

fvCell operator*(double, const fvCell&);
fvCell operator*(const fvCell&, double);

#endif
