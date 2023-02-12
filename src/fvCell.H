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
		double& operator[](int);
		double  operator[](int) const;
		bool    operator==(const fvCell&) const;
		fvCell  toPrim(const double&) const;
		fvCell  toCons(const double&);
		double calc_e(const double&) const;
		void displayCell() const;
	private:
		// Private member variables
		std::array<double,3> valC;
		bool m_isConservative;
};

fvCell operator*(double, const fvCell&);
fvCell operator*(const fvCell&, double);

#endif
