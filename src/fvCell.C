// Author - Oliver Marx ojm40@cam.ac.uk

#include <algorithm>
#include <iostream>
#include <cstddef>
#include <vector>
#include <array>

#include "fvCell.H"

fvCell::fvCell()
{
	for(int i=0;i<nVars;i++)
	{
		valC[i] = 0.0;
	}
	m_isConservative = true;
}

fvCell::fvCell(const std::array<double,nVars>& inCell)
{
	for(int i=0;i<nVars;i++)
	{
		valC[i] = inCell[i];
	}
	m_isConservative = false;
}

fvCell::fvCell(const std::array<double,nVars>& inCell, bool form)
{
	for(int i=0;i<nVars;i++)
	{
		valC[i] = inCell[i];
	}
	m_isConservative = form;
}

// Overloaded + operator for cells, will not deal with discrepencies between
// forms of the two cells
fvCell fvCell::operator+(const fvCell& a) const
{
	fvCell result;

	for(int i=0;i<nVars;i++)
		result.valC[i] = valC[i] + a.valC[i];

	if(m_isConservative==a.m_isConservative)
		result.m_isConservative = m_isConservative;
	else
		std::cout << "ADDING CELLS IN DIFFERENT FORMS" << std::endl;

	return result;
}

// Overloaded - operator for cells
fvCell fvCell::operator-(const fvCell& a) const
{
	fvCell result;

	for(int i=0;i<nVars;i++)
		result.valC[i] = valC[i] - a.valC[i];

	if(m_isConservative==a.m_isConservative)
		result.m_isConservative = m_isConservative;
	else
		std::cout << "SUBTRACTING CELLS IN DIFFERENT FORMS" << std::endl;

	return result;
}

// Overloaded - operator for cells
bool fvCell::operator==(const fvCell& a) const
{
	return std::equal(valC.begin(),valC.end(),a.valC.begin());
}

void fvCell::operator=(const fvCell& a){
	valC = a.valC;
	m_isConservative = a.m_isConservative;
}

double& fvCell::operator[](int i){
	return valC[i];
}

double fvCell::operator[](int i) const{
	return valC[i];
}

fvCell fvCell::toPrim(const double& gamma) const{
	double rho  = valC[0], 
	       u    = valC[1]/valC[0],
	       v    = valC[2]/valC[0];
	double e    = calc_e(gamma);
	double p    = rho * e * (gamma - 1.0);

	return fvCell({rho, u, v, p}, false);
}

fvCell fvCell::toCons(const double& gamma){
	double rho  = valC[0], 
	       u    = valC[1],
	       v    = valC[2];
	double kE   = 0.5 * rho * u * u; 
	double e    = calc_e(gamma);
	double E    = rho * e + kE;

	return fvCell({rho, rho*u, rho*v, E}, true);
}

double fvCell::calc_e(const double& gamma) const
{
	if(m_isConservative)
	{
		return (valC[3] - 0.5 * (valC[1] * valC[1] + valC[2] * valC[2]) / valC[0])
			/ valC[0];
	}
	else 
		return valC[3] / (valC[0] * (gamma - 1.0));
}

void fvCell::displayCell() const
{
	std::cout << valC[0] << " " << valC[1] << " " << valC[2] << " " << valC[3] << std::endl; 
}

// The following are defined operators for scalar multiplication and are not
// member funtions of fvCell, therefore, variables like nVars need to be updated
// seperately

fvCell operator*(const double s, const fvCell& a)
{
	int nVars = 4;
	fvCell result = a;

	for(int i=0;i<nVars;i++)
		result[i] = s * a[i];

	return result;
}

fvCell operator*(const fvCell a, const double& s)
{
	int nVars = 4;
	fvCell result = a;

	for(int i=0;i<nVars;i++)
		result[i] = s * a[i];

	return result;
}
