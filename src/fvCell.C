// Author - Oliver Marx ojm40@cam.ac.uk

#include <iostream>
#include <cstddef>
#include <vector>
#include <array>

#include "fvCell.H"

fvCell::fvCell(){
	valC = {0.0, 0.0, 0.0};
}

fvCell::fvCell(const std::array<double,3>& inCell){
	valC = {inCell[0],inCell[1],inCell[2]};
	m_isConservative = false;
}

fvCell::fvCell(const std::array<double,3>& inCell, bool form){
	valC = {inCell[0],inCell[1],inCell[2]};
	m_isConservative = form;
}

// Overloaded + operator for cells, will not deal with discrepencies between
// forms of the two cells
fvCell fvCell::operator+(const fvCell& a) const{
	fvCell result;
	for(std::size_t i=0;i<valC.size();i++)
		result.valC[i] = valC[i] + a.valC[i];
	return result;
}

// Overloaded - operator for cells
fvCell fvCell::operator-(const fvCell& a) const{
	fvCell result;
	for(std::size_t i=0;i<valC.size();i++)
		result.valC[i] = valC[i] - a.valC[i];
	return result;
}

// Overloaded - operator for cells
bool fvCell::operator==(const fvCell& a) const{
	return (valC[0]==a.valC[0]
		&& valC[1]==a.valC[1]
		&& valC[2]==a.valC[2]);
}

fvCell operator*(const double s, const fvCell& a){
	fvCell result;
	for(std::size_t i=0;i<3;i++)
		result[i] = s * a[i];
	return result;
}

fvCell operator*(const fvCell a, const double& s){
	fvCell result;
	for(std::size_t i=0;i<3;i++)
		result[i] = s * a[i];
	return result;
}

void fvCell::operator=(const fvCell& a){
	valC = a.valC;
	m_isConservative = a.m_isConservative;
}

double& fvCell::operator[](std::size_t i){
	return valC[i];
}

double fvCell::operator[](std::size_t i) const{
	return valC[i];
}

fvCell fvCell::toPrim(const double& gamma) const{
	double rho  = valC[0], 
	       u    = valC[1]/valC[0];
	double e    = calc_e(gamma);
	double p    = rho * e * (gamma - 1.0);

	return fvCell({rho, u, p}, false);
}

fvCell fvCell::toCons(const double& gamma){
	double rho  = valC[0], 
	       u    = valC[1];
	double kE   = 0.5 * rho * u * u; 
	double e    = calc_e(gamma);
	double E    = rho * e + kE;

	return fvCell({rho, rho*u, E}, true);
}

double fvCell::calc_e(const double& gamma) const{
	if(m_isConservative)
		return (valC[2] - 0.5 * valC[1] * valC[1] 
				/ valC[0]) / valC[0];
	else 
		return valC[2] / (valC[0] * (gamma - 1.0));
}
