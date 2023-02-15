// AUTHOR - Oliver Marx ojm40@cam.ac.uk

#include <array>
#include <cmath>

#include "fvCell.H"
#include "DataReconstruction.H"

std::array<fvCell,2> constantDataReconstruction(const fvCell& Ql, const fvCell& Qi, const fvCell& QR)
{
	std::array<fvCell,2> result;

	result[0] = Qi;
	result[1] = Qi;
	
	return result;
}

double superbee(const double& r)
{
	if(r<=0)
		return 0.0;
	else if(r>0 && r<=0.5)
		return 2 * r;
	else if(r>0.5 && r<=1)
		return 1.0;
	else{
		double xiR = 2.0 / (1+r);
		return std::min(r, std::min(2.0, xiR));
	}
}

double minbee(const double& r)
{
	if(r<=0)
		return 0.0;
	else if(r>0 && r<=1)
		return r;
	else{
		double xiR = 2.0 / (1+r);
		return std::min(r, std::min(1.0, xiR));
	}
}


std::array<fvCell,2> linearDataReconstruction(const fvCell& QL, const fvCell& Qi, const fvCell& QR){
	fvCell dL = Qi - QL; 
	double dL_E = (fabs(dL[3]) <= 1.0e-8) ? 1.0e-8 : dL[3];
	fvCell dR = QR - Qi; 
	double dR_E = (fabs(dR[3]) <= 1.0e-8) ? 1.0e-8 : dR[3];

	double r = dL_E / dR_E;

	double xi = minbee(r);

	fvCell di = 0.5 * (dL + dR);

	std::array<fvCell,2> result;
	result[0] = Qi - 0.5 * xi * di;
	result[1] = Qi + 0.5 * xi * di;
	
	return result;
}
