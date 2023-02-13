// Author - Oliver Marx ojm40@cam.ac.uk

// Standard packages
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstddef>
#include <iterator>
#include <string>
#include <vector>
#include <cmath>
#include <functional>

// Libs
#include <libconfig.h++>

// Headers
#include "fvSim.H"
#include "fvCell.H"

// Constructor, will take in settings file and test number to construct full sim
// with initial conditions and configs for certain tests
fvSim::fvSim(const char* configFileName, int testNum, std::string solver){
	m_solver = solver;
	gamma = 1.4;
	m_ghost = 2;

	libconfig::Config cfg;

	try{
		cfg.readFile(configFileName);
	}
	catch(const libconfig::FileIOException &fioex){
		std::cerr << "I/O error while reading file." << std::endl;
	}
	
	const libconfig::Setting& root = cfg.getRoot();

	try{
		const libconfig::Setting& tests = root["simulation"]["tests"];
		const libconfig::Setting& test = tests[testNum];
		if(!(test.lookupValue("nCells", nCells)
			&& test.lookupValue("x0", x0)
			&& test.lookupValue("x1", x1)
			&& test.lookupValue("t0", t0)
			&& test.lookupValue("t1", t1)
			&& test.lookupValue("CFL", CFL)
			&& test.lookupValue("gamma", gamma)))
			std::cout << "Settings for test #" << testNum+1 << " read in."
				<< std::endl;
		const libconfig::Setting& range = tests[testNum]["inits"];
		std::cout << "Solv: " << m_solver << " Test: " << testNum << " x0 " 
			<< x0 << " x1 " << x1 << " t0 " << t0 << " t1 " << t1 << " CFL " << CFL << std::endl;

		// Init storage arrays
		int problemSize = nCells + 2 * m_ghost;

		Q.resize(nCells);
		Q_new.resize(nCells);
		Q_bc.resize(problemSize);
		f_half.resize(nCells+1);

		xCentroids.resize(nCells);

		// Initialize constant run-time variables
		dx = (x1 - x0) / nCells; 

		// Create output file name and convert to char*
		std::string s = "/home/ojm40/Documents/MPhil_MHD/dat/output_" 
			+ m_solver + "_test" + std::to_string(testNum) + ".dat";
		const char* outputFileName = s.c_str();

		outputFile.open(outputFileName);

		// Initial conditions based on test case
		init(range);

	} catch(const libconfig::SettingNotFoundException &nfex){
		//Ignore
	}

	// Write initial data to file
	output();
}

// Destructor
fvSim::~fvSim(){
	outputFile.close();
}

void fvSim::run(){

	// Create a function pointer based on chosen flux function
	// https://stackoverflow.com/questions/7582546/using-generic-stdfunction-objects-with-member-functions-in-one-class
	std::function<const fvCell(const fvCell&,const fvCell&)> f; 
	std::function<std::array<fvCell,2>(const fvCell&,const fvCell&,const fvCell&)> g; 
	if(m_solver == "LF")
	{
		f = std::bind(&fvSim::LF_Flux, this, std::placeholders::_1, std::placeholders::_2);
		g = std::bind(&fvSim::constantDataReconstruction, this, std::placeholders::_1,
				std::placeholders::_2, std::placeholders::_3);
	}
	else if(m_solver == "FORCE")
	{
		f = std::bind(&fvSim::FORCE_Flux, this, std::placeholders::_1, std::placeholders::_2);
		g = std::bind(&fvSim::constantDataReconstruction, this, std::placeholders::_1,
				std::placeholders::_2, std::placeholders::_3);
	}
	else if(m_solver == "HLL")
	{
		f = std::bind(&fvSim::HLL_Flux, this, std::placeholders::_1, std::placeholders::_2);
		g = std::bind(&fvSim::constantDataReconstruction, this, std::placeholders::_1,
				std::placeholders::_2, std::placeholders::_3);
	}
	else if(m_solver == "SLIC")
	{
		f = std::bind(&fvSim::FORCE_Flux, this, std::placeholders::_1, std::placeholders::_2);
		g = std::bind(&fvSim::linearDataReconstruction, this, std::placeholders::_1,
				std::placeholders::_2, std::placeholders::_3);
	}
	else if(m_solver == "HLLC")
	{
		f = std::bind(&fvSim::HLLC_Flux, this, std::placeholders::_1, std::placeholders::_2);
		g = std::bind(&fvSim::constantDataReconstruction, this, std::placeholders::_1,
				std::placeholders::_2, std::placeholders::_3);
	}

	double t = t0;
	do
	{
		computeTimeStep();
		std::cout << "dt " << dt << std::endl;
		t = t + dt;

		// Copy the actual domain into the ghost cell domain
		for(int i=0;i<nCells;i++)
		{
			Q_bc[i+m_ghost] = Q[i];
		}

		// Do boundary conditions here 
		// Transmissive
		for(int i=0;i<m_ghost;i++)
		{
			Q_bc[i] = Q[0];
			Q_bc[nCells + m_ghost + i] = Q[nCells-1];
		}

		// Qi stores left and right states
		std::array<fvCell,2> Qi;

		std::array<fvCell,2> Qiplus1 = reconstructData(Q_bc[0], Q_bc[1], Q_bc[2], g);
				
		for(int i=1;i<nCells+m_ghost;i++)
		{
			Qi = Qiplus1;

			Qiplus1 = reconstructData(Q_bc[i], Q_bc[i+1], Q_bc[i+2], g);

			const fvCell QL = Qi[1] - (0.5 * dt / dx) * (F(Qi[1]) - F(Qi[0]));
			const fvCell QR = Qiplus1[0] - (0.5 * dt / dx) * (F(Qiplus1[1]) - F(Qiplus1[0]));

			f_half[i-1] = calculateFlux(QL, QR, f);
		}
		
		fullTimeStepUpdate();
		outputFile << "\n\n";
		outputFile << "t=" << t << "\n";
		output();

	} while (t < t1);
}

// This function will initialize the grid based on the test case chosen
void fvSim::init(const libconfig::Setting& ranges)
{
	int count = ranges.getLength();

	for(int i=0;i<count;i++)
	{
		const libconfig::Setting& range = ranges[i];
		double xi, xf, rho, u, v, p;

		if(!(range.lookupValue("xi", xi)
			&& range.lookupValue("xf", xf)
			&& range.lookupValue("rho", rho)
			&& range.lookupValue("u", u)
			&& range.lookupValue("v", v)
			&& range.lookupValue("p", p)))
			continue;

		for(int j=0;j<nCells;j++)
		{
			double x = x0 + (j + 0.5) * dx;
			double y = 0.0;

			xCentroids[j] = x;

			if(x <= xf && x >= xi)
				Q[j] = fvCell({rho, u, v, p}).toCons(gamma);
		}
	}
}

// Utility to output to a ofstream when needed
void fvSim::output()
{
	for(int i=0;i<nCells;i++)
	{
		fvCell Wi = Q[i].toPrim(gamma);

		outputFile << xCentroids[i] << " " 
			   << Wi[0] << " " 
			   << Wi[1] << " " 
			   << Wi[2] << " "
			   << Wi[3] << " "
			   << Wi.calc_e(gamma)
			   << std::endl;
	}
}

// compute a full time step update with flux values at cell boundaries
// This is not specific to a certain type of flux, and can therefore be defined
// using a centred scheme or RP based scheme
// Uses Q_i_n, Q_i+1_n, f_i+1/2_n and f_i-1/2_n to give Q_i_n+1
void fvSim::fullTimeStepUpdate()
{
	for(int i=0;i<nCells;i++)
	{
		Q_new[i] = Q_bc[i+m_ghost] - (dt / dx) * (f_half[i+1] - f_half[i]);
	}
	Q = Q_new;
}

// Compute the stable time step based on the grid size and max velocity
void fvSim::computeTimeStep()
{
	// Initialize temp Cell to store primitive form
	double SL, SR;
	double max = 0.0;

	for(int i=0;i<nCells-1;i++)
	{
		waveSpeedEstimates(Q[i], Q[i+1], SL, SR);
		SL = fabs(SL);
		SR = fabs(SR);
		max = std::max(max, std::max(SL,SR));;
	}

	dt = CFL * dx / max;
}

// The flux function for the PDE in conservation form. In this case it is the
// Burgers' Flux ---- To be replaced by Euler when ready
const fvCell fvSim::F(const fvCell& Qi)
{
	const fvCell Wi = Qi.toPrim(gamma);
	
	double f0 = Wi[0] * Wi[1];
	double f1 = Wi[0] * Wi[1] * Wi[1] + Wi[3];
	double f2 = Wi[0] * Wi[1] * Wi[2];
	double f3 = (Qi[3] + Wi[3]) * Wi[1];

	return fvCell({f0, f1, f2, f3}, true);
}

// Compute a half time step update of the data in Q_i_n giving Q_i+1/2_n+1/2
const fvCell fvSim::halfTimeStepUpdate(const fvCell& QL, const fvCell& QR){
	return 0.5 * (QR + QL) - (0.5 * dt / dx) * (F(QR) - F(QL));
}

// Lax-Friedrichs Flux calculation, takes left amd right states at cell boundary
// gives the flux at a cell boundary f_i+1/2_n
const fvCell fvSim::LF_Flux(const fvCell& QL, const fvCell& QR){
	return 0.5 * (dx / dt) * (QL - QR) + 0.5 * (F(QR) + F(QL));
}

// Richtmyer Flux calculation, takes left and right cell boundary values of
// half-time step updated data and simply calculates the flux function
const fvCell fvSim::Richt_Flux(const fvCell& Q_half){
	return F(Q_half);
}

// FORCE Flux is the average of a LF and Richt Flux. So, we find those at
// Q_i+1/2_n and then calculate tye flux.
const fvCell fvSim::FORCE_Flux(const fvCell& QL, const fvCell& QR){
	const fvCell LF = LF_Flux(QL, QR);
	const fvCell Q_half = halfTimeStepUpdate(QL, QR);
	const fvCell Richt = F(Q_half);
	return 0.5 * (LF + Richt);
}

const fvCell fvSim::HLL_Flux(const fvCell& QL, const fvCell& QR)
{
	double SL, SR;
	waveSpeedEstimates(QL, QR, SL, SR);

	const fvCell fL = F(QL);                                                        
	const fvCell fR = F(QR);                                                        
	const fvCell fHLL = (1.0 / (SR-SL)) * (SR*fL - SL*fR + SL * SR * (QR - QL));

       	if(SL >= 0)
	{
		return fL;
	}else if(SR <= 0)
	{
       	        return fR; 
	}else
	{
       	        return fHLL;
	}
       	                                                                          
}                    

const fvCell fvSim::HLLC_Flux(const fvCell& QL, const fvCell& QR)
{
	double SL, SR, Sstar;
	waveSpeedEstimates(QL, QR, SL, SR);

	const fvCell WL = QL.toPrim(gamma);
	const fvCell WR = QR.toPrim(gamma);

	const double& rhoL  = QL[0],
		      rhoR  = QR[0],
		      uL    = WL[1],
		      uR    = WR[1],
		      vL    = WL[2],
		      vR    = WR[2],
		      pL    = WL[3],
		      pR    = WR[3],
		      rhouL = QL[1],
		      rhouR = QR[1],
		      EL    = QL[3],
		      ER    = QR[3];

	Sstar = (pR - pL + rhouL * (SL - uL) - rhouR * (SR - uR)) 
		/ (rhoL * (SL - uL) - rhoR * (SR - uR));

	double prefixL = rhoL * (SL - uL) / (SL - Sstar);
	double prefixR = rhoR * (SR - uR) / (SR - Sstar);

	double ELstar = EL / rhoL + (Sstar - uL) * (Sstar + pL / (rhoL * (SL - uL)));
	double ERstar = ER / rhoR + (Sstar - uR) * (Sstar + pR / (rhoR * (SR - uR)));

	const fvCell QstarL = prefixL * fvCell({1.0, Sstar, vL, ELstar}, true);
	const fvCell QstarR = prefixR * fvCell({1.0, Sstar, vL, ERstar}, true);

	const fvCell fL = F(QL);                                                        
	const fvCell fR = F(QR);                                                        
	const fvCell fstarL = fL + SL * (QstarL - QL);
	const fvCell fstarR = fR + SR * (QstarR - QR);

       	if(SL >= 0)
	{
		return fL; 
	}else if (SL <= 0.0 && Sstar >= 0.0)
	{
       	        return fstarL;
	}else if (Sstar <= 0.0 && SR >= 0.0)
	{
		return fstarR;
	}else
	{
		return fR;
	}
}

std::array<fvCell,2> fvSim::constantDataReconstruction(const fvCell& Ql, const fvCell& Qi, const fvCell& QR)
{
	std::array<fvCell,2> result;

	result[0] = Qi;
	result[1] = Qi;
	
	return result;
}

std::array<fvCell,2> fvSim::linearDataReconstruction(const fvCell& QL, const fvCell& Qi, const fvCell& QR){
	fvCell dL = Qi - QL; 
	double dL_E = fabs(dL[3]) <= 1.0e-8 ? 1.0e-8 : dL[2];
	fvCell dR = QR - Qi; 
	double dR_E = fabs(dR[3]) <= 1.0e-8 ? 1.0e-8 : dR[2];

	double r = dL_E / dR_E;

	double xi = superbee(r);

	fvCell di = 0.5 * (dL + dR);

	std::array<fvCell,2> result;
	result[0] = Qi - 0.5 * xi * di;
	result[1] = Qi + 0.5 * xi * di;
	
	return result;
}

double fvSim::superbee(const double& r)
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

void fvSim::waveSpeedEstimates(const fvCell& QL, const fvCell& QR, double& SL, double& SR)
{
       fvCell WL = QL.toPrim(gamma);                                             
       fvCell WR = QR.toPrim(gamma);                                             

       double aL  = sqrt(gamma * WL[3] / WL[0]);                                
       double aR  = sqrt(gamma * WR[3] / WR[0]);                                

       const double& rhoL  = QL[0],
		     rhoR  = QR[0],
		     uL    = WL[1],
		     uR    = WR[1],
		     vL    = WL[2],
		     vR    = WR[2],
		     pL    = WL[3],
		     pR    = WR[3];

	double rhobar = 0.5 * (rhoL + rhoR);
	double abar   = 0.5 * (aL + aR);
	double ppvrs_x  = 0.5 * (pL + pR) - 0.5 * (uR - uL) * rhobar * abar;
	double ppvrs_y  = 0.5 * (pL + pR) - 0.5 * (vR - vL) * rhobar * abar;
	double pstar  = std::max(0.0, std::max(ppvrs_x, ppvrs_y));

	double qL, qR;

	//std::cout << "pL " << pL << " pstar " << pstar << " pR " << pR << std::endl; 

	if(pstar <= pL)
	{
		qL = 1.0;
	}
	else if(pstar > pL)
	{
		qL = sqrt(1.0 + ((gamma + 1.0) / (2.0 * gamma)) * (pstar / pL - 1.0));
	}

	if(pstar <= pR)
	{
		qR = 1.0;
	}
	else if(pstar > pR)
	{
		qR = sqrt(1.0 + ((gamma + 1.0) / (2.0 * gamma)) * (pstar / pR - 1.0));
	}

	SL = std::max(uL - aL * qL, vL - aL * qL);
	SR = std::max(uR + aR * qR, vR + aR * qR);
}

// Calculate flux is used to parse any chosen flux function in based on user
// choice or testing, to reduce repetition of code.
const fvCell fvSim::calculateFlux(const fvCell& QL, const fvCell& QR,
		std::function<const fvCell(const fvCell&,const fvCell&)> func){
	return func(QL, QR);
}

// Reconstruct Data is used to activate linear data reconstruction or not.
std::array<fvCell,2> fvSim::reconstructData(const fvCell& QL, const fvCell& Qi, const fvCell& QR,
		std::function<std::array<fvCell,2>(const fvCell&,const fvCell&, const fvCell&)> func){
	return func(QL, Qi, QR);
}
