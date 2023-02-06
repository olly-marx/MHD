// Author - Oliver Marx ojm40@cam.ac.uk

// Standard packages
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
		Q.resize(nCells+2);
		Q_i_nplus1.resize(nCells+2);
		f_iplushalf_n.resize(nCells+1);

		// Initialize constant run-time variables
		dx = (x1 - x0) / nCells; 
		std::string s = "/home/ojm40/Documents/MPhil_MHD/dat/output_" 
			+ m_solver + "_test" + std::to_string(testNum) + ".dat";
		const char* outputFileName = s.c_str();
		outputFile.open(outputFileName);

		gamma = 1.4;

		// Initial conditions based on test case
		init(range);

	} catch(const libconfig::SettingNotFoundException &nfex){
		//Ignore
	}

	// Write initial data to file
	output();
	outputFile << "\n\n";
}

// Destructor
fvSim::~fvSim(){
	outputFile.close();
}

// This function will initialize the grid based on the test case chosen
void fvSim::init(const libconfig::Setting& ranges){
	int count = ranges.getLength();

	for(int i=0;i<count;i++){
		const libconfig::Setting& range = ranges[i];
		double xi, xf, rho, u, p;

		if(!(range.lookupValue("xi", xi)
			&& range.lookupValue("xf", xf)
			&& range.lookupValue("rho", rho)
			&& range.lookupValue("u", u)
			&& range.lookupValue("p", p)))
			continue;

		for(std::size_t j=0;j<Q.size();j++){
			double x = x0 + (j-1) * dx;
			if(x <= xf && x >= xi)
				Q[j] = fvCell({rho, u, p}).toCons(gamma);
		}
	}
}

// Utility to output to a ofstream when needed
void fvSim::output(){
	for(std::size_t i=1;i<nCells+1;i++){
		fvCell Wi = Q[i].toPrim(gamma);
		outputFile << (x0 + (i-1)*dx) << " " 
			   << Wi[0] << " " 
			   << Wi[1] << " " 
			   << Wi[2] << std::endl;
	}
}

// Compute the stable time step based on the grid size and max velocity
void fvSim::computeTimeStep(){
	fvCell Wi = Q[1].toPrim(gamma);
	double cs  = sqrt(gamma * Wi[2] / Wi[0]);
	double max = fabs(Wi[1]) + cs;
	for(std::size_t i=2;i<nCells+1;i++){
		Wi = Q[i].toPrim(gamma);
		cs  = sqrt(gamma * Wi[2] / Wi[0]);
		double s = fabs(Wi[1]) + cs;
		if(s >= max) max = s;
	}
	dt = CFL * dx / max;
}

// The flux function for the PDE in conservation form. In this case it is the
// Burgers' Flux ---- To be replaced by Euler when ready
const fvCell fvSim::F(const fvCell& Qi){
	const fvCell Wi = Qi.toPrim(gamma);
	fvCell result;
	result[0] = Qi[1];
	result[1] = Qi[1] * Wi[1] + Wi[2];
	result[2] = (Qi[2] + Wi[2]) * Wi[1];
	return result;
}

// compute a full time step update with flux values at cell boundaries
// This is not specific to a certain type of flux, and can therefore be defined
// using a centred scheme or RP based scheme
// Uses Q_i_n, Q_i+1_n, f_i+1/2_n and f_i-1/2_n to give Q_i_n+1
void fvSim::fullTimeStepUpdate(){
	for(std::size_t i=1;i<nCells+1;i++)
		Q_i_nplus1[i] = Q[i] - (dt / dx) * (f_iplushalf_n[i] - f_iplushalf_n[i-1]);
	Q = Q_i_nplus1;
}

// Compute a half time step update of the data in Q_i_n giving Q_i+1/2_n+1/2
const fvCell fvSim::halfTimeStepUpdate(const fvCell& QL, const fvCell& QR){
	return 0.5 * (QR + QL) - 0.5 * (dt / dx) * (F(QR) - F(QL));
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
	const fvCell Richt = Richt_Flux(Q_half);
	return 0.5 * (LF + Richt);
}

std::array<fvCell,2> fvSim::DataReconstruction(const fvCell& QL, const fvCell& Qi, const fvCell& QR){
	fvCell dL = Qi - QL; 
	double dL_E = fabs(dL[2]) <= 1.0e-8 ? 1.0e-8 : dL[2];
	fvCell dR = QR - Qi; 
	double dR_E = fabs(dR[2]) <= 1.0e-8 ? 1.0e-8 : dR[2];

	double r = dL_E / dR_E;

	double xi = superbee(r);

	fvCell di = 0.5 * (dL + dR);

	std::array<fvCell,2> result;
	result[0] = Qi - 0.5 * xi * di;
	result[1] = Qi + 0.5 * xi * di;
	
	return result;
}

double fvSim::superbee(const double& r){
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

// Calculate flux is used to parse any chosen flux function in based on user
// choice or testing, to reduce repetition of code.
const fvCell fvSim::calculateFlux(const fvCell& QL, const fvCell& QR,
		std::function<const fvCell(const fvCell&,const fvCell&)> func){
	return func(QL, QR);
}

void fvSim::run(){

	// Create a function pointer based on chosen flux function
	// https://stackoverflow.com/questions/7582546/using-generic-stdfunction-objects-with-member-functions-in-one-class
	std::function<const fvCell(const fvCell&,const fvCell&)> f; 
	if(m_solver == "LF")
		f = std::bind(&fvSim::LF_Flux, this, std::placeholders::_1, std::placeholders::_2);
	else if(m_solver == "FORCE")
		f = std::bind(&fvSim::FORCE_Flux, this, std::placeholders::_1, std::placeholders::_2);

	double t = t0;
	do{
		computeTimeStep();
		//std::cout << "dt " << dt << std::endl;
		t = t + dt;

		// Do boundary conditions here 
		// Transmissive
		Q[0] = Q[1];
		Q[nCells+1] = Q[nCells];

		std::array<fvCell,2> Qi;
		Qi[0] = Q[0];
		Qi[1] = Q[1];
		std::array<fvCell,2> Qiplus1;
				
		for(std::size_t i=1;i<nCells+2;i++){
			// Reconstruct data at cell i+1, and then set data at i
			// equal to data at i+1 for next iteration
			Qiplus1 = DataReconstruction(Q[i-1], Q[i], Q[i+1]);

			// We set the right value of the cell interface problem
			// equal to the left reconstrution of the right hand
			// cell
			// The left value of the interface problem is the right
			// reconstruction of the left cell
			fvCell QL = Qi[1] - 0.5 * (dt / dx) * (F(Qi[1]) - F(Qi[0]));
			fvCell QR = Qiplus1[0] - 0.5 * (dt / dx) * (F(Qiplus1[1]) - F(Qiplus1[0]));

			f_iplushalf_n[i-1] = calculateFlux(QL, QR, f );

			Qi = Qiplus1;
		}
		
		fullTimeStepUpdate();
		//outputFile << "t=" << t << "\n";
		output();
		outputFile << "\n\n";
	} while (t < t1);
}
