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
			&& test.lookupValue("CFL", CFL)))
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
		double xi, xf, Qin;

		if(!(range.lookupValue("xi", xi)
			&& range.lookupValue("xf", xf)
			&& range.lookupValue("Q", Qin)))
			continue;

		for(std::size_t j=0;j<Q.size();j++){
			double x = x0 + (j-1) * dx;
			if(x < xf && x >= xi)
				Q[j] = Qin;
		}
	}
}

// Utility to output to a ofstream when needed
void fvSim::output(){
	for(std::size_t i=1;i<nCells+1;i++){
		outputFile << (x0 + (i-1)*dx) << " " << Q[i] << std::endl;
	}
}

// Compute the stable time step based on the grid size and max velocity
void fvSim::computeTimeStep(){
	double max = Q[1];
	for(std::size_t i=1;i<nCells+1;i++){
		if(fabs(Q[i]) >= max) max = fabs(Q[i]);
	}
	dt = CFL * dx / max;
}

// The flux function for the PDE in conservation form. In this case it is the
// Burgers' Flux ---- To be replaced by Euler when ready
double fvSim::F(double Qi){
	return 0.5 * Qi * Qi;
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
double fvSim::halfTimeStepUpdate(double QL, double QR){
	return 0.5 * (QR + QL) - 0.5 * (dt / dx) * (F(QR) - F(QR));
}

// Lax-Friedrichs Flux calculation, takes left amd right states at cell boundary
// gives the flux at a cell boundary f_i+1/2_n
double fvSim::LF_Flux(double QL, double QR){
	return 0.5 * (dx / dt) * (QL - QR) + 0.5 * (F(QR) - F(QL));
}

// Richtmyer Flux calculation, takes left and right cell boundary values of
// half-time step updated data and simply calculates the flux function
double fvSim::Richt_Flux(double Q_half){
	return F(Q_half);
}

// FORCE Flux is the average of a LF and Richt Flux. So, we find those at
// Q_i+1/2_n and then calculate tye flux.
double fvSim::FORCE_Flux(double QL, double QR){
	double LF = LF_Flux(QL, QR);
	double Q_half = halfTimeStepUpdate(QL, QR);
	double Richt = Richt_Flux(Q_half);
	return 0.5 * (LF + Richt);
}

// Calculate flux is used to parse any chosen flux function in based on user
// choice or testing, to reduce repetition of code.
double fvSim::calculateFlux(double QL, double QR, std::function<double(double,double)> func){
	return func(QL, QR);
}

void fvSim::run(){

	// Create a function pointer based on chosen flux function
	// https://stackoverflow.com/questions/7582546/using-generic-stdfunction-objects-with-member-functions-in-one-class
	std::function<double(double,double)> f; 
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
		
		for(std::size_t i=0;i<nCells+1;i++)
			f_iplushalf_n[i] = calculateFlux(Q[i], Q[i+1], f );
		
		fullTimeStepUpdate();
		output();
		outputFile << "\n\n";
	} while (t < t1);
}
