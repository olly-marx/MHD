// Author - Oliver Marx ojm40@cam.ac.uk

// Standard packages
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstddef>
#include <iterator>
#include <ostream>
#include <string>
#include <vector>
#include <cmath>
#include <functional>

// Libs
#include <libconfig.h++>

// Headers
#include "fvSim.H"
#include "fvCell.H"
#include "DataReconstruction.H"
#include "FluxFunctions.H"

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
	catch (libconfig::ParseException &e){
        /*inform user about the parse exception*/
        std::cerr << "Parse error at " << e.getFile() << ":" << e.getLine()
                  << " - " << e.getError() << std::endl;
	}
	
	const libconfig::Setting& root = cfg.getRoot();

	try{
		const libconfig::Setting& tests = root["simulation"]["tests"];
		const libconfig::Setting& test = tests[testNum];

		if(!(test.lookupValue("nCellsX", nCellsX)
			&& test.lookupValue("nCellsY", nCellsY)
			&& test.lookupValue("x0", x0)
			&& test.lookupValue("x1", x1)
			&& test.lookupValue("y0", y0)
			&& test.lookupValue("y1", y1)
			&& test.lookupValue("t0", t0)
			&& test.lookupValue("t1", t1)
			&& test.lookupValue("CFL", CFL)
			&& test.lookupValue("test", m_test)
			&& test.lookupValue("gamma", gamma)))
			std::cout << "Settings for test #" << testNum+1 << " read in."
				<< std::endl;


		const libconfig::Setting& range = tests[testNum]["inits"];
		std::cout << "Solv: " << m_solver << " Test: " << (testNum+1) 
			<< " nCellsX " << nCellsX 
			<< " nCellsY " << nCellsY 
			<< " x0 " << x0 
			<< " x1 " << x1 
			<< " y0 " << y0 
			<< " y1 " << y1 
			<< " t0 " << t0 
			<< " t1 " << t1 
			<< " CFL " << CFL 
			<< std::endl;

		// Create output file name and convert to char*
		std::string s = "/home/ojm40/Documents/MPhil_MHD/dat/output_" 
			+ m_solver + "_test" + std::to_string(testNum+1) + "_1D.dat";
		const char* outputFileName = s.c_str();

		outputFile.open(outputFileName);

		// Initial conditions based on test case
		init(range);

	} catch(const libconfig::SettingNotFoundException &nfex){
		//Ignore
	}

	// Write initial data to file
	output();
	std::cout << "Running:  #Beg" << std::flush;
}

// Destructor
fvSim::~fvSim(){
	outputFile.close();
}

void fvSim::run(){

	// Create a function pointer based on chosen flux function
	// https://stackoverflow.com/questions/7582546/using-generic-stdfunction-objects-with-member-functions-in-one-class
	std::function<const fvCell(const fvCell&,const fvCell&, const double&, const double&, const double&, bool)> f; 
	std::function< std::array<fvCell,2> (const fvCell&,const fvCell&,const fvCell&) > g; 
	
	if(m_solver == "LF")
	{
		f = LF_Flux;
		g = constantDataReconstruction;
	}
	else if(m_solver == "FORCE")
	{
		f = FORCE_Flux;
		g = constantDataReconstruction;
	}
	else if(m_solver == "SLIC")
	{
		f = FORCE_Flux;
		g = linearDataReconstruction;
	}
	else if(m_solver == "HLLC")
	{
		f = HLLC_Flux;
		g = linearDataReconstruction;
	}
	else if(m_solver == "HLL")
	{
		f = HLL_Flux;
		g = linearDataReconstruction;
	}

	double printNumber = 1.0;
	bool printedAtTime = false;
	double nPrints = 10.0;
	double fracGap = 1.0/nPrints;
	t = t0;

	double CFLTemp = CFL;
	CFL = 0.2;
	int slowStart = 0;
	do
	{
		std::cout << "+" << std::flush;
		if(slowStart++ >= 5)
			CFL = CFLTemp;

		computeTimeStep();
		t = t + dt;
		std::cout << "dt " << dt << std::endl;

		// Copy the actual domain into the ghost cell domain
		for(int i=0;i<nCellsX;i++)
		{
			for(int j=0;j<nCellsY;j++)
			{
				Q_bc[i+m_ghost][j+m_ghost] = Q[i][j];
			}
		}

		// Do boundary conditions here 
		// Transmissive
		//
		// LENGTHS
		for(int i=0;i<m_ghost;i++)
		{
			for(int j=0;j<nCellsY;j++)
			{
				Q_bc[i][j+m_ghost] = Q[0][j];
				Q_bc[i + m_ghost + nCellsX][j + m_ghost] = Q[nCellsX-1][j];
			}
		}

		// Qi stores left and right states
		// Initially just for the x direction reconstruction and flux
		// calculations
		std::array<fvCell,2> Qi;

		bool x_dir = true;

		for(int j=0;j<nCellsY;j++)
		{
			std::array<fvCell,2> Qiplus1 = reconstructData(Q_bc[0][j+m_ghost],
								       Q_bc[1][j+m_ghost], 
								       Q_bc[2][j+m_ghost],
								       g );
				
			for(int i=1;i<nCellsX+m_ghost;i++)
			{

				Qi = Qiplus1;

				Qiplus1 = reconstructData(Q_bc[i][j+m_ghost],
						          Q_bc[i+1][j+m_ghost],
							  Q_bc[i+2][j+m_ghost], 
							  g);

				const fvCell dFL = (F(Qi[1], gamma, x_dir) - F(Qi[0], gamma, x_dir));
				const fvCell dFR = (F(Qiplus1[1], gamma, x_dir) - F(Qiplus1[0], gamma, x_dir));

				const fvCell QL = Qi[1] - (0.5 * dt / dx) * dFL;
				const fvCell QR = Qiplus1[0] - (0.5 * dt / dx) * dFR;

				f_half[i-1][j] = calculateFlux(QL, QR, dx, dt, gamma, x_dir, f);

			}
		}
		
		fullTimeStepUpdate(x_dir);

		//x_dir = false;

		//// Copy the actual domain into the ghost cell domain
		//for(int i=0;i<nCells;i++)
		//{
		//	for(int j=0;j<nCells;j++)
		//	{
		//		Q_bc[i+m_ghost][j+m_ghost] = Q[i][j];
		//	}
		//}

		//// BREADTHS
		//for(int i=0;i<nCells;i++)
		//{
		//	for(int j=0;j<m_ghost;j++)
		//	{
		//		Q_bc[i+m_ghost][j] = Q[i][0];
		//		Q_bc[i + m_ghost][j + m_ghost + nCells] = Q[i][nCells-1];
		//	}
		//}

		//for(int i=0;i<nCells;i++)
		//{
		//	std::array<fvCell,2> Qiplus1 = reconstructData(Q_bc[i+m_ghost][0],
		//						       Q_bc[i+m_ghost][1], 
		//						       Q_bc[i+m_ghost][2],
		//						       g );
		//		
		//	for(int j=1;j<nCells+m_ghost;j++)
		//	{

		//		Qi = Qiplus1;

		//		Qiplus1 = reconstructData(Q_bc[i+m_ghost][j],
		//				          Q_bc[i+m_ghost][j+1],
		//					  Q_bc[i+m_ghost][j+2], 
		//					  g);

		//		const fvCell dFL = (F(Qi[1], gamma, x_dir) - F(Qi[0], gamma, x_dir));
		//		const fvCell dFR = (F(Qiplus1[1], gamma, x_dir) - F(Qiplus1[0], gamma, x_dir));

		//		const fvCell QL = Qi[1] - (0.5 * dt / dx) * dFL;
		//		const fvCell QR = Qiplus1[0] - (0.5 * dt / dx) * dFR;

		//		f_half[i][j-1] = calculateFlux(QL, QR, dx, dt, gamma, x_dir, f);

		//	}
		//}
		//
		//fullTimeStepUpdate(x_dir);

		double tFrac = t/t1;
		double scale = pow(10.0, ceil(log10(fabs(tFrac))) + 2);
		tFrac = round(tFrac * scale) / scale;

		if(tFrac >= printNumber * fracGap && !printedAtTime)
		{
			outputFile << "\n\n";
			outputFile << "t=" << t << "\n";
			output();

			printedAtTime = true;
			std::cout << "#" << printNumber << std::endl;
		}
		else if(printedAtTime && tFrac > printNumber * fracGap)
		{
			printNumber+=1.0;
			printedAtTime = false;
		}

	} while (t < t1);

	outputFile << "\n\n";
	outputFile << "t=" << t << "\n";
	output();
	
	std::cout << "#Fin" << std::flush;
	std::cout << std::endl;
}

// compute a full time step update with flux values at cell boundaries
// This is not specific to a certain type of flux, and can therefore be defined
// using a centred scheme or RP based scheme
// Uses Q_i_n, Q_i+1_n, f_i+1/2_n and f_i-1/2_n to give Q_i_n+1
void fvSim::fullTimeStepUpdate(bool x_dir)
{
	int x = (x_dir==true);
	int y = (x_dir==false);

	for(int j=0;j<nCellsY;j++)
	{
		for(int i=0;i<nCellsX;i++)
		{
			Q_new[i][j] = Q_bc[i+m_ghost][j+m_ghost] - (dt / dx) * (f_half[i+x][j+y] - f_half[i][j]);
		}
	}
	Q = Q_new;
}

// Compute the stable time step based on the grid size and max velocity
void fvSim::computeTimeStep()
{
	// Initialize temp Cell to store primitive form
	double SL, SR;
	double max = 0.0;

	for(int j=0;j<nCellsY;j++)
	{
		for(int i=0;i<nCellsX-1;i++)
		{
			waveSpeedEstimatesMHD(Q[i][j], Q[i+1][j], SL, SR, gamma, true);
			SL = fabs(SL);
			SR = fabs(SR);
			max = std::max(max, std::max(SL,SR));;
		}
	}

	//for(int i=0;i<nCellsX;i++)
	//{
	//	for(int j=0;j<nCellsY-1;j++)
	//	{
	//		waveSpeedEstimatesMHD(Q[i][j], Q[i][j+1], SL, SR, gamma, false);
	//		SL = fabs(SL);
	//		SR = fabs(SR);
	//		max = std::max(max, std::max(SL,SR));;
	//	}
	//}

	dt = std::min(CFL * dx / max, t1 - t);
}


// This function will initialize the grid based on the test case chosen
void fvSim::init(const libconfig::Setting& ranges)
{
	// Init storage arrays
	int problemSizeX = nCellsX + 2 * m_ghost;
	int problemSizeY = nCellsY + 2 * m_ghost;

	Q.resize(nCellsX);
	Q_new.resize(nCellsX);
	Q_bc.resize(problemSizeX);
	f_half.resize(nCellsX+1);

	for(int i=0;i<problemSizeX;i++)
	{
		Q_bc[i].resize(problemSizeY);
		if(i<nCellsX)
		{
			Q[i].resize(nCellsY);
			Q_new[i].resize(nCellsY);
		}
		if(i<nCellsX+1)
		{
			f_half[i].resize(nCellsY+1);
		}
	}

	xCentroids.resize(nCellsX);
	yCentroids.resize(nCellsY);

	// Initialize constant run-time variables
	dx = (x1 - x0) / nCellsX; 
	dy = (y1 - y0) / nCellsY; 
	std::cout << dx << " " << dy << std::endl;

	int count = ranges.getLength();

	std::vector<std::array<double,8>> states;
	states.resize(count);

	for(int n=0;n<count;n++)
	{
		const libconfig::Setting& range = ranges[n];
		double rho, u, v, w, p, Bx, By, Bz;

		if(!(      range.lookupValue("rho", rho)
			&& range.lookupValue("u", u)
			&& range.lookupValue("v", v)
			&& range.lookupValue("w", w)
			&& range.lookupValue("p", p)
			&& range.lookupValue("Bx", Bx)
			&& range.lookupValue("By", By)
			&& range.lookupValue("Bz", Bz)))
			continue;

		states[n] = {rho, u, v, w, p, Bx, By, Bz};

	}

	if(m_test==3 || m_test==5)
	{
		double norm = sqrt(4.0*M_PI);

		states[0][5] = states[0][5] / norm;
		states[0][6] = states[0][6] / norm;
		states[0][7] = states[0][7] / norm;

		states[1][5] = states[1][5] / norm;
		states[1][6] = states[1][6] / norm;
		states[1][7] = states[1][7] / norm;
	}

	for(int i=0;i<nCellsX;i++)
	{
		double x = x0 + (i + 0.5) * dx;
		xCentroids[i] = x;

		for(int j=0;j<nCellsY;j++)
		{
			double y = y0 + (j + 0.5) * dy;
			yCentroids[j] = y;

			switch (m_test) 
			{
				case 1 : 
					{
						if(x < (x1-x0)/2 && x >= x0)
							Q[i][j] = fvCell(states[0]).toCons(gamma);
						else if(x < x1 && x >= (x1-x0)/2)
							Q[i][j] = fvCell(states[1]).toCons(gamma);
					}
					break;
				case 2 : 
					{
						if(x < (x1-x0)/2 && x >= x0)
							Q[i][j] = fvCell(states[0]).toCons(gamma);
						else if(x < x1 && x >= (x1-x0)/2)
							Q[i][j] = fvCell(states[1]).toCons(gamma);
					}
					break;
				case 3 : 
					{
						if(x < (x1-x0)/2 && x >= x0)
							Q[i][j] = fvCell(states[0]).toCons(gamma);
						else if(x < x1 && x >= (x1-x0)/2)
							Q[i][j] = fvCell(states[1]).toCons(gamma);
					}
					break;
				case 4 : 
					{
						if(x < (x1-x0)/2 && x >= x0)
							Q[i][j] = fvCell(states[0]).toCons(gamma);
						else if(x < x1 && x >= (x1-x0)/2)
							Q[i][j] = fvCell(states[1]).toCons(gamma);
					}
					break;
				case 5 : 
					{
						if(x < (x1-x0)/2 && x >= x0)
							Q[i][j] = fvCell(states[0]).toCons(gamma);
						else if(x < x1 && x >= (x1-x0)/2)
							Q[i][j] = fvCell(states[1]).toCons(gamma);
					}
					break;
			}
					

		}
	}
}

// Utility to output to a ofstream when needed
void fvSim::output()
{
	for(int i=0;i<nCellsX;i++)
	{
		for(int j=0;j<nCellsY;j++)
		{
			fvCell Wi = Q[i][j].toPrim(gamma);

			outputFile << xCentroids[i] << " " 
				   << yCentroids[j] << " "
				   << Wi[0] << " "
				   << Wi[1] << " "
				   << Wi[2] << " "
				   << Wi[3] << " "
				   << (Wi[4] + 0.5*(Wi[5]*Wi[5]+Wi[6]*Wi[6]+Wi[7]*Wi[7])) << " "
				   << Wi[5] << " "
				   << Wi[6] << " "
				   << Wi[7] << " "
				   << Wi[4] / (Wi[0]*(gamma-1.0)) << " "
				   << Q[i][j][4]
				   << std::endl;
		}
		//outputFile << "\n";
	}
}

// Calculate flux is used to parse any chosen flux function in based on user
// choice or testing, to reduce repetition of code.
const fvCell fvSim::calculateFlux(const fvCell& QL, const fvCell& QR, const double& dx,
		const double& dt, const double& gamma, bool x_dir,
		std::function<const fvCell(const fvCell&,const fvCell&, const double&, const double&, const double&, bool)> func){
	return func(QL, QR, dx, dt, gamma, x_dir);
}

// Reconstruct Data is used to activate linear data reconstruction or not.
std::array<fvCell,2> fvSim::reconstructData(const fvCell& QL, const fvCell& Qi, const fvCell& QR,
		std::function<std::array<fvCell,2>(const fvCell&,const fvCell&, const fvCell&)> func){
	return func(QL, Qi, QR);
}
