// Author - Oliver Marx ojm40@cam.ac.uk

// Standard packages
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstddef>
#include <iterator>
#include <math.h>
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
fvSim::fvSim(const char* configFileName, int testNum, std::string solver)
{
	m_solver = solver;
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
			&& test.lookupValue("test", m_testName)
			&& test.lookupValue("gamma", gamma)))
			std::cout << "Settings for test #" << testNum+1 << " read in."
				<< std::endl;

		m_test = testNum+1;

		const libconfig::Setting& range = tests[testNum]["inits"];
		std::cout << "Solv: " << m_solver << " Test: " << m_testName 
			<< " nCellsX " << nCellsX 
			<< " nCellsY " << nCellsY 
			<< " x0 " << x0 
			<< " x1 " << x1 
			<< " y0 " << y0 
			<< " y1 " << y1 
			<< " t0 " << t0 
			<< " t1 " << t1 
			<< " CFL " << CFL 
			<< " gamma " << gamma
			<< " test " << m_test 
			<< std::endl;

		// Create output file name and convert to char*
		std::string s = "/home/ojm40/Documents/MPhil_MHD/dat/output_" + m_testName + ".dat";
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

	do
	{
		computeTimeStep();
		t = t + dt;
		std::cout << "+" << std::flush;
		std::cout << "dt " << dt << std::endl;

		bool x_dir = true;

		// Copy the actual domain into the ghost cell domain
		copyDomain();

		// LENGTHS
		setBoundaryConditions(x_dir);

		// Qi stores left and right states
		// Initially just for the x direction reconstruction and flux
		// calculations
		std::array<fvCell,2> Qi;

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

/************************************** y-dir ************************************/

		x_dir = false;

		// Copy the actual domain into the ghost cell domain
		copyDomain();

		// BREADTHS
		setBoundaryConditions(x_dir);

		for(int i=0;i<nCellsX;i++)
		{
			std::array<fvCell,2> Qiplus1 = reconstructData(Q_bc[i+m_ghost][0],
								       Q_bc[i+m_ghost][1], 
								       Q_bc[i+m_ghost][2],
								       g );
				
			for(int j=1;j<nCellsY+m_ghost;j++)
			{

				Qi = Qiplus1;

				Qiplus1 = reconstructData(Q_bc[i+m_ghost][j],
						          Q_bc[i+m_ghost][j+1],
							  Q_bc[i+m_ghost][j+2], 
							  g);

				const fvCell dFL = (F(Qi[1], gamma, x_dir) - F(Qi[0], gamma, x_dir));
				const fvCell dFR = (F(Qiplus1[1], gamma, x_dir) - F(Qiplus1[0], gamma, x_dir));

				const fvCell QL = Qi[1] - (0.5 * dt / dx) * dFL;
				const fvCell QR = Qiplus1[0] - (0.5 * dt / dx) * dFR;

				f_half[i][j-1] = calculateFlux(QL, QR, dx, dt, gamma, x_dir, f);

			}
		}
		
		fullTimeStepUpdate(x_dir);

		double tFrac = t/t1;
		double scale = pow(10.0, ceil(log10(fabs(tFrac))) + 2);
		tFrac = round(tFrac * scale) / scale;

		if(tFrac >= printNumber * fracGap && !printedAtTime)
		{
			//outputFile << "\n\n";
			//outputFile << "t=" << t << "\n";
			//output();

			printedAtTime = true;
			std::cout << printNumber*10 << "%" << std::endl;
		}
		else if(printedAtTime && tFrac > printNumber * fracGap)
		{
			printNumber+=1.0;
			printedAtTime = false;
		}

	} while (t < t1);

	outputFile << "\n\n";
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

	if(nCellsX>1)
	{
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
	}

	if(nCellsY>1)
	{
		for(int i=0;i<nCellsX;i++)
		{
			for(int j=0;j<nCellsY-1;j++)
			{
				waveSpeedEstimatesMHD(Q[i][j], Q[i][j+1], SL, SR, gamma, false);
				SL = fabs(SL);
				SR = fabs(SR);
				max = std::max(max, std::max(SL,SR));;
			}
		}
	}

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
	std::cout << "dx = " << dx << " dy = " << dy << std::endl;

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
				// discontinuity halfway along x-axis
				case 1: case 2: case 6:
				{
					if(x < (x1-x0)/2 && x >= x0)
						Q[i][j] = fvCell(states[0]).toCons(gamma);
					else if(x < x1 && x >= (x1-x0)/2)
						Q[i][j] = fvCell(states[1]).toCons(gamma);
					break;
				}
				// discontinuity halfway along y-axis
				case 3 : case 7:
				{
					if(y < (y1-y0)/2 && y >= x0)
						Q[i][j] = fvCell(states[0]).toCons(gamma);
					else if(y < y1 && y >= (y1-y0)/2)
						Q[i][j] = fvCell(states[1]).toCons(gamma);
					break;
				}
				// discontinuity halfway along y=1-x
				case 4 :
				{
					if((y1-x) >= y)
						Q[i][j] = fvCell(states[0]).toCons(gamma);
					else if((y1-x) < y)
						Q[i][j] = fvCell(states[1]).toCons(gamma);
					break;
				}
				// circular discontinuity at radial distance
				// r=0.4
				case 5 : 
				{
					double r = sqrt(x*x + y*y);
					if(r <= 0.4)
						Q[i][j] = fvCell(states[0]).toCons(gamma);
					else
						Q[i][j] = fvCell(states[1]).toCons(gamma);
					break;
				}
				// diagonal Brio and Wu test
				case 8 :
				{
					double ang = cos(M_PI/4);

					states[0][5] = -ang + 0.75*ang;
					states[0][6] = ang + 0.75*ang;

					states[1][5] = ang + 0.75*ang;
					states[1][6] = -ang + 0.75*ang;

					if((y1-x) >= y)
						Q[i][j] = fvCell(states[0]).toCons(gamma);
					else if((y1-x) < y)
						Q[i][j] = fvCell(states[1]).toCons(gamma);
					break;
				}
				// Orszag-Tang vortex
				case 9 :
				{
					states[0][0] = gamma*gamma;
					states[0][1] = -sin(2*M_PI*y);
					states[0][2] = sin(2*M_PI*x);
					states[0][3] = 0.0;
					states[0][4] = gamma;
					states[0][5] = -sin(2*M_PI*y);
					states[0][6] = sin(4*M_PI*x);
					states[0][7] = 0.0;

					Q[i][j] = fvCell(states[0]).toCons(gamma);
					break;
				}
				// Kelvin-Helmholtz Instability
				case 10 :
				{
					double theta = M_PI/3;
					double M     = 1.0;
					double ca    = 0.1;
					double sigma = 0.1;

					states[0][0] = 1.0;
					states[0][1] = M/(2*tanh(20*y));
					states[0][2] = 0.1*sin(2*M_PI*x)*exp(-y*y/(theta*theta));
					states[0][3] = 0.0;
					states[0][4] = 1/gamma;
					states[0][5] = ca*sqrt(1.0)*cos(theta);
					states[0][6] = 0.0;
					states[0][7] = ca*sqrt(1.0)*sin(theta);

					Q[i][j] = fvCell(states[0]).toCons(gamma);
					break;
				}
			}
					

		}
	}
}

void fvSim::copyDomain()
{
	for(int i=0;i<nCellsX;i++)
	{
		for(int j=0;j<nCellsY;j++)
		{
			Q_bc[i+m_ghost][j+m_ghost] = Q[i][j];
		}
	}
}

void fvSim::setBoundaryConditions(bool x_dir)
{
	switch (m_test)
	{
		case 1 : case 2: case 3 : case 4 : case 5 : case 6 : case 7 : case 8 :
		{
			if(x_dir)
			{
				// LENGTHS
				for(int i=0;i<m_ghost;i++)
				{
					for(int j=0;j<nCellsY;j++)
					{
						Q_bc[i][j+m_ghost] = Q[0][j];
						Q_bc[i + m_ghost + nCellsX][j + m_ghost] = Q[nCellsX-1][j];
					}
				}
			}
			else {
				// BREADTHS
				for(int i=0;i<nCellsX;i++)
				{
					for(int j=0;j<m_ghost;j++)
					{
						Q_bc[i+m_ghost][j] = Q[i][0];
						Q_bc[i + m_ghost][j + m_ghost + nCellsY] = Q[i][nCellsY-1];
					}
				}
			}
			break;
		}
		case 9 :
		{
			if(x_dir)
			{
				// LENGTHS
				for(int i=0;i<m_ghost;i++)
				{
					for(int j=0;j<nCellsY;j++)
					{
						Q_bc[i][j+m_ghost] = Q[nCellsX-1][j];
						Q_bc[i + m_ghost + nCellsX][j + m_ghost] = Q[0][j];
					}
				}
			}
			else {
				// BREADTHS
				for(int i=0;i<nCellsX;i++)
				{
					for(int j=0;j<m_ghost;j++)
					{
						Q_bc[i+m_ghost][j] = Q[i][nCellsY-1];
						Q_bc[i + m_ghost][j + m_ghost + nCellsY] = Q[i][0];
					}
				}
			}
			break;
		}
	}
	
}

// Utility to output to a ofstream when needed
void fvSim::output()
{
	switch (m_test)
	{
		// Brio & Wu 1D test -- Print along x-axis, j=0
		case 1 : 
		{
			for(int i=0;i<nCellsX;i++)
			{
				fvCell Wi = Q[i][0].toPrim(gamma);

				outputFile << xCentroids[i] << " " 
					   << yCentroids[0] << " "
					   << Wi[0] << " "
					   << Wi[1] << " "
					   << Wi[2] << " "
					   << Wi[3] << " "
					   //<< Wi[4] << " "
					   << (Wi[4] + 0.5*(Wi[5]*Wi[5]+Wi[6]*Wi[6]+Wi[7]*Wi[7])) << " "
					   << Wi[5] << " "
					   << Wi[6] << " "
					   << Wi[7] << " "
					   << Wi[4] / (Wi[0]*(gamma-1.0)) << " "
					   << Q[i][0][4]
					   << std::endl;
			}
			break;
		}
		// Sod 2D test with x discontinuity -- Print along y=y1-y0/2
		case 2 : case 6 :
		{
			for(int i=0;i<nCellsX;i++)
			{
				fvCell Wi = Q[i][nCellsY/2].toPrim(gamma);

				outputFile << xCentroids[i] << " " 
					   << yCentroids[nCellsY/2] << " "
					   << Wi[0] << " "
					   << Wi[1] << " "
					   << Wi[2] << " "
					   << Wi[3] << " "
					   << (Wi[4] + 0.5*(Wi[5]*Wi[5]+Wi[6]*Wi[6]+Wi[7]*Wi[7])) << " "
					   << Wi[5] << " "
					   << Wi[6] << " "
					   << Wi[7] << " "
					   << Wi[4] / (Wi[0]*(gamma-1.0)) << " "
					   << Q[i][nCellsY/2][4]
					   << std::endl;
			}
			break;
		}
		// Sod 2D test with x discontinuity -- Print along x=x1-x0/2
		case 3 : case 7 :
		{
			for(int j=0;j<nCellsY;j++)
			{
				fvCell Wi = Q[nCellsX/2][j].toPrim(gamma);

				outputFile << xCentroids[nCellsX/2] << " " 
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
					   << Q[nCellsX/2][j][4]
					   << std::endl;
			}
			break;
		}
		// Sod 2D test with x=-y -- Print along x=y
		case 4 : case 8 :
		{
			for(int i=0;i<nCellsX;i++)
			{
				fvCell Wi = Q[i][i].toPrim(gamma);

				outputFile << xCentroids[i] << " " 
					   << yCentroids[i] << " "
					   << Wi[0] << " "
					   << Wi[1] << " "
					   << Wi[2] << " "
					   << Wi[3] << " "
					   << (Wi[4] + 0.5*(Wi[5]*Wi[5]+Wi[6]*Wi[6]+Wi[7]*Wi[7])) << " "
					   << Wi[5] << " "
					   << Wi[6] << " "
					   << Wi[7] << " "
					   << Wi[4] / (Wi[0]*(gamma-1.0)) << " "
					   << Q[i][i][4] << " "
					   << cos(M_PI/4)*Wi[1] + cos(M_PI/4)*Wi[2] << " "
					   << cos(M_PI/4)*Wi[2] - cos(M_PI/4)*Wi[1] << " "
					   << sqrt(2*xCentroids[i]*xCentroids[i]) << " "
					   << cos(M_PI/4)*Wi[6] - cos(M_PI/4)*Wi[5] << " "
					   << std::endl;
			}
			break;
		}
		// Sod 2D test -- Cynlindrical Explosion
		case 5 :
		{
			for(int i=0;i<nCellsX;i++)
			{
				for(int j=0;j<nCellsY;j++)
				{
					fvCell Wi = Q[i][j].toPrim(gamma);

					outputFile << xCentroids[i] << " " 
						   << yCentroids[j] << " "
						   << Wi[0] << " "
						   << sqrt(Wi[1]*Wi[1] + Wi[2]*Wi[2]) << " "
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
				outputFile << "\n";
			}
			break;
		}
		// Orszag Tang Vortex
		case 9 :
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
						   << sqrt(Wi[1]*Wi[1] + Wi[2]*Wi[2]) << " "
						   << std::endl;
				}
				outputFile << "\n";
			}
			break;
		}
		// Kelvin-Helmholtz Instability
		case 10 :
		{
			for(int i=0;i<nCellsX;i++)
			{
				for(int j=0;j<nCellsY;j++)
				{
					fvCell Wi = Q[i][j].toPrim(gamma);

					outputFile << xCentroids[i] << " " 
						   << yCentroids[j] << " "
						   << Wi[0] << " "
						   << sqrt(Wi[1]*Wi[1] + Wi[2]*Wi[2]) << " "
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
				outputFile << "\n";
			}
			break;
		}
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
