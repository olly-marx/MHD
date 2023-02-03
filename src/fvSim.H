// Author - Oliver Marx ojm40@cam.ac.uk

#include <functional>
#include <iostream>
#include <cstddef>
#include <string>
#include <vector>
#include <fstream>

#include <libconfig.h++>

class fvSim{
	public:
		fvSim(const char*, int, std::string);
		~fvSim();
		void run();
	private:
		//Private member functions
		void init(const libconfig::Setting&);
		void output();
		void computeTimeStep();
		void fullTimeStepUpdate();
		double halfTimeStepUpdate(double, double);
		double F(double);
		double LF_Flux(double, double);
		double Richt_Flux(double);
		double FORCE_Flux(double, double);
		double calculateFlux(double, double, std::function<double(double,double)>);

		// Private member variables
		unsigned int    nCells;
		double          x0;
		double          x1;
		double          t0;
		double          t1;
		double          a;
		double          CFL;
		double          dx;
		double          dt;
		std::string     m_solver;

		// Main Grid and Temps
		std::vector<double> Q;
		std::vector<double> Q_i_nplus1;
		std::vector<double> f_iplushalf_n;

		// Other files
		std::ofstream outputFile;
};
