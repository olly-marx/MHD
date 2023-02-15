// Author - Oliver Marx ojm40@cam.ac.uk

#include <functional>
#include <iostream>
#include <cstddef>
#include <string>
#include <vector>
#include <fstream>

#include <libconfig.h++>

#include "fvCell.H"

#ifndef FVSIM_HEADER
#define FVSIM_HEADER

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
		void fullTimeStepUpdate(bool);
		std::array<fvCell,2> reconstructData(const fvCell&, const fvCell&, const fvCell&,
			std::function<std::array<fvCell,2>(const fvCell&,const fvCell&, const fvCell&)>);
		const fvCell calculateFlux(const fvCell&, const fvCell&, const double&, const double&, const double&, bool, 
					std::function<const fvCell(const fvCell&, const fvCell&, const double&, const double&, 
					const double&, bool)>);

		// Private member variables
		int             nCells;
		int             m_ghost;
		int             m_test;
		double          x0;
		double          x1;
		double          y0;
		double          y1;
		double          t0;
		double          t1;
		double          gamma;
		double          CFL;
		double          dx;
		double          dy;
		double          dt;
		double          t;
		std::string     m_solver;

		// Main Grid and Temps
		std::vector<std::vector<fvCell>> Q;
		std::vector<std::vector<fvCell>> Q_new;
		std::vector<std::vector<fvCell>> Q_bc;
		std::vector<std::vector<fvCell>> f_half;

		std::vector<double> xCentroids;
		std::vector<double> yCentroids;

		// Other files
		std::ofstream outputFile;
};

#endif
