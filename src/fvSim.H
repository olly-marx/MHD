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
		void fullTimeStepUpdate();
		const fvCell halfTimeStepUpdate(const fvCell&, const fvCell&);
		const fvCell F(const fvCell&);
		const fvCell LF_Flux(const fvCell&, const fvCell&);
		const fvCell Richt_Flux(const fvCell&);
		const fvCell FORCE_Flux(const fvCell&, const fvCell&);
		const fvCell HLL_Flux(const fvCell&, const fvCell&);
		void waveSpeedEstimates(const fvCell&, const fvCell&, double&, double&);
		std::array<fvCell,2> linearDataReconstruction(const fvCell&, const fvCell&,
				const fvCell&);
		std::array<fvCell,2> constantDataReconstruction(const fvCell&, const fvCell&,
				const fvCell&);
		std::array<fvCell,2> reconstructData(const fvCell&, const fvCell&, const fvCell&,
			std::function<std::array<fvCell,2>(const fvCell&,const fvCell&, const fvCell&)>);
		const fvCell calculateFlux(const fvCell&, const fvCell&, 
				std::function<const fvCell(const fvCell&,const fvCell&)>);
		double superbee(const double&);

		// Private member variables
		unsigned int    nCells;
		double          x0;
		double          x1;
		double          t0;
		double          t1;
		double          gamma;
		double          CFL;
		double          dx;
		double          dt;
		std::string     m_solver;

		// Main Grid and Temps
		std::vector<fvCell> Q;
		std::vector<fvCell> Q_i_nplus1;
		std::vector<fvCell> f_iplushalf_n;

		// Other files
		std::ofstream outputFile;
};

#endif
