//AUTHOR - Oliver Marx ojm40@cam.ac.uk

#include <iostream>
#include <cstddef>
#include <vector>

// Headers
#include "fvSim.H"

int main(int argc, char* argv[]){
	if(argc==1){
		std::cout << "No test case OR solver chosen. Exit." << std::endl;
	}
	else if(argc==2){
		std::string solver = argv[1];
		for(int test=1; test<6; test++){
			const char* fileName = "/home/ojm40/Documents/MPhil_MHD/stg/config.cfg";
			fvSim Sim(fileName, test-1, solver);
			Sim.run();
		}
	}
	else if(argc==3){
		std::string solver = argv[1];
		int test = strtol(argv[2], nullptr, 0);
		const char* fileName = "/home/ojm40/Documents/MPhil_MHD/stg/config.cfg";
		fvSim Sim(fileName, test-1, solver);
		Sim.run();
		
	}
	else{
		std::cout << "Too many arguments parsed. Fatal error. Exit." << std::endl;
	}
}
