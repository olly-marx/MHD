// AUTHOR - Oliver Marx ojm40@cam.ac.uk

#include <array>

#include "fvCell.H"

std::array<fvCell,2> constantDataReconstruction(const fvCell&, const fvCell&, const fvCell&);
std::array<fvCell,2> linearDataReconstruction(const fvCell&, const fvCell&, const fvCell&);
double superbee(const double&);
