// AUTHOR - Oliver Marx ojm40@cam.ac.uk

#include <array>

#include "fvCell.H"

#ifndef FLUXFUNCTIONS_HEADER
#define FLUXFUNTIONS_HEADER

const fvCell halfTimeStepUpdate(const fvCell& QL, const fvCell& QR, const double& dx, const double& dt, bool x_dir);
const fvCell F(const fvCell&, const double&, bool);
const fvCell LF_Flux(const fvCell& QL, const fvCell& QR, const double& dx, const double& dt, const double& gamma, bool x_dir);
const fvCell FORCE_Flux(const fvCell& QL, const fvCell& QR, const double& dx, const double& dt, const double& gamma, bool x_dir);
const fvCell HLL_Flux(const fvCell& QL, const fvCell& QR, const double& dx, const double& dt, const double& gamma, bool x_dir);
const fvCell HLLC_Flux(const fvCell& QL, const fvCell& QR, const double& dx, const double& dt, const double& gamma, bool x_dir);
void waveSpeedEstimates(const fvCell&, const fvCell&, double&, double&, const double&, bool);

#endif
