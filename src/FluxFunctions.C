// AUTHOR - Oliver Marx ojm40@cam.ac.uk

#include <array>
#include <cmath>

#include "fvCell.H"
#include "FluxFunctions.H"

// The flux function for the PDE in conservation form. In this case it is the
// Burgers' Flux ---- To be replaced by Euler when ready
const fvCell F(const fvCell& Qi, const double& gamma, bool x_dir)
{
	const fvCell Wi = Qi.toPrim(gamma);
	
	double f0 = Wi[0] * (x_dir ? Wi[1] : Wi[2]);
	double f1 = x_dir ? 
		Wi[0] * Wi[1] * Wi[1] + Wi[3] :
		Wi[0] * Wi[1] * Wi[2];
	double f2 = x_dir ?
		Wi[0] * Wi[1] * Wi[2] :
		Wi[0] * Wi[2] * Wi[2] + Wi[3];
	double f3 = (Qi[3] + Wi[3]) * (x_dir ? Wi[1] : Wi[2]);

	return fvCell({f0, f1, f2, f3}, true);
}

// Compute a half time step update of the data in Q_i_n giving Q_i+1/2_n+1/2
const fvCell halfTimeStepUpdate (const fvCell& QL, const fvCell& QR, const double& dx,
		const double& dt, const double& gamma, bool x_dir)
{
	return 0.5 * (QR + QL) - (0.5 * dt / dx) * (F(QR, gamma, x_dir) - F(QL, gamma, x_dir));
}

// Lax-Friedrichs Flux calculation, takes left amd right states at cell boundary
// gives the flux at a cell boundary f_i+1/2_n
const fvCell LF_Flux (const fvCell& QL, const fvCell& QR, const double& dx,
		const double& dt, const double& gamma, bool x_dir)
{
	return 0.5 * (dx / dt) * (QL - QR) + 0.5 * (F(QR, gamma, x_dir) + F(QL, gamma, x_dir));
}

// FORCE Flux is the average of a LF and Richt Flux. So, we find those at
// Q_i+1/2_n and then calculate tye flux.
const fvCell FORCE_Flux(const fvCell& QL, const fvCell& QR, const double& dx, const double& dt, const double& gamma, bool x_dir)
{
	const fvCell LF = LF_Flux(QL, QR, dx, dt, gamma, x_dir);
	const fvCell Q_half = halfTimeStepUpdate(QL, QR, dx, dt, gamma, x_dir);
	const fvCell Richt = F(Q_half, gamma, x_dir);
	return 0.5 * (LF + Richt);
}

const fvCell HLL_Flux(const fvCell& QL, const fvCell& QR, const double& dx, const double& dt, const double& gamma, bool x_dir)
{
	double SL, SR;
	waveSpeedEstimates(QL, QR, SL, SR, gamma, x_dir);

	const fvCell fL = F(QL, gamma, x_dir);                                                        
	const fvCell fR = F(QR, gamma, x_dir);                                                        
	const fvCell fHLL = (1.0 / (SR-SL)) * (SR*fL - SL*fR + SL * SR * (QR - QL));

       	if(SL >= 0)
	{
		return fL;
	}else if(SR <= 0)
	{
       	        return fR; 
	}else
	{
       	        return fHLL;
	}
       	                                                                          
}                    

const fvCell HLLC_Flux(const fvCell& QL, const fvCell& QR, const double& dx, const double& dt, const double& gamma, bool x_dir)
{
	double SL, SR, Sstar;
	waveSpeedEstimates(QL, QR, SL, SR, gamma, x_dir);

	const fvCell WL = QL.toPrim(gamma);
	const fvCell WR = QR.toPrim(gamma);

	const double& rhoL  = QL[0],
		      rhoR  = QR[0],
		      L    = x_dir ? WL[1] : WL[2],
		      R    = x_dir ? WR[1] : WR[2],
		      pL    = WL[3],
		      pR    = WR[3],
		      rhouL = QL[1],
		      rhouR = QR[1],
		      EL    = QL[3],
		      ER    = QR[3];

	Sstar = (pR - pL + rhouL * (SL - L) - rhouR * (SR - R)) 
		/ (rhoL * (SL - L) - rhoR * (SR - R));

	double prefixL = rhoL * (SL - L) / (SL - Sstar);
	double prefixR = rhoR * (SR - R) / (SR - Sstar);

	double ELstar = EL / rhoL + (Sstar - L) * (Sstar + pL / (rhoL * (SL - L)));
	double ERstar = ER / rhoR + (Sstar - R) * (Sstar + pR / (rhoR * (SR - R)));

	double uL    = x_dir ? Sstar : WL[1],
	       uR    = x_dir ? Sstar : WR[1],
	       vL    = x_dir ? WL[2] : Sstar,
	       vR    = x_dir ? WR[2] : Sstar;

	const fvCell QstarL = prefixL * fvCell({1.0, uL, vL, ELstar}, true);
	const fvCell QstarR = prefixR * fvCell({1.0, uR, vR, ERstar}, true);

	const fvCell fL = F(QL, gamma, x_dir);                                                        
	const fvCell fR = F(QR, gamma, x_dir);                                                        
	const fvCell fstarL = fL + SL * (QstarL - QL);
	const fvCell fstarR = fR + SR * (QstarR - QR);

       	if(SL >= 0)
	{
		return fL; 
	}else if (SL <= 0.0 && Sstar >= 0.0)
	{
       	        return fstarL;
	}else if (Sstar <= 0.0 && SR >= 0.0)
	{
		return fstarR;
	}else
	{
		return fR;
	}
}

void waveSpeedEstimates(const fvCell& QL, const fvCell& QR, double& SL, double& SR, const double& gamma, bool x_dir)
{
       fvCell WL = QL.toPrim(gamma);                                             
       fvCell WR = QR.toPrim(gamma);                                             

       double aL  = sqrt(gamma * WL[3] / WL[0]);                                
       double aR  = sqrt(gamma * WR[3] / WR[0]);                                

       const double& rhoL  = QL[0],
		     rhoR  = QR[0],
		     uL    = WL[1],
		     uR    = WR[1],
		     vL    = WL[2],
		     vR    = WR[2],
		     pL    = WL[3],
		     pR    = WR[3];

	double L = x_dir ? uL : vL;
	double R = x_dir ? uR : vR;

	double rhobar = 0.5 * (rhoL + rhoR);
	double abar   = 0.5 * (aL + aR);
	double ppvrs  = 0.5 * (pL + pR) - 0.5 * (R - L) * rhobar * abar;
	double pstar  = std::max(0.0, ppvrs);

	double qL, qR;

	//std::cout << "pL " << pL << " pstar " << pstar << " pR " << pR << std::endl; 

	if(pstar <= pL)
	{
		qL = 1.0;
	}
	else if(pstar > pL)
	{
		qL = sqrt(1.0 + ((gamma + 1.0) / (2.0 * gamma)) * (pstar / pL - 1.0));
	}

	if(pstar <= pR)
	{
		qR = 1.0;
	}
	else if(pstar > pR)
	{
		qR = sqrt(1.0 + ((gamma + 1.0) / (2.0 * gamma)) * (pstar / pR - 1.0));
	}

	SL = L - aL * qL;
	SR = R + aR * qR;
}
