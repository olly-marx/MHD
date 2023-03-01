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

	const double& rho = Wi[0],
	              u   = Wi[1],
		      v   = Wi[2],
		      w   = Wi[3],
		      p   = Wi[4],
		      E   = Qi[4],
		      Bx  = Wi[5],
		      By  = Wi[6],
		      Bz  = Wi[7];

	double q = x_dir ? u : v;
	double Bn = x_dir ? Bx : By;

	std::array<double,3> n = {(x_dir ? 1.0 : 0.0), (x_dir ? 0.0 : 1.0), 0.0};
	
	double f0 = rho * q;

	double f1 = rho*q*u + n[0] * (p + 0.5*(Bx*Bx + By*By + Bz*Bz)) - Bn*Bx;
	double f2 = rho*q*v + n[1] * (p + 0.5*(Bx*Bx + By*By + Bz*Bz)) - Bn*By;
	double f3 = rho*q*w + n[2] * (p + 0.5*(Bx*Bx + By*By + Bz*Bz)) - Bn*Bz;

	double f4 = (E + p + 0.5*(Bx*Bx + By*By + Bz*Bz))*q - (u*Bx + v*By + w*Bz)*Bn;

	double f5 = x_dir ?
		0.0 :
		Bx*v - By*u;
	double f6 = x_dir ?
		By*u - Bx*v :
		0.0;
	double f7 = x_dir ?
		Bz*u - Bx*w :
		Bz*v - By*w;

	return fvCell({f0, f1, f2, f3, f4, f5, f6, f7}, true);
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

const fvCell HLLC_Flux(const fvCell& QL, const fvCell& QR, const double& dx, const double& dt, const double& gamma, bool x_dir)
{
	double SL=0.0, SR=0.0;
	waveSpeedEstimatesMHD(QL, QR, SL, SR, gamma, x_dir);

	const fvCell fL = F(QL, gamma, x_dir);                                                        
	const fvCell fR = F(QR, gamma, x_dir);                                                        

	if(SL >= 0)
		return fL;
	else if(SR <= 0)
		return fR;

	const fvCell WL = QL.toPrim(gamma);
	const fvCell WR = QR.toPrim(gamma);

	const double& rhoL  = WL[0],
		      rhoR  = WR[0],
		      uL    = WL[1],
		      uR    = WR[1],
		      vL    = WL[2],
		      vR    = WR[2],
		      wL    = WL[3],
		      wR    = WR[3],
		      EL    = QL[4],
		      ER    = QR[4],
		      pLgas = WL[4],
		      pRgas = WR[4],
		      BxL   = WL[5],
		      BxR   = WR[5],
		      ByL   = WL[6],
		      ByR   = WR[6],
		      BzL   = WL[7],
		      BzR   = WR[7];

	double pL = pLgas + 0.5 * (BxL*BxL + ByL*ByL + BzL*BzL);
	double pR = pRgas + 0.5 * (BxR*BxR + ByR*ByR + BzR*BzR);
	
	const double& qL = x_dir ? uL : vL,
	              qR = x_dir ? uR : vR;
	const double& BnL = x_dir ? BxL : ByL,
	              BnR = x_dir ? BxR : ByR;

	double qstar = (rhoR*qR*(SR - qR) - rhoL*qL*(SL - qL) + pL - pR - BnL*BnL + BnR*BnR) / (rhoR*(SR - qR) - rhoL*(SL - qL));

	const fvCell QHLL = (1.0 / (SR - SL)) * (SR*QR - SL*QL - fR + fL);

	double BxHLL  = QHLL[5];
	double ByHLL  = QHLL[6];
	double BzHLL  = QHLL[7];

	double BnHLL = x_dir ? BxHLL : ByHLL;

	double pstar   = rhoL*(SL - qL)*(qstar - qL) + pL - BnL*BnL + BnHLL*BnHLL;

	double rho, u, v, w, p, E, Bx, By, Bz, SK, q, Bn;
	fvCell F, U;

	if(qstar >= 0.0)
	{
		rho  = rhoL,
		u    = uL,
		v    = vL,
		w    = wL,
		p    = pL,
		E    = EL,
		Bx   = BxL,
		By   = ByL,
		Bz   = BzL,
		SK   = SL,
		q    = qL,
		Bn   = BnL,
		U    = QL,
		F    = fL;
	}
	else if(qstar < 0.0)
	{
		rho  = rhoR,
		u    = uR,
		v    = vR,
		w    = wR,
		p    = pR,
		E    = ER,
		Bx   = BxR,
		By   = ByR,
		Bz   = BzR,
		SK   = SR,
		q    = qR,
		Bn   = BnR,
		U    = QR,
		F    = fR;
	}

	double rhostar = rho*(SK - q) / (SK - qstar);

	double rhoUstar = x_dir ? 
			  rhostar * qstar :
			  (rho*u)*(SK-q)/(SK-qstar) - (BnHLL*BxHLL - Bn*Bx)/(SK-qstar);

	double rhoVstar = x_dir ? 
			  (rho*v)*(SK-q)/(SK-qstar) - (BnHLL*ByHLL - Bn*By)/(SK-qstar):
			  rhostar * qstar;
	
	double rhoWstar = (rho*w)*(SK-q)/(SK-qstar) - (BnHLL*BzHLL - Bn*Bz)/(SK-qstar);

	double rhoHLL   = QHLL[0];

	double rhouHLL  = QHLL[1];
	double rhovHLL  = QHLL[2];
	double rhowHLL  = QHLL[3];

	double uHLL  = rhouHLL / rhoHLL;
	double vHLL  = rhovHLL / rhoHLL;
	double wHLL  = rhowHLL / rhoHLL;

	double BdotUStar   = BxHLL*uHLL + ByHLL*vHLL + BzHLL*wHLL;

	double BdotU   = Bx*u + By*v + Bz*w;

	double Estar    = E*(SK-q)/(SK-qstar) + ((pstar*qstar - p*q) - (BnHLL*BdotUStar - Bn*BdotU)) / (SK - qstar);

	fvCell Ustar = fvCell({rhostar, rhoUstar, rhoVstar, rhoWstar, Estar, BxHLL, ByHLL, BzHLL}, true);

	fvCell Fstar = F + SK*(Ustar - U);

	return Fstar;
}

const fvCell HLL_Flux(const fvCell& QL, const fvCell& QR, const double& dx, const double& dt, const double& gamma, bool x_dir)
{
	double SL=0.0, SR=0.0;
	waveSpeedEstimatesMHD(QL, QR, SL, SR, gamma, x_dir);

	const fvCell fL = F(QL, gamma, x_dir);                                                        
	const fvCell fR = F(QR, gamma, x_dir);                                                        

	const fvCell fHLL = (1.0 / (SR - SL)) * (SR*fR - SL*fL + SL*SR*(QR - QL));

	if(SL >= 0.0)
		return fL;
	else if(SR > 0.0)
		return fHLL;
	else 
		return fR;
}

void waveSpeedEstimatesMHD(const fvCell& QL, const fvCell& QR, double& SL, double& SR, const double& gamma, bool x_dir)
{
       fvCell WL = QL.toPrim(gamma);                                             
       fvCell WR = QR.toPrim(gamma);                                             

       const double& rhoL  = WL[0],
        	     rhoR  = WR[0],
        	     pLgas = WL[4],
        	     pRgas = WR[4],
        	     BxL   = WL[5],
        	     BxR   = WR[5],
        	     ByL   = WL[6],
        	     ByR   = WR[6],
        	     BzL   = WL[7],
        	     BzR   = WR[7];

       double pL = pLgas + 0.5 * (BxL*BxL + ByL*ByL + BzL*BzL);
       double pR = pRgas + 0.5 * (BxR*BxR + ByR*ByR + BzR*BzR);

       double BdotBL = BxL*BxL + ByL*ByL + BzL*BzL;
       double BdotBR = BxR*BxR + ByR*ByR + BzR*BzR;

       double qL = x_dir ? WL[1] : WL[2];
       double qR = x_dir ? WR[1] : WR[2];

       double BnL = x_dir ? BxL : ByL;
       double BnR = x_dir ? BxR : ByR;

       double GL = (gamma*pL + BdotBL) / rhoL ;
       double GR = (gamma*pR + BdotBR) / rhoR ;

       double cfL = sqrt(0.5*(GL + sqrt(GL*GL - 4.0*gamma*pL*BnL*BnL / (rhoL*rhoL))));
       double cfR = sqrt(0.5*(GR + sqrt(GR*GR - 4.0*gamma*pR*BnR*BnR / (rhoR*rhoR))));

       double SWL = fabs(qL + cfL);
       double SWR = fabs(qR + cfR);

       SR = std::max(SWL, SWR);
       SL = -SR;
}
