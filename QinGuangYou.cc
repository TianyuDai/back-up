/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 * 
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#include "QinGuangYou.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include <string>

#include "tinyxml2.h"
#include<iostream>

#include "FluidDynamics.h"
#define hbarc 0.197327053

#define MAGENTA "\033[35m"

using namespace Jetscape;

using std::ofstream;
using std::ifstream;
using std::ostream;
using std::ios;


const double QS = 1.0;

QinGuangYou::QinGuangYou()
{
  SetId("QinQuangYou");
  VERBOSE(8);

  //vectors for elastic rates:
  dGamma_qq = new vector<double>;
  dGamma_qg = new vector<double>;
  dGamma_qq_q = new vector<double>;
  dGamma_qg_q = new vector<double>;
}

QinGuangYou::~QinGuangYou()
{
  VERBOSE(8);
}

void QinGuangYou::Init()
{
  INFO<<"Intialize QinGuangYou ...";

  // Redundant (get this from Base) quick fix here for now
  tinyxml2::XMLElement *eloss= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" );
  if ( !eloss )     throw std::runtime_error("Eloss not properly initialized in XML file ...");

  tinyxml2::XMLElement *qinguangyou=eloss->FirstChildElement("QinGuangYou");
  // check that all is there
  if ( !qinguangyou )     throw std::runtime_error("QinGuangYou not properly initialized in XML file ...");
  if ( !qinguangyou->FirstChildElement( "name" ) )     throw std::runtime_error("QinGuangYou not properly initialized in XML file ...");
  if ( !qinguangyou->FirstChildElement( "alpha_s" ) )     throw std::runtime_error("QinGuangYou not properly initialized in XML file ...");
  if ( !qinguangyou->FirstChildElement( "pcut" ) )     throw std::runtime_error("QinGuangYou not properly initialized in XML file ...");

  string s = qinguangyou->FirstChildElement( "name" )->GetText();
  JSDEBUG << s << " to be initilizied ...";

  alpha_s = 0.3;
  qinguangyou->FirstChildElement("alpha_s")->QueryDoubleText(&alpha_s);
    
  pcut = 2.0;
  qinguangyou->FirstChildElement("pcut")->QueryDoubleText(&pcut);

  g = sqrt(4.*M_PI*alpha_s);
  alpha_em = 1./137.;

  // Path to additional data
  if ( !qinguangyou->FirstChildElement( "path" ) )     throw std::runtime_error("QinGuangYou not properly initialized in XML file ...");
  PathToTables=qinguangyou->FirstChildElement( "path" )->GetText();

  // Initialize random number distribution
  ZeroOneDistribution = uniform_real_distribution<double> { 0.0, 1.0 };
  
  readRadiativeRate(&dat, &Gam);
  readElasticRateOmega();
  readElasticRateQ();
}

void QinGuangYou::DoEnergyLoss(double deltaT, double Time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
{
  VERBOSESHOWER(5)<< MAGENTA << "SentInPartons Signal received : "<<deltaT<<" "<<Q2<<" "<< pIn.size();

  // particle info
  int Id, newId;
  double pAbs, px, py, pz;   // momentum for initial parton (pIn)
  double pRest, pxRest;      // momentum in the rest frame of fluid cell (pIn)
  double pyRest, pzRest;
  double k, kRest;           // momentum for radiated parton (pOut)
  double pNew, pxNew;        // momentum for final parton (pOut)
  double pyNew, pzNew;
  double pNewRest;           // momentum in the rest frame of fluid cell (pOut)
  double q;           // transferred energy/momentum for scattering
  double xx, yy, zz;         // position of initial parton (pIn)
  FourVector pVec, pVecNew;  // 4 vectors for momenta before & after process
  FourVector pVecRest;       // 4 vector in the rest frame of fluid cell
  FourVector pVecNewRest;
  FourVector kVec;           // 4 vector for momentum of radiated particle
  FourVector xVec;           // 4 vector for position (for next time step!)
  double eta;                // pseudo-rapidity
  
  // flow info
  double vx, vy, vz;         // 3 components of flow velocity
  double T;                  // Temperature of fluid cell
  double beta, gamma;        // flow velocity & gamma factor
  double cosPhi;             // angle between flow and particle
  double cosPhiRest;         // angle between flow and particle in rest frame
  double boostBack;          // factor for boosting back to lab frame
  double cosPhiRestEl;       // angle between flow and scat. particle in rest frame
  double boostBackEl;
  
  
  for (int i=0;i<pIn.size();i++) {
    // Only accept low t particles
    if (pIn[i].t() > QS*QS + rounding_error) continue;
    TakeResponsibilityFor ( pIn[i] ); // Generate error if another module already has responsibility.
    
    Id = pIn[i].pid();

    px = pIn[i].px();
    py = pIn[i].py();
    pz = pIn[i].pz();
    // In QinGuangYou, particles are all massless and on-shell
    pAbs = sqrt(px*px+py*py+pz*pz);
    pVec = FourVector ( px, py, pz, pAbs );

    xx = pIn[i].x_in().x();
    yy = pIn[i].x_in().y();
    zz = pIn[i].x_in().z();

    eta = pIn[i].eta();

    std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;
    GetHydroCellSignal(Time, xx, yy, zz, check_fluid_info_ptr);
    VERBOSE(8)<< MAGENTA<<"Temperature from Brick (Signal) = "
	      <<check_fluid_info_ptr->temperature;

    vx = check_fluid_info_ptr->vx;
    vy = check_fluid_info_ptr->vy;
    vz = check_fluid_info_ptr->vz;
    T = check_fluid_info_ptr->temperature;

    beta = sqrt( vx*vx + vy*vy + vz*vz );

    // Set momentum in fluid cell's frame
    // 1: for brick
    if (beta < 1e-10)
      {
	gamma = 1.;
	cosPhi = 1.;
	pRest = pAbs;
	pVecRest = pVec;

	cosPhiRest = 1.;
	boostBack = 1.;
      }
    // 2: for evolving medium
    else
      {
	gamma  = 1./sqrt( 1. - beta*beta );
	cosPhi = ( px*vx + py*vy + pz*vz )/( pAbs*beta );

	// boost particle to the local rest frame of fluid cell
	pRest  = pAbs*gamma*( 1. - beta*cosPhi );

	pxRest = -vx*gamma*pAbs
	  + (1.+(gamma-1.)*vx*vx/(beta*beta))*px
	  + (gamma-1.)*vx*vy/(beta*beta)*py
	  + (gamma-1.)*vx*vz/(beta*beta)*pz;
	pyRest = -vy*gamma*pAbs
	  + (1.+(gamma-1.)*vy*vy/(beta*beta))*py
	  + (gamma-1.)*vx*vy/(beta*beta)*px
	  + (gamma-1.)*vy*vz/(beta*beta)*pz;
	pzRest = -vz*gamma*pAbs
	  + (1.+(gamma-1.)*vz*vz/(beta*beta))*pz
	  + (gamma-1.)*vx*vz/(beta*beta)*px
	  + (gamma-1.)*vy*vz/(beta*beta)*py;

	pVecRest = FourVector ( pxRest, pyRest, pzRest, pRest );

	cosPhiRest = ( pxRest*vx + pyRest*vy + pzRest*vz )/( pRest*beta );
	boostBack = gamma*( 1. + beta*cosPhiRest );
      }

    if (pRest < pcut) continue;

    xVec = FourVector( xx+px/pAbs*deltaT, yy+py/pAbs*deltaT, zz+pz/pAbs*deltaT,
		       Time+deltaT );

    int sign = DetermineSign(pRest, T, Id); 
    if (sign < 0) omega = -0.02; 
    else omega = 0.02; 
    int process = DetermineProcess(pRest, T, deltaT, Id, sign);
    // WARN<<"PROCESS "<<process; 
    VERBOSE(8)<< MAGENTA
	      << "Time = " << Time << " Id = " << Id
	      << " process = " << process << " T = " << T
	      << " pAbs = " << pAbs << " " << px << " " << py << " " << pz 
	      << " | position = " << xx << " " << yy << " " << zz;


    // Do nothing for this parton at this timestep
    if (process == 0) 
      {
	pOut.push_back(Parton(0, Id, 0, pVec, xVec));
	pOut[pOut.size()-1].set_form_time(0.);
	return;
      }
    if (std::abs(Id) == 1 || std::abs(Id) == 2 || std::abs(Id) == 3)
      {
	// quark radiating gluon (q->qg)
	if (process == 1)
	  {
	    if (pRest/T < AMYpCut) return;

	    // sample radiated parton's momentum
	    kRest = getNewMomentumRad(pRest, T, process);
            if(kRest > pRest) return;

	    // final state parton's momentum
	    pNewRest = pRest - kRest;

	    // if pNew is smaller than pcut, final state parton is
	    // absorbed into medium
	    if (pNewRest > pcut)
	      {
		pNew = pNewRest*boostBack;
		pVecNew.Set( (px/pAbs)*pNew, (py/pAbs)*pNew, (pz/pAbs)*pNew, pNew );
		pOut.push_back(Parton(0, Id, 0, pVecNew, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
	      }

	    if (kRest > pcut)
	      {
		k = kRest*boostBack;
		kVec.Set( (px/pAbs)*k, (py/pAbs)*k, (pz/pAbs)*k, k );
		pOut.push_back(Parton(0, Id, 0, kVec, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
	      }

	    return;
	  }
	// quark radiating photon (q->qgamma)
	else if (process == 2)
	  {
	    if (pRest/T < AMYpCut) return;

	    // sample radiated parton's momentum
	    kRest = getNewMomentumRad(pRest, T, process);
            if(kRest > pRest) return;

	    // final state parton's momentum
	    pNewRest = pRest - kRest;

	    // if pNew is smaller than pcut, final state parton is
	    // absorbed into medium
	    if (pNewRest > pcut)
	      {
		pNew = pNewRest*boostBack;
		pVecNew.Set( (px/pAbs)*pNew, (py/pAbs)*pNew, (pz/pAbs)*pNew, pNew );
		pOut.push_back(Parton(0, Id, 0, pVecNew, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
	      }

	    // photon doesn't have energy threshold; No absorption into medium
	    // However, we only keep positive energy photons
	    if (kRest > 0.)
	      {
		k = kRest*boostBack;
		kVec.Set( (px/pAbs)*k, (py/pAbs)*k, (pz/pAbs)*k, k );
		pOut.push_back(Parton(0, Id, 0, kVec, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
	      }

	    return;
	  }
	// quark scattering with either quark (qq->qq) or gluon (qg->qg)
	else if (process == 5 || process == 6)
	  {
	    //omega = getEnergyTransfer(pRest, T, process);
	    q = getMomentumTransfer(pRest, T, process);
	    pVecNewRest = getNewMomentumElas(pVecRest, q);

	    pNewRest = pVecNewRest.t();

	    // if pNew is smaller than pcut, final state parton is
	    // absorbed into medium
	    if (pNewRest > pcut)
	      {
		// Boost scattered particle to lab frame
		// 1: for brick
		if (beta < 1e-10)
		  {
		    pVecNew = pVecNewRest;
		  }
		// 2: for evolving medium
		else
		  {
		    pxNew = vx*gamma*pVecNewRest.t() 
		      + (1.+(gamma-1.)*vx*vx/(beta*beta))*pVecNewRest.x() 
		      + (gamma-1.)*vx*vy/(beta*beta)*pVecNewRest.y()
		      + (gamma-1.)*vx*vz/(beta*beta)*pVecNewRest.z();
		    pyNew = vy*gamma*pVecNewRest.t() 
		      + (1.+(gamma-1.)*vy*vy/(beta*beta))*pVecNewRest.y()
		      + (gamma-1.)*vx*vy/(beta*beta)*pVecNewRest.x()
		      + (gamma-1.)*vy*vz/(beta*beta)*pVecNewRest.z();
		    pzNew = vz*gamma*pVecNewRest.t() 
		      + (1.+(gamma-1.)*vz*vz/(beta*beta))*pVecNewRest.z()
		      + (gamma-1.)*vx*vz/(beta*beta)*pVecNewRest.x()
		      + (gamma-1.)*vy*vz/(beta*beta)*pVecNewRest.y();
			  
		    pNew = sqrt( pxNew*pxNew + pyNew*pyNew + pzNew*pzNew );
		    pVecNew.Set( pxNew, pyNew, pzNew, pNew );
		  }

		pOut.push_back(Parton(0, Id, 0, pVecNew, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
	      }

	    return;
	  }
	// quark converting to gluon
	else if (process == 9)
	  {
	    pOut.push_back(Parton(0, 21, 0, pVec, xVec));
	    pOut[pOut.size()-1].set_form_time(0.);

	    return;
	  }
	// quark converting to photon
	else if (process == 10)
	  {
	    pOut.push_back(Parton(0, 22, 0, pVec, xVec));
	    pOut[pOut.size()-1].set_form_time(0.);

	    return;
      }
    else if (Id == 21)
      {
	// gluon radiating gluon (g->gg)
	if (process == 3)
	  {
	    if (pRest/T < AMYpCut) return;

	    // sample radiated parton's momentum
	    kRest = getNewMomentumRad(pRest, T, process);
            if(kRest > pRest) return;

	    // final state parton's momentum
	    pNewRest = pRest - kRest;

	    // if pNew is smaller than pcut, final state parton is
	    // absorbed into medium
	    if (pNewRest > pcut)
	      {
		pNew = pNewRest*boostBack;
		pVecNew.Set( (px/pAbs)*pNew, (py/pAbs)*pNew, (pz/pAbs)*pNew, pNew );
		pOut.push_back(Parton(0, Id, 0, pVecNew, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
	      }

	    if (kRest > pcut)
	      {
		k = kRest*boostBack;
		kVec.Set( (px/pAbs)*k, (py/pAbs)*k, (pz/pAbs)*k, k );
		pOut.push_back(Parton(0, Id, 0, kVec, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
	      }

	    return;
	  }
	}
	// gluon split into quark-antiquark pair (g->qqbar)
	if (process == 4)
	  {
	    if (pRest/T < AMYpCut) return;

	    // sample radiated parton's momentum
	    kRest = getNewMomentumRad(pRest, T, process);
            if(kRest > pRest) return;

	    // final state parton's momentum
	    pNewRest = pRest - kRest;

	    // choose the Id of new qqbar pair. Note that we only deal with nf = 3
	    double r = ZeroOneDistribution(*GetMt19937Generator());
	    if (r < 1./3.) newId = 1;
	    else if (r < 2./3.) newId = 2;
	    else newId = 3;

	    // if pNew is smaller than pcut, final state parton is
	    // absorbed into medium
	    if (pNewRest > pcut)
	      {
		// *momentum of quark is usually larger than that of anti-quark
		pNew = pNewRest*boostBack;
		pVecNew.Set( (px/pAbs)*pNew, (py/pAbs)*pNew, (pz/pAbs)*pNew, pNew );
		pOut.push_back(Parton(0, Id, 0, pVecNew, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
	      }

	    if (kRest > pcut)
	      {
		k = kRest*boostBack;
		kVec.Set( (px/pAbs)*k, (py/pAbs)*k, (pz/pAbs)*k, k );
		pOut.push_back(Parton(0, Id, 0, kVec, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
	      }

	    return;
	  }
	// gluon scattering with either quark (gq->gq) or gluon (gg->gg)
	else if (process == 7 || process == 8)
	  {
	    //omega = getEnergyTransfer(pRest, T, process);
	    q = getMomentumTransfer(pRest, T, process);
	    pVecNewRest = getNewMomentumElas(pVecRest, q);

	    pNewRest = pVecNewRest.t();

	    // if pNew is smaller than pcut, final state parton is
	    // absorbed into medium
	    if (pNewRest > pcut)
	      {
		// Boost scattered particle to lab frame
		// 1: for brick
		if (beta < 1e-10)
		  {
		    pVecNew = pVecNewRest;
		  }
		// 2: for evolving medium
		else
		  {
		    pxNew = vx*gamma*pVecNewRest.t() 
		      + (1.+(gamma-1.)*vx*vx/(beta*beta))*pVecNewRest.x() 
		      + (gamma-1.)*vx*vy/(beta*beta)*pVecNewRest.y()
		      + (gamma-1.)*vx*vz/(beta*beta)*pVecNewRest.z();
		    pyNew = vy*gamma*pVecNewRest.t() 
		      + (1.+(gamma-1.)*vy*vy/(beta*beta))*pVecNewRest.y()
		      + (gamma-1.)*vx*vy/(beta*beta)*pVecNewRest.x()
		      + (gamma-1.)*vy*vz/(beta*beta)*pVecNewRest.z();
		    pzNew = vz*gamma*pVecNewRest.t() 
		      + (1.+(gamma-1.)*vz*vz/(beta*beta))*pVecNewRest.z()
		      + (gamma-1.)*vx*vz/(beta*beta)*pVecNewRest.x()
		      + (gamma-1.)*vy*vz/(beta*beta)*pVecNewRest.y();
			  
		    pNew = sqrt( pxNew*pxNew + pyNew*pyNew + pzNew*pzNew );
		    pVecNew.Set( pxNew, pyNew, pzNew, pNew );
		  }

		pOut.push_back(Parton(0, Id, 0, pVecNew, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
	      }

	    return;
	  }
	// gluon converting to quark
	else if (process == 11)
	  {
	    // choose the Id of new qqbar pair. Note that we only deal with nf = 3
	    double r = ZeroOneDistribution(*GetMt19937Generator());
	    if (r < 1./3.) newId = 1;
	    else if (r < 2./3.) newId = 2;
	    else newId = 3;

	    pOut.push_back(Parton(0, newId, 0, pVec, xVec));
	    pOut[pOut.size()-1].set_form_time(0.);

	    return;
	  }
      } // Id==21
  } // particle loop
}

int QinGuangYou::DetermineSign(double pRest, double T, int Id)
{
	RateElastic rateElasP, rateElasN; 
	rateElasP = getRateElasTotal(pRest, T, 1); 
	rateElasN = getRateElasTotal(pRest, T, -1); 
	double accumProb = 0.; 
	double randProb = ZeroOneDistribution(*GetMt19937Generator()); 
	double ProbP = 0.; 
	if (std::abs(Id) == 1 || std::abs(Id) == 2 || std::abs(Id) == 3)
	{
		accumProb = rateElasP.qq + rateElasP.qg + rateElasN.qq + rateElasN.qg; 
		ProbP = (rateElasP.qq + rateElasP.qg) / accumProb; 
		if (randProb > ProbP) return -1; 
		else return 1; 
	}
	if (std::abs(Id) == 21)
	{
		accumProb = rateElasP.gg + rateElasP.gq + rateElasN.gg + rateElasN.gq; 
		ProbP = (rateElasP.gg + rateElasP.gq) / accumProb; 
		if (randProb > ProbP) return -1; 
		else return 1; 
	}
	return 1; 
}
int QinGuangYou::DetermineProcess(double pRest, double T, double deltaT, int Id, int sign)
{

  double dT = deltaT/hbarc;   // time step in [GeV^(-1)]

  // get the rates for each process
  // total Probability = dT*Rate
  RateRadiative rateRad;
  rateRad = getRateRadTotal(pRest, T);
  RateElastic rateElas;
  rateElas = getRateElasTotal(pRest, T, sign);
  RateConversion rateConv;
  rateConv = getRateConv(pRest, T);

  // evolution for quark (u, d, s)
  if (std::abs(Id) == 1 || std::abs(Id) == 2 || std::abs(Id) == 3)
    {
      // multiplying by (ef/e)^2
      if (abs(Id) == 1)
	{
	  rateRad.qqgamma *= 4./9.;
	  rateConv.qgamma *= 4./9.;
	}
      else
	{
	  rateRad.qqgamma *= 1./9.;
	  rateConv.qgamma *= 1./9.;
	}

      double totalQuarkProb = 0.;

      //if (pRest > pcut) totalQuarkProb += (rateRad.qqg + rateRad.qqgamma)*dT;
      //totalQuarkProb += (rateElas.qq + rateElas.qg + rateConv.qg + rateConv.qgamma)*dT;
      if (pRest > pcut) totalQuarkProb += rateRad.qqg*dT;
      totalQuarkProb += (rateElas.qq + rateElas.qg + rateConv.qg)*dT;
      // WARN<<"omega "<<omega<<" rate_Rad "<<rateRad.qqg<<" rate_Elas "<<rateElas.qq<<" "<<rateElas.qg<<" rate_Conv "<<rateConv.qg<< "dT "<<dT; 
      // warn if total probability exceeds 1
      if (totalQuarkProb > 1.){
	WARN << " : Total Probability for quark processes exceeds 1 ("
	     << totalQuarkProb << "). "
	     << " : Most likely this means you should choose a smaller deltaT in the xml (e.g. 0.01).";
	throw std::runtime_error ("QinGuangYou probability problem.");
      }

      double accumProb = 0.;
      double nextProb = 0.;
      double Prob = 0.;

      if (ZeroOneDistribution(*GetMt19937Generator()) < totalQuarkProb)
	{
	  /* label for process
	     [1-4  : Radiation ] 1: q->qg , 2 : q->qgamma, 3 : g->gg , 4: g->qqbar
	     [5-8  : Elastic   ] 5: qq->qq, 6 : qg->qg   , 7 : gq->gq, 8: gg->gg
	     [9-11 : Conversion] 9: q->g  , 10: q->gamma , 11: g->q                */
	  double randProb = ZeroOneDistribution(*GetMt19937Generator());

	  // AMY radiation only happens if energy scale is above certain threshold.
	  // but elastic/conversion processes doesn't have threshold.
	  if(pRest > pcut)
	    {
	      Prob = rateRad.qqg*dT/totalQuarkProb;
	      if (accumProb <= randProb && randProb < (accumProb + Prob))
		return 1;

	      //accumProb += Prob;
	      //Prob = rateRad.qqgamma*dT/totalQuarkProb;
	      //if (accumProb <= randProb && randProb < (accumProb + Prob))
	      //  return 2;
	    }

	  accumProb += Prob;
	  Prob = rateElas.qq*dT/totalQuarkProb;
	  if (accumProb <= randProb && randProb < (accumProb + Prob))
	    return 5;

	  accumProb += Prob;
	  Prob = rateElas.qg*dT/totalQuarkProb;
	  if (accumProb <= randProb && randProb < (accumProb + Prob))
	    return 6;

	  accumProb += Prob;
	  Prob = rateConv.qg*dT/totalQuarkProb;
	  if (accumProb <= randProb && randProb < (accumProb + Prob))
	    return 9;

	  //accumProb += Prob;
	  //Prob = rateConv.qgamma*dT/totalQuarkProb;
	  //if (accumProb <= randProb && randProb < (accumProb + Prob))
	  //  return 10;
	}
      else
	{
	  // nothing happens to quark
	  return 0;
	}
    }
  // evolution for gluon
  else if (Id == 21)
    {
      double totalGluonProb = 0.;

      if (pRest > pcut) totalGluonProb += (rateRad.ggg + rateRad.gqq)*dT;
      totalGluonProb += (rateElas.gq + rateElas.gg + rateConv.gq)*dT;
      // WARN<<"rate_Rad "<<rateRad.ggg<<" "<<rateRad.gqq<<" rate_Elas "<<rateElas.gq<<" "<<rateElas.gg<<" rate_Conv "<<rateConv.gq; 
      // warn if total probability exceeds 1
      if (totalGluonProb > 1.)
	WARN << " : Total Probability for quark processes exceeds 1 ("
	     << totalGluonProb << ")";

      double accumProb = 0.;
      double nextProb = 0.;
      double Prob = 0.;

      if (ZeroOneDistribution(*GetMt19937Generator()) < totalGluonProb)
	{
	  /* label for process
	     [1-4  : Radiation ] 1: q->qg, 2 : q->qgamma, 3 : g->gg, 4: g->qq
	     [5-8  : Elastic   ] 5: q->q , 6 : q->g     , 7 : g->q , 8: g->g
	     [9-11 : Conversion] 9: q->g , 10: q->gamma , 11: g->q            */
	  double randProb = ZeroOneDistribution(*GetMt19937Generator());

	  // AMY radiation only happens if energy scale is above certain threshold.
	  // but elastic/conversion processes doesn't have threshold.
	  if (pRest > pcut)
	    {
	      Prob = rateRad.ggg*dT/totalGluonProb;
	      if (accumProb <= randProb && randProb < (accumProb + Prob))
		return 3;

	      accumProb += Prob;
	      Prob = rateRad.gqq*dT/totalGluonProb;
	      if (accumProb <= randProb && randProb < (accumProb + Prob))
		return 4;
	    }

	  accumProb += Prob;
	  Prob = rateElas.gq*dT/totalGluonProb;
	  if (accumProb <= randProb && randProb < (accumProb + Prob))
	    return 7;

	  accumProb += Prob;
	  Prob = rateElas.gg*dT/totalGluonProb;
	  if (accumProb <= randProb && randProb < (accumProb + Prob))
	    return 8;

	  accumProb += Prob;
	  Prob = rateConv.gq*dT/totalGluonProb;
	  if (accumProb <= randProb && randProb < (accumProb + Prob))
	    return 11;
	}
      else
	{
	  // nothing happens to gluon
	  return 0;
	}
    }

  // if parton is other than u,d,s,g, do nothing
  return 0;
}

RateRadiative QinGuangYou::getRateRadTotal(double pRest, double T)
{
  RateRadiative rate;
/*
  double u = pRest/T;  // making arguments in log to be dimensionless

  rate.qqg = (0.8616 - 3.2913/(u*u) + 2.1102/u - 0.9485/sqrt(u))*pow(g, 4.)*T;
  rate.ggg = (1.9463 + 61.7856/(u*u*u) - 30.7877/(u*u) + 8.0409/u - 2.6249/sqrt(u))
    *pow(g, 4.)*T;
  rate.gqq = (2.5830/(u*u*u) - 1.7010/(u*u) + 1.4977/u - 1.1961/pow(u,0.8) + 0.1807/sqrt(u))
    *pow(g, 4.)*T*nf; 
  rate.qqgamma  = (0.0053056 + 2.3279/pow(u,3.) - 0.6676/u + 0.3223/sqrt(u))
    *pow(g, 4.)*alpha_em/alpha_s*T;

  double runningFactor = log(g*T*pow(10., 0.25)/.175)/log(g*T*pow(u, 0.25)/.175);
  if (runningFactor < 1.)
    {
      rate.qqg *= runningFactor;
      rate.gqq *= runningFactor;
      rate.ggg *= runningFactor;
      rate.qqgamma *= runningFactor;
    }
*/
  rate.qqg = 0.; 
  rate.gqq = 0.; 
  rate.ggg = 0.; 
  rate.qqgamma = 0.; 
  return rate;
}

RateRadiative QinGuangYou::getRateRadPos(double u, double T)
{
  RateRadiative rate;
/*
  rate.qqg = (0.5322 - 3.1037/(u*u) + 2.0139/u - 0.9417/sqrt(u))*pow(g, 4.)*T;
  rate.ggg = (1.1923 - 11.5250/(u*u*u) + 3.3010/u - 1.9049/sqrt(u))*pow(g, 4.)*T;
  rate.gqq = (0.0004656 - 0.04621/(u*u) + 0.0999/u - 0.08171/pow(u,0.8) 
              + 0.008090/pow(u,0.2) - 1.2525*pow(10.,-8.)*u)*pow(g, 4.)*T*nf; 
  rate.qqgamma = 0.;
*/
  rate.qqg = 0.; 
  rate.ggg = 0.; 
  rate.gqq = 0.; 
  rate.qqgamma = 0.; 
  return rate;
}

RateRadiative QinGuangYou::getRateRadNeg(double u, double T)
{
  RateRadiative rate;
/*
  rate.qqg = (0.3292 - 0.6759/(u*u) + 0.4871/pow(u,1.5) - 0.05393/u + 0.007878/sqrt(u))
    *pow(g, 4.)*T;
  rate.ggg = (0.7409 + 1.8608/(u*u*u) - 0.1353/(u*u) + 0.1401/u)*pow(g, 4.)*T;
  rate.gqq = (0.03215/(u*u*u) + 0.01419/(u*u) + 0.004338/u - 0.00001246/sqrt(u))
    *pow(g, 4.)*T*nf;
  rate.qqgamma = 0.;
*/
  rate.qqg = 0.; 
  rate.ggg = 0.; 
  rate.gqq = 0.; 
  rate.qqgamma = 0.; 
  return rate;
}

double QinGuangYou::getNewMomentumRad(double pRest, double T, int process)
{
  double newp = 0.;
  double randA;
  double x, y;
  double fy, fyAct;

  RateRadiative Pos, Neg;
  double u = pRest/T;  // making arguments in log to be dimensionless

  Pos = getRateRadPos(u, T);
  Neg = getRateRadNeg(u, T);

  // this switch will hold the decision whether k is positive or negative:
  // 0 : negative, 1 : positive
  int posNegSwitch = 1; 

  /* process == 1 : quark radiating gluon
     process == 2 : quark radiating photon
     process == 3 : gluon radiating gluon
     process == 4 : gluon split into quark-antiquark pair */

  if (process == 1)
    {
      // decide whether k shall be positive or negative 
      // if x (uniform on [0,1]) < area(k<0)/area(all k) then k < 0
      if (ZeroOneDistribution(*GetMt19937Generator()) < Neg.qqg/(Neg.qqg+Pos.qqg))
	posNegSwitch = 0;

      if (posNegSwitch == 1) // if k > 0
	{
	  do
	    {
	      //randA is a uniform random number on [0, Area under the envelop function]
	      randA = ZeroOneDistribution(*GetMt19937Generator())*area(u+12., u, posNegSwitch, 1);
	      y = 2.5/(LambertWM(2.59235*pow(10.,23.)*exp(-100.*randA)));

	      fy = 0.025/(y*y)+0.01/y;             // total area under the envelop function
	      fyAct = function(u, y, process);     // actual rate

	      x = ZeroOneDistribution(*GetMt19937Generator());       // random number, uniform on [0,1]

	    } while (x > fyAct/fy); 
	  // reject if x is larger than the ratio fyAct/fy
	  newp = y;
	}
      else // if k < 0
	{
	  do
	    {
	      //randA is a uniform random number on [0, Area under the envelop function]
	      randA = ZeroOneDistribution(*GetMt19937Generator())*area(-0.05, u, posNegSwitch, 1);
	      y = -12./(1.+480.*randA);

	      fy = 0.025/(y*y);                    // total area under the envelop function
	      fyAct = function(u, y, process);     // actual rate

	      x = ZeroOneDistribution(*GetMt19937Generator());       // random number, uniform on [0,1]

	    } while (x > fyAct/fy); 
	  // reject if x is larger than the ratio fyAct/fy
	  newp = y;
	}
    }
  else if (process == 2)
    {
      do
	{
	  randA = ZeroOneDistribution(*GetMt19937Generator())*area(1.15*u, u, posNegSwitch, 2);
	  y = 83895.3*pow(pow(u, 0.5)*randA, 10./3.);

	  fy = (0.01/(pow(y, 0.7)))/pow(u, 0.5); // total area under the envelop function
	  fyAct = function(u, y, process);       // actual rate

	  x = ZeroOneDistribution(*GetMt19937Generator());         // random number, uniform on [0,1]

	} while (x > fyAct/fy); 
      // reject if x is larger than the ratio fyAct/fy
      newp = y;
    }
  else if (process == 3)
    {
      // decide whether k shall be positive or negative 
      // if x (uniform on [0,1]) < area(k<0)/area(all k) then k < 0
      if (ZeroOneDistribution(*GetMt19937Generator()) < Neg.ggg/(Neg.ggg+Pos.ggg))
	posNegSwitch = 0;

      if( posNegSwitch == 1 ) // if k > 0
	{
	  do
	    {
	      //randA is a uniform random number on [0, Area under the envelop function]
	      randA = ZeroOneDistribution(*GetMt19937Generator())*area(u/2., u, posNegSwitch, 3);
	      y = 5./(LambertWM(2.68812*pow(10., 45.)*exp(-50.*randA)));

	      fy = 0.1/(y*y)+0.02/y;               // total area under the envelop function
	      fyAct = function(u, y, process);     // actual rate

	      x = ZeroOneDistribution(*GetMt19937Generator());       // random number, uniform on [0,1]

	    } while (x > fyAct/fy); 
	  // reject if x is larger than the ratio fyAct/fy
	  newp = y;
	}
      else // if k < 0
	{
	  do
	    {
	      //randA is a uniform random number on [0, Area under the envelop function]
	      randA = ZeroOneDistribution(*GetMt19937Generator())*area(-0.05, u, posNegSwitch, 3);
	      y = -12./(1. + 120.*randA);

	      fy = 0.1/(y*y);                      // total area under the envelop function
	      fyAct = function(u, y, process);     // actual rate

	      x = ZeroOneDistribution(*GetMt19937Generator());       // random number, uniform on [0,1]

	    } while(x > fyAct/fy); 
	  // reject if x is larger than the ratio fyAct/fy
	  newp = y;
	}
    }
  else if (process == 4)
    {
      // decide whether k shall be positive or negative 
      // if x (uniform on [0,1]) < area(k<0)/area(all k) then k < 0
      if (ZeroOneDistribution(*GetMt19937Generator()) < Neg.gqq/(Neg.gqq+Pos.gqq))
	posNegSwitch = 0;

      if( posNegSwitch == 1 ) // if k > 0
	{
	  do
	    {
	      //randA is a uniform random number on [0, Area under the envelop function]
	      randA = ZeroOneDistribution(*GetMt19937Generator())*area(u/2., u, posNegSwitch, 4);
	      y = 0.83333*(0.06*function(u, 0.05, process)+randA)/function(u, 0.05, process);

	      fy = 1.2*function(u, 0.05, process); // total area under the envelop function
	      fyAct = function(u, y, process);     // actual rate

	      x = ZeroOneDistribution(*GetMt19937Generator());       // random number, uniform on [0,1]

	    } while (x > fyAct/fy); 
	  // reject if x is larger than the ratio fyAct/fy
	  newp = y;
	}
      else // if k < 0
	{
	  do
	    {
	      //randA is a uniform random number on [0, Area under the envelop function]
	      randA = ZeroOneDistribution(*GetMt19937Generator())*area(-0.05, u, posNegSwitch, 4);
	      y = (2.5-u*log(7.81082*pow(10., -6.)*exp(14.5/u)+(-115.883+113.566*u)*randA))/(1.-0.98*u);

	      fy = 0.98*exp((1.-1./u)*(-2.5+y))/u; // total area under the envelop function
	      fyAct = function(u, y, process);     // actual rate

	      x = ZeroOneDistribution(*GetMt19937Generator());       // random number, uniform on [0,1]

	    } while (x > fyAct/fy); 
	  // reject if x is larger than the ratio fyAct/fy
	  newp = y;
	}
    }
  else
    {
      WARN << "Invalid process number (" << process << ")";
    }

  return newp;
}

// calculates the area under the envelope function when using the rejection method
// integrals had been solved analytically before
double QinGuangYou::area(double y, double u, int posNegSwitch, int process)
{
  if (process == 1)
    {
      if (posNegSwitch == 1)
	return (0.5299 - 0.025/y + 0.01*log(y));
      else 
	return (-0.002083-0.025/y);
    }
  else if (process == 2)
    {
      return ((0.03333*pow(y,0.3))/pow(u,0.5));
    }
  else if (process == 3)
    {
      if (posNegSwitch == 1)
	return (2.05991 - 0.1/y + 0.02*log(y));
      else 
	return (-0.008333 - 0.1/y);
    }      
  else if (process == 4)
    {
      if (posNegSwitch == 1)
	return (1.2*function(u, 0.05, process)*(y-0.05));
      else 
	return ((6.8778*pow(10., -8.)*exp(14.5/u)
		 -0.008805*exp((2.5-y+0.98*u*y)/u))/(1.0204-u));
    }      

  return 0.;
}

double QinGuangYou::function(double u, double y, int process)
{
  if (process == 1)      return getRate_qqg(u, y);
  else if (process == 2) return getRate_qqgamma(u, y);
  else if (process == 3) return getRate_ggg(u, y);
  else if (process == 4) return getRate_gqq(u, y);

  return 0.;
}

RateElastic QinGuangYou::getRateElasTotal(double pRest, double T, int pn)
{
    RateElastic rate;
    double domega = std::abs(omega); 
    double fb = 1./(exp(domega/T)-1.); 
    if (pn > 0)
    {
    	rate.qq = (1+fb)/domega*getRateElas(pRest, T).qq; 
    	rate.qg = (1+fb)/domega*getRateElas(pRest, T).qg; 
    	rate.gq = (1+fb)/domega*getRateElas(pRest, T).gq; 
    	rate.gg = (1+fb)/domega*getRateElas(pRest, T).gg; 
    	//WARN<<"THE SIGN IS "<<pn<<" rate "<<rate.qq<<" "<<rate.qg; 
    }
    else
    {
    	rate.qq = fb/domega*getRateElas(pRest, T).qq; 
    	rate.qg = fb/domega*getRateElas(pRest, T).qg; 
    	rate.gq = fb/domega*getRateElas(pRest, T).gq; 
    	rate.gg = fb/domega*getRateElas(pRest, T).gg; 
    	//WARN<<"THE SIGN IS "<<pn<<" rate "<<rate.qq<<" "<<rate.qg; 
    }
    return rate;
}

RateElastic QinGuangYou::getRateElas(double pRest, double T)
{
  RateElastic rate; 
  double cb = -0.577215664901532 - 0.569961; 
  double cf = cb + log(2.); 
  double cs = -1.66246; 
  double mg2 = 2.*M_PI*alpha_s*T*T*(1+nf/6.);  
  rate.qq = 2./9.*nf*M_PI*alpha_s*alpha_s*T*T*(log(pRest*T/mg2)+cf+23./12.+cs); 
  rate.qg = 4./3.*M_PI*alpha_s*alpha_s*T*T*(log(pRest*T/mg2)+cb+13./6.+cs); 
  rate.gq = 1./2.*nf*M_PI*alpha_s*alpha_s*T*T*(log(pRest*T/mg2)+cf+13./6.+cs); 
  rate.gg = 3.*M_PI*alpha_s*alpha_s*T*T*(log(pRest*T/mg2)+cb+131./48.+cs); 
  //rate.qgamma = 2.*M_PI*alpha_s*alpha_em*T*T/(3.*p)*(0.5*log(p*T/((1./6.)*pow(g*T, 2.)))-0.36149);
  return rate;
}

RateConversion QinGuangYou::getRateConv(double p, double T)
{
  RateConversion rate;
/*
  rate.qg     = 4./3.*2.*M_PI*alpha_s*alpha_s*T*T/(3.*p)
    *(0.5*log(p*T/((1./6.)*pow(g*T, 2.)))-0.36149);
  rate.gq     = nf*3./8.*4./3.*2.*M_PI*alpha_s*alpha_s*T*T/(3.*p)
    *(0.5*log(p*T/((1./6.)*pow(g*T, 2.)))-0.36149);
  rate.qgamma = 2.*M_PI*alpha_s*alpha_em*T*T/(3.*p)
    *(0.5*log(p*T/((1./6.)*pow(g*T, 2.)))-0.36149);
*/
  rate.qg = 0.; 
  rate.gq = 0.; 
  rate.qgamma = 0.; 
  return rate;
}
double QinGuangYou::getMomentumTransfer(double pRest, double T, int process)
{
/*  double q = 0.;
  double randA;
  double A, B;
  double x, y;
  double fy, fyAct;

  double u = pRest/T;  // making arguments in log to be dimensionless
  double omegaT = omega / T;
*/
  /* process == 5 : qq -> qq (first q is hard, second q is soft) 
     process == 6 : qg -> qg 
     process == 7 : gq -> gq 
     process == 8 : gg -> gg  */

  // small omega using the rejection method
/*  if (omegaT < 10. && omegaT > -3.)
    {
      if (process == 5 || process == 7) // for qq or gq
	{
	  A = (0.7+alpha_s)*0.0014*(1000.+40./sqrt(omegaT*omegaT)+10.5*pow(omegaT, 4.))*alpha_s;
	  B = 2.*sqrt(omegaT*omegaT)+0.01;
	}
      else if (process == 6 || process == 8) // for qg or gg
	{
	  A = (0.7+alpha_s)*0.0022*(1000.+40./sqrt(omegaT*omegaT)+10.5* pow(omegaT, 4.))*alpha_s;
	  B = 2.*sqrt(omegaT*omegaT)+0.002;
	}
      else
	{
	  WARN << "Invalid process number (" << process << ")";

	  A = 0.;
	  B = 0.;
	}

      do
	{
	  //randA is a uniform random number on [0, Area under the envelop function]
	  randA = ZeroOneDistribution(*GetMt19937Generator())*areaQ(u, process, T);
	  y = pow(B, 0.25)*sqrt(tan((2.*sqrt(B)*randA+A*atan(omegaT*omegaT/sqrt(B)))/A));

	  fy = A*y/(pow(y, 4.)+B);                       // total area under the envelop function
	  fyAct = functionQ(u, y, process, T);  // actual rate

	  x = ZeroOneDistribution(*GetMt19937Generator());                 // random number, uniform on [0,1]

	} while (x > fyAct/fy); 
      // reject if x is larger than the ratio fyAct/fy
      q = y;
    }
  // large omega using the Metropolis method
  else
    {
      double g = 0, g_new = 0;
      double ratio;
      double y_new;
    
      // the ranges in which the variables u and phi need to be sampled
      const double y_min = sqrt(omegaT*omegaT);
      const double y_max = u;

      // randomly select initial values of q=y, such that
      do
	{
	  y = y_min+ZeroOneDistribution(*GetMt19937Generator())*(y_max-y_min);
	  g = functionQ(u, y, process, T);

	} while (g == 0. || g != g);
      // number of steps in the Markov chain
      const int n_steps = 500;
    
      for(int i=0; i<n_steps; i++)
	{        
	  do
	    {
	      y_new = y_min+ZeroOneDistribution(*GetMt19937Generator())*(y_max-y_min); 
	    }
	  while (y_new < y_min || y_new > y_max);                      
	  // check that the new value is in range

	  g_new = functionQ(u, y_new, process, T); // calculate the function at the proposed point
	  if (g_new != g_new) g_new = 0.; 
	  ratio = g_new/g;                             // ratio of g(y_new)/g(y)

	  // accept if probability g(y_new)/g(y) is larger than randon number
	  if (ZeroOneDistribution(*GetMt19937Generator()) < ratio)
	    {
	      y = y_new;
	      g = g_new;
	    }
	}
      q = y;
      //WARN<<" u "<<u<<" omega "<<omega<<" q "<<q<<" process "<<process<<" functionQ "<<functionQ(u, omega, q, process); 
    }
  // WARN<<"omega "<<omega<<" q "<<q; 
  return q*T;  // q*T is in [GeV]
*/
	return omega; 
}

// calculates the area under the envelope function when using the rejection method
// integrals had been solved analytically before
double QinGuangYou::areaOmega(double u, int posNegSwitch, int process)
{
  if (process == 5 || process == 7)
    {
      if (posNegSwitch == 1)
	return (0.0333333*alpha_s*(7.+ 4.*alpha_s)*(2.99573+log(u)));
      else 
	return (-0.133333*alpha_s*(1.75+alpha_s)*log(-0.0833333*u));
    }
  else if (process == 6 || process == 8)
    {
      if (posNegSwitch == 1)
	return (0.05*alpha_s*(7.+ 4.*alpha_s)*(2.99573+log(u)));
      else 
	return (-0.2*alpha_s*(1.75+alpha_s)*log(-0.0833333*u));
    }
  else
    {
      WARN << "Invalid process number (" << process << ")";
    }

  return 0.;
}

double QinGuangYou::getElasticRateQ(double u, double q, int process, double T)
{
  double omegaT = omega / T; 
  if (q > std::abs(omegaT) && ((omegaT > 0 && omegaT < u) || (omegaT < 0 && omegaT > -u)))
    return use_elastic_table_q(q, process, T);

  return 0.;
}

double QinGuangYou::areaQ(double u, int process, double T)
{
  double A, B;
  double areaQ;
  double omegaT = omega / T; 

  if (process == 5 || process == 7)
    {
      A = (0.7+alpha_s-0.3)*0.0014*alpha_s
        *(1000.+40./sqrt(omegaT*omegaT)+10.5*pow(omegaT, 4.));
      B = 2.*sqrt(omegaT*omegaT)+0.01;
    }
  else if (process == 6 || process == 8)
    {
      A = (0.7+alpha_s-0.3)*0.0022*alpha_s
        *(1000.+40./sqrt(omegaT*omegaT)+10.5*pow(omegaT, 4.));
      B = 2.*sqrt(omegaT*omegaT)+0.002;
    }
  else
    {
      WARN << "Invalid process number (" << process << ")";

      A = 0.;
      B = 0.;
    }
  areaQ = (0.5*A*(atan(u*u/sqrt(B))-atan(omegaT*omegaT/sqrt(B))))/sqrt(B);

  return areaQ;
}

FourVector QinGuangYou::getNewMomentumElas(FourVector pVec, double q)
{
  FourVector pVecNew, pVecNewTemp;
  FourVector etVec, qtVec, qlVec;
  FourVector r;
  double qt, ql;
  double cosTheta_pq;
  double pAbs=pVec.t();
  double phi;
  double M[3][3]; //rotation matrix
  double u;
  double xx, yy, zz ,tt;

  if (omega == q)
    {
      xx = pVec.x()*(pAbs-omega)/pAbs;
      yy = pVec.y()*(pAbs-omega)/pAbs;
      zz = pVec.z()*(pAbs-omega)/pAbs;
      tt = pVec.t()*(pAbs-omega)/pAbs;
    
      pVecNew.Set(xx, yy, zz ,tt);
      // WARN<<"energy "<<pVecNew.t(); 
      return pVecNew;
    }

  cosTheta_pq = (-omega*omega+2.*pAbs*omega+q*q)/(2.*pAbs*q);
  qt = q*sqrt(1.-cosTheta_pq*cosTheta_pq);          // transverse momentum transfer
  ql = q*cosTheta_pq;
  
  if (pVec.y()*pVec.y() > pVec.x()*pVec.x())
    {
      xx = 0.;
      yy = -pVec.z()/sqrt(pVec.y()*pVec.y()+pVec.z()*pVec.z());
      zz = pVec.y()/sqrt(pVec.y()*pVec.y()+pVec.z()*pVec.z());
      tt = sqrt(yy*yy+zz*zz);

      etVec.Set(xx, yy, zz, tt);
    }
  else
    {
      xx = pVec.z()/sqrt(pVec.x()*pVec.x()+pVec.z()*pVec.z());
      yy = 0.;
      zz = -pVec.x()/sqrt(pVec.x()*pVec.x()+pVec.z()*pVec.z());
      tt = sqrt(xx*xx+zz*zz);

      etVec.Set(xx, yy, zz, tt);
    }

  // the transverse transferred momentum vector
  qtVec.Set(etVec.x()*qt, etVec.y()*qt, etVec.z()*qt, etVec.t()*qt);
  // the longuitudinal transferred momentum vector
  qlVec.Set(pVec.x()/pAbs*ql, pVec.y()/pAbs*ql,
            pVec.z()/pAbs*ql, pVec.t()/pAbs*ql);

  pVecNewTemp = pVec;
  pVecNewTemp -= qtVec;  // change transverse momentum
  pVecNewTemp -= qlVec;  // change longitudinal momentum
  
  // WARN<<"etVec.t "<<etVec.t(); 

  phi = 2.*M_PI*ZeroOneDistribution(*GetMt19937Generator());
  r.Set(pVec.x()/pVec.t(), pVec.y()/pVec.t(), pVec.z()/pVec.t(), 1.);
  u = 1.-cos(phi);
 
  // define the rotation matrix for rotations around pvecRest
  M[0][0]=r.x()*r.x()*u+cos(phi);
  M[1][0]=r.x()*r.y()*u-r.z()*sin(phi);
  M[2][0]=r.x()*r.z()*u+r.y()*sin(phi);

  M[0][1]=r.y()*r.x()*u+r.z()*sin(phi);
  M[1][1]=r.y()*r.y()*u+cos(phi);
  M[2][1]=r.y()*r.z()*u-r.x()*sin(phi);

  M[0][2]=r.z()*r.x()*u-r.y()*sin(phi);
  M[1][2]=r.z()*r.y()*u+r.x()*sin(phi);
  M[2][2]=r.z()*r.z()*u+cos(phi);

  xx = M[0][0]*pVecNewTemp.x()+M[0][1]*pVecNewTemp.y()+M[0][2]*pVecNewTemp.z();
  yy = M[1][0]*pVecNewTemp.x()+M[1][1]*pVecNewTemp.y()+M[1][2]*pVecNewTemp.z();
  zz = M[2][0]*pVecNewTemp.x()+M[2][1]*pVecNewTemp.y()+M[2][2]*pVecNewTemp.z();
  tt = sqrt(xx*xx+yy*yy+zz*zz);

  pVecNew.Set(xx, yy, zz, tt);
  // WARN<<"energy transfer "<<omega<<" "<<sqrt(qt*qt+ql*ql)<<" "<<sqrt(qtVec.t()*qtVec.t()+qlVec.t()*qlVec.t())<<" "<<pVecNew.t(); 
  return pVecNew;
}

// Reads in the binary stored file of dGamma values
void QinGuangYou::readRadiativeRate(Gamma_info *dat, dGammas *Gam)
{
  FILE *rfile;
  string filename;
  filename = PathToTables+"radgamma";

  cout << "Reading rates of inelastic collisions from file " << endl;
  cout << filename.c_str() << " ... " << endl;
  size_t bytes_read;

  rfile = fopen(filename.c_str(), "rb"); 
  bytes_read = fread((char *)(&dat->ddf), sizeof(double), 1, rfile);
  bytes_read = fread((char *)(&dat->dda), sizeof(double), 1, rfile);
  bytes_read = fread((char *)(&dat->dcf), sizeof(double), 1, rfile);
  bytes_read = fread((char *)(&dat->dca), sizeof(double), 1, rfile);
  bytes_read = fread((char *)(&dat->Nc), sizeof(int), 1, rfile);
  bytes_read = fread((char *)(&dat->Nf), sizeof(int), 1, rfile);
  bytes_read = fread((char *)(&dat->BetheHeitler),sizeof(int) , 1, rfile);
  bytes_read = fread((char *)(&dat->BDMPS), sizeof(int), 1, rfile);
  bytes_read = fread((char *)(&dat->include_gluons), sizeof(int), 1, rfile);
  bytes_read = fread((char *)Gam->qqg, sizeof(double), NP*NK, rfile);
  bytes_read = fread((char *)Gam->tau_qqg, sizeof(double), NP*NK, rfile);
  bytes_read = fread((char *)Gam->gqq, sizeof(double), NP*NK , rfile);
  bytes_read = fread((char *)Gam->tau_gqq, sizeof(double), NP*NK, rfile);
  bytes_read = fread((char *)Gam->ggg, sizeof(double), NP*NK, rfile);
  bytes_read = fread((char *)Gam->tau_ggg, sizeof(double), NP*NK, rfile);
  bytes_read = fread((char *)Gam->qqgamma, sizeof(double), NP*NK, rfile);
  bytes_read = fread((char *)Gam->tau_qqgamma, sizeof(double), NP*NK, rfile);
  fclose (rfile);
  cout << " ok." << endl;

  dat->Nf = nf;
  dat->dp = 0.05;
  dat->p_max = 20;
  dat->p_min = 0;                 // exp(LogEmin) set to zero because my array starts at zero!...
  dat->n_p = static_cast<int>(1.001+dat->p_max/dat->dp);    // np = int(0.4 + 121 / 0.5) = 241
  dat->p_max = dat->dp*dat->n_p;                            // p_max = 0.5*241 = 120.5
  dat->n_pmin = static_cast<int>(1.001+dat->p_min/dat->dp); // np_min = int(0.4 + 3.3 / 0.5) = 7
  dat->n_p -= dat->n_pmin-1;                                // n_p = 241 - (7 - 1) = 235
  dat->p_min = dat->dp * dat->n_pmin;                       // p_min = 0.5 * 7 = 3.5
  dat->n_kmin = 1+2*(static_cast<int>(2./dat->dp));
  dat->k_min = -dat->dp*dat->n_kmin;
  dat->n_k = static_cast<int>((8+dat->p_max)/(2*dat->dp));
  dat->k_max = 2*dat->dp*(dat->n_k-1)+dat->k_min;
}

void QinGuangYou::readElasticRateOmega()
{
  ifstream fin;
  string filename[2];

  double as, iomega;
  double dGamma;
      
  // open files with data to read in:
  filename[0] = PathToTables + "logEnDtrqq";
  filename[1] = PathToTables + "logEnDtrqg";
  
  cout << "Reading rates of elastic collisions from files" << endl;
  cout << filename[0] << endl;
  cout << filename[1] << " ..." << endl;

  fin.open(filename[0].c_str(), ios::in);
  if(!fin)
    {
      cerr << "[readElasticRateOmega]: ERROR: Unable to open file " << filename[0] << endl;
      exit(1);
    }

  int ik = 0;
  while (!fin.eof())
    {
      fin >> as;
      fin >> iomega;
      fin >> dGamma;
      dGamma_qq->push_back(dGamma);

      ik++;
    }
  fin.close();

  fin.open(filename[1].c_str(), ios::in);
  if(!fin)
    {
      cerr << "[readElasticRateOmega]: ERROR: Unable to open file " << filename[1] << endl;
      exit(1);
    }

  ik = 0;
  while (!fin.eof())
    {
      fin >> as;
      fin >> iomega;
      fin >> dGamma;
      dGamma_qg->push_back(dGamma);

      ik++;
    }
  fin.close();

  cout << " ok." << endl;
}

void QinGuangYou::readElasticRateQ()
{
  ifstream fin;
  string filename[2];

  double as, iomega, q;
  double dGamma;
      
  // open files with data to read in:
  filename[0] = PathToTables + "logEnDqtrqq";
  filename[1] = PathToTables + "logEnDqtrqg";
  
  cout << "Reading rates of elastic collisions from files" << endl;
  cout << filename[0] << endl;
  cout << filename[1] << " ..." << endl;

  fin.open(filename[0].c_str(), ios::in);
  if(!fin)
    {
      cerr << "[readElasticRateQ]: ERROR: Unable to open file " << filename[0] << endl;
      exit(1);
    }

  int ik = 0;
  while (!fin.eof())
    {
      fin >> as;
      fin >> iomega;
      fin >> q;
      fin >> dGamma;
      dGamma_qq_q->push_back(dGamma);

      ik++;
    }
  fin.close();

  fin.open(filename[1].c_str(),ios::in);
  if(!fin)
    {
      cerr << "[readElasticRateQ]: ERROR: Unable to open file " << filename[1] << endl;
      exit(1);
    }
  
  ik = 0;
  while (!fin.eof())
    {
      fin >> as;
      fin >> iomega;
      fin >> q;
      fin >> dGamma;
      dGamma_qg_q->push_back(dGamma);

      ik++;
    }
  fin.close();
  
  cout << " ok." << endl;
}

double QinGuangYou::getRate_qqg(double p, double k)
{
  return use_table(p, k, Gam.qqg, 0);
}

double QinGuangYou::getRate_gqq(double p, double k)
{
  if (k < p/2.) return use_table(p, k, Gam.gqq, 1);
  else return 0.;
}

double QinGuangYou::getRate_ggg(double p, double k)
{
  if (k < p/2.) return use_table(p, k, Gam.ggg, 2);
  else return 0.;
}

double QinGuangYou::getRate_qqgamma(double p, double k)
{
  return use_table(p, k, Gam.qqgamma, 3);
}

double QinGuangYou::use_table(double p, double k, double dGamma[NP][NK], int which_kind)
/* Uses the lookup table and simple interpolation to get the value
   of dGamma/dkdx at some value of p,k.
   This works by inverting the relations between (p,k) and (n_p,n_k)
   used in building the table, to find out what continuous values
   of n_p, n_k should be considered; then linearly interpolates.     */
{
  double a, b, result;     // fraction of way from corner of box
  int    n_p, n_k;         // location of corner of box

  // out of range
  if ((p < 4.01) || (p > 46000.) || (k < -12.) || (k > p+12.))
    return 0.;

  if ((which_kind % 3) && (k > p/2))
    k = p - k;  // Take advantage of symmetry in these cases

  a = 24.7743737154026*log(p*0.2493765586034912718l);
  n_p = (int)a;
  a -= n_p;
  if (k < 2.)
    {
      if (k < -1)
	{
	  if (k < -2) b = 60.+5.*k;
	  else b = 70.+10.*k;
	}
      else
	{
	  if (k < 1.) b = 80. + 20.*k;
	  else b = 90.+10.*k;
	}
    }
  else if ( k < p-2. )
    { /* This is that tricky middle ground. */
      b = 190.-10.*log(1.000670700260932956l/ 
		       (0.0003353501304664781l+(k-2.)/(p-4.))-1.);
    }
  else
    {
      if (k < p+1.)
	{
	  if (k < p-1.) b = 290. + 10.*(k-p);
	  else  b = 300. + 20.*(k-p);
	}
      else
	{
	  if (k < p+2.) b = 310. + 10.*(k-p);
	  else b = 320. + 5.*(k-p);
	}
    }

  n_k = (int)b;
  b -= n_k;
  result = (1.-a)*((1.-b)*dGamma[n_p][n_k]+b*dGamma[n_p][n_k+1])
    +a*((1.-b)*dGamma[n_p+1][n_k]+b*dGamma[n_p+1][n_k+1]);

  if (std::abs(k) > 0.001) // Avoid division by 0, should never get asked for
    {
      switch (which_kind)
	{
	case 0:
	  result /= k;
	  if (k < 20.)
	    result /= 1.-exp(-k);
	  if (k > p-20.)
	    result /= 1. + exp(k-p);
	  break;
	case 1:
	  result /= p;
	  if (k < 20.)
	    result /= 1 + exp(-k);
	  if (k > p-20.)
	    result /= 1. + exp(k-p);
	  break;
	case 2:
	  result /= k*(p-k)/p;
	  if (k < 20.)
	    result /= 1.-exp(-k);
	  if (k > p-20.)
	    result /= 1.-exp(k-p);
	  break;
	case 3:
	  result /= k;
	  if (k < 0) result = 0.;
	  if (k > p-20.)
	    result /= 1. + exp(k-p);
	  break;
	}
    }

  return result;
}

double QinGuangYou::use_elastic_table_q(double q, int which_kind, double T)
/* Uses the lookup table and simple interpolation to get the value
   of dGamma/domegadq at some value of omega, q, and alpha_s. */
{
  double result;
  double alphaFrac, omegaFrac, qFrac;
  int iOmega, iAlphas, iQ;
  int position, positionAlphaUp, positionOmegaUp, positionAlphaUpOmegaUp;
  int positionQUp, positionAlphaUpQUp, positionOmegaUpQUp, positionAlphaUpOmegaUpQUp;
  int position2QUp;
  double rate, rateAlphaUp, rateOmegaUp, rateAlphaUpOmegaUp;
  double rateQUp, rateAlphaUpQUp, rateOmegaUpQUp, rateAlphaUpOmegaUpQUp;
  double rate2QUp, rateAlphaUp2QUp, rateOmegaUp2QUp, rateAlphaUpOmegaUp2QUp;
  double rateOmegaAv, rateAlphaUpOmegaAv, rateQUpOmegaAv, rateAlphaUpQUpOmegaAv;
  double rate2QUpOmegaAv, rateAlphaUp2QUpOmegaAv;
  double rateQAv, rateAlphaUpQAv;
  double slope, slopeAlphaUp;
  double omegaT = omega / T; 

  rate2QUp = 0.;
  rateAlphaUp2QUp = 0.;
  rateOmegaUp2QUp = 0.;
  rateAlphaUpOmegaUp2QUp = 0.;
  rateOmegaAv = 0.;
  rateAlphaUpOmegaAv = 0.;
  rateQUpOmegaAv = 0.;
  rateAlphaUpQUpOmegaAv = 0.;
  rate2QUpOmegaAv = 0.;
  rateAlphaUp2QUpOmegaAv = 0.;

  if (omegaT > 0.) iOmega = Nomega/2+floor((log(omegaT)+5)/omegaStep);
  else iOmega = Nomega/2-ceil((log(-omegaT)+5)/omegaStep)-1;
  iQ = floor((log(q)+5)/qStep+0.0001);
  iAlphas = floor((alpha_s-0.15)/alphaStep+0.0001);

  position = iQ + Nq*(iOmega+Nomega*(iAlphas));
  positionAlphaUp = iQ + Nq*(iOmega+Nomega*(iAlphas+1));
  positionOmegaUp = iQ + Nq*(iOmega+1+Nomega*(iAlphas));
  positionQUp = iQ+1 + Nq*(iOmega+Nomega*(iAlphas));
  position2QUp = iQ+2 + Nq*(iOmega+Nomega*(iAlphas));
  positionAlphaUpOmegaUp = iQ + Nq*(iOmega+1+Nomega*(iAlphas+1));
  positionAlphaUpQUp = iQ+1 + Nq*(iOmega+Nomega*(iAlphas+1));
  positionOmegaUpQUp = iQ+1 + Nq*(iOmega+1+Nomega*(iAlphas));
  positionAlphaUpOmegaUpQUp = iQ+1 + Nq*(iOmega+1+Nomega*(iAlphas+1));

  alphaFrac = (alpha_s-(floor((alpha_s-alphaMin)/alphaStep)*alphaStep+alphaMin))/alphaStep;
  if (omegaT > 0.) 
    {
      if (exp(ceil((log(omegaT)+5)/omegaStep)*omegaStep-5)
	  != exp(floor((log(omegaT)+5)/omegaStep)*omegaStep-5))
	omegaFrac = (omegaT-(exp(floor((log(omegaT)+5)/omegaStep)*omegaStep-5)))
	  /((exp(ceil((log(omegaT)+5)/omegaStep)*omegaStep-5))
	    -exp(floor((log(omegaT)+5)/omegaStep)*omegaStep-5));
      else omegaFrac = 0.;
    }
  else
    {
      if (exp(ceil((log(-omegaT)+5)/omegaStep)*omegaStep-5)
	  != exp(floor((log(-omegaT)+5)/omegaStep)*omegaStep-5))
	omegaFrac = (-omegaT-(exp(floor((log(-omegaT)+5)/omegaStep)*omegaStep-5)))
	  /((exp(ceil((log(-omegaT)+5)/omegaStep)*omegaStep-5))
	    -exp(floor((log(-omegaT)+5)/omegaStep)*omegaStep-5));
      else omegaFrac = 0.;
    }

  // interpolate the logs linearly for large omegas 
  // since there the spectrum is dropping exp in q
  if (omegaT > 20.) 
    {
      qFrac = (log(q)-(floor((log(q)+5.)/qStep)*qStep-5.))/qStep;
    }
  // linear interpolation in q for small omegas
  else 
    {
      if (exp(ceil((log(q)+5.)/qStep)*qStep-5.)
	  !=exp(floor((log(q)+5.)/qStep)*qStep-5.))
	qFrac = (q-(exp(floor((log(q)+5.)/qStep)*qStep-5.)))
	  /((exp(ceil((log(q)+5.)/qStep)*qStep-5.))
	    -exp(floor((log(q)+5.)/qStep)*qStep-5.));
      else qFrac = 0.;
    }

  if (which_kind == 5 || which_kind == 7)
    {
      if (position >= 0 && iAlphas<Nalphas && iOmega<Nomega && iQ < Nq )
	rate = dGamma_qq_q->at(position);
      else
	rate = 0.;

      if (iAlphas+1 < Nalphas)
	rateAlphaUp = dGamma_qq_q->at(positionAlphaUp);
      else
	rateAlphaUp = rate;

      if (iOmega+1 < Nomega)
	rateOmegaUp = dGamma_qq_q->at(positionOmegaUp);
      else
	rateOmegaUp = rate;

      if (iQ+1 < Nq)
	rateQUp = dGamma_qq_q->at(positionQUp);
      else
	rateQUp = rate;

      if (iAlphas < Nalphas && iOmega < Nomega)
	rateAlphaUpOmegaUp = dGamma_qq_q->at(positionAlphaUpOmegaUp);
      else
	rateAlphaUpOmegaUp = rate;

      if (iAlphas < Nalphas && iQ < Nq)
	rateAlphaUpQUp = dGamma_qq_q->at(positionAlphaUpQUp);
      else
	rateAlphaUpQUp = rate;

      if (iOmega+1 < Nomega && iQ+1 < Nq)
	rateOmegaUpQUp = dGamma_qq_q->at(positionOmegaUpQUp);
      else
	rateOmegaUpQUp = rate;

      if (iAlphas < Nalphas && iOmega < Nomega && iQ < Nq)
	rateAlphaUpOmegaUpQUp = dGamma_qq_q->at(positionAlphaUpOmegaUpQUp);
      else
	rateAlphaUpOmegaUpQUp = rate;

      // used for extrapolation when the data points are too far apart
      if (omegaT > 20.)
	{ 
	  if (iQ+2 < Nq )
	    rate2QUp = dGamma_qq_q->at(position2QUp);
	  else
	    rate2QUp = rateQUp;

	  if (iAlphas < Nalphas && iQ+2 < Nq )
	    rateAlphaUp2QUp = dGamma_qq_q->at(positionAlphaUpQUp+1);
	  else
	    rateAlphaUp2QUp = rateAlphaUpQUp;

	  if (iOmega < Nomega && iQ+2 < Nq )
	    rateOmegaUp2QUp = dGamma_qq_q->at(positionOmegaUpQUp+1);
	  else
	    rateOmegaUp2QUp = rateOmegaUpQUp;

	  if (iAlphas < Nalphas && iOmega < Nomega && iQ+2 < Nq )
	    rateAlphaUpOmegaUp2QUp = dGamma_qq_q->at(positionAlphaUpOmegaUpQUp+1);
	  else
	    rateAlphaUpOmegaUp2QUp = rateAlphaUpOmegaUpQUp;
	}
    }
  else
    {
      if (position > 0 && iAlphas < Nalphas && iOmega < Nomega && iQ < Nq )
	rate = dGamma_qg_q->at(position);
      else
	rate = 0.;

      if (iAlphas+1 < Nalphas)
	rateAlphaUp = dGamma_qg_q->at(positionAlphaUp);
      else
	rateAlphaUp = rate;

      if (iOmega+1 < Nomega)
	rateOmegaUp = dGamma_qg_q->at(positionOmegaUp);
      else
	rateOmegaUp = rate;

      if (iQ+1 < Nq)
	rateQUp = dGamma_qg_q->at(positionQUp);
      else
	rateQUp = rate;

      if (iAlphas < Nalphas && iOmega < Nomega)
	rateAlphaUpOmegaUp = dGamma_qg_q->at(positionAlphaUpOmegaUp);
      else
	rateAlphaUpOmegaUp = rate;

      if (iAlphas < Nalphas && iQ < Nq)
	rateAlphaUpQUp = dGamma_qg_q->at(positionAlphaUpQUp);
      else
	rateAlphaUpQUp = rate;

      if (iOmega+1 < Nomega && iQ+1 < Nq)
	rateOmegaUpQUp = dGamma_qg_q->at(positionOmegaUpQUp);
      else
	rateOmegaUpQUp = rate;

      if (iAlphas < Nalphas && iOmega < Nomega && iQ < Nq)
	rateAlphaUpOmegaUpQUp = dGamma_qg_q->at(positionAlphaUpOmegaUpQUp);
      else
	rateAlphaUpOmegaUpQUp = rate;

      // used for extrapolation when the data points are too far apart
      if (omegaT > 20.)
	{ 
	  if (iQ+2 < Nq )
	    rate2QUp = dGamma_qg_q->at(position2QUp);
	  else
	    rate2QUp = rateQUp;

	  if (iAlphas < Nalphas && iQ+2 < Nq )
	    rateAlphaUp2QUp = dGamma_qg_q->at(positionAlphaUpQUp+1);
	  else
	    rateAlphaUp2QUp = rateAlphaUpQUp;

	  if (iOmega < Nomega && iQ+2 < Nq )
	    rateOmegaUp2QUp = dGamma_qg_q->at(positionOmegaUpQUp+1);
	  else
	    rateOmegaUp2QUp = rateOmegaUpQUp;

	  if (iAlphas < Nalphas && iOmega < Nomega && iQ+2 < Nq )
	    rateAlphaUpOmegaUp2QUp = dGamma_qg_q->at(positionAlphaUpOmegaUpQUp+1);
	  else
	    rateAlphaUpOmegaUp2QUp = rateAlphaUpOmegaUpQUp;
	}
    }
  
  if (omegaT > 0. && omegaT <= 20.)
    {
      rateOmegaAv = (1.-omegaFrac)*rate + omegaFrac*rateOmegaUp;
      rateAlphaUpOmegaAv = (1.-omegaFrac)*rateAlphaUp + omegaFrac*rateAlphaUpOmegaUp;
      rateQUpOmegaAv = (1.-omegaFrac)*rateQUp + omegaFrac*rateOmegaUpQUp;
      rateAlphaUpQUpOmegaAv = (1.-omegaFrac)*rateAlphaUpQUp + omegaFrac*rateAlphaUpOmegaUpQUp;
    }
  else if (omegaT > 20.)
    {
      if (rate != 0. && rateOmegaUp != 0.)
	rateOmegaAv = exp((1.-omegaFrac)*log(rate)+omegaFrac*log(rateOmegaUp));
      else if (rate == 0.)
	rateOmegaAv = rateOmegaUp;
      else if (rateOmegaUp == 0.)
	rateOmegaAv = rate;
      else 
	rateOmegaAv = 0.;

      if (rateAlphaUpOmegaUp != 0. && rateAlphaUp != 0.)
	rateAlphaUpOmegaAv = exp((1.-omegaFrac)*log(rateAlphaUp)
				 +omegaFrac*log(rateAlphaUpOmegaUp));
      else if (rateAlphaUp == 0.)
	rateAlphaUpOmegaAv = rateAlphaUpOmegaUp;
      else if (rateAlphaUpOmegaUp == 0.)
	rateAlphaUpOmegaAv = rateAlphaUp;
      else 
	rateAlphaUpOmegaAv = 0.;

      if (rateOmegaUpQUp != 0. && rateQUp != 0.)
	rateQUpOmegaAv = exp((1.-omegaFrac)*log(rateQUp)
			     +omegaFrac*log(rateOmegaUpQUp));
      else if (rateOmegaUpQUp == 0.)
	rateQUpOmegaAv = rateQUp;
      else if (rateQUp == 0.)
	rateQUpOmegaAv = rateOmegaUpQUp;
      else 
	rateQUpOmegaAv = 0.;

      if (rateAlphaUpOmegaUpQUp != 0. && rateAlphaUpQUp != 0.)
	rateAlphaUpQUpOmegaAv = exp((1.-omegaFrac)*log(rateAlphaUpQUp)
				    +omegaFrac*log(rateAlphaUpOmegaUpQUp));
      else if (rateAlphaUpQUp == 0.)
	rateAlphaUpQUpOmegaAv = rateAlphaUpOmegaUpQUp;
      else if (rateAlphaUpOmegaUpQUp == 0.)
	rateAlphaUpQUpOmegaAv = rateAlphaUpQUp;
      else 
	rateAlphaUpQUpOmegaAv = 0.;

      rate2QUpOmegaAv = exp((1.-omegaFrac)*log(rate2QUp)+omegaFrac*log(rateOmegaUp2QUp));
      rateAlphaUp2QUpOmegaAv = exp((1.-omegaFrac)*log(rateAlphaUp2QUp)
				   +omegaFrac*log(rateAlphaUpOmegaUp2QUp));
    }
  else if (omegaT < 0.)
    {
      rateOmegaAv = (omegaFrac)*rate+(1.-omegaFrac)*rateOmegaUp;
      rateQUpOmegaAv = (omegaFrac)*rateQUp+(1.-omegaFrac)*rateOmegaUpQUp;
      rateAlphaUpOmegaAv = (omegaFrac)*rateAlphaUp+(1.-omegaFrac)*rateAlphaUpOmegaUp;
      rateAlphaUpQUpOmegaAv = (omegaFrac)*rateAlphaUpQUp+(1.-omegaFrac)*rateAlphaUpOmegaUpQUp;
    }      

  // interpolate logs for large omega
  if (omegaT > 20.)
    {
      if (rateOmegaAv > 0.)
	{
	  rateQAv = exp((1.-qFrac)*log(rateOmegaAv)+qFrac*log(rateQUpOmegaAv));
	}
      else // use extrapolation
	{
	  slope = (log(rate2QUpOmegaAv)-log(rateQUpOmegaAv))/qStep;
	  rateQAv = exp(log(rateQUpOmegaAv)-slope*((1.-qFrac)*qStep));
	}
      if (rateAlphaUpOmegaAv > 0.)
	{
	  rateAlphaUpQAv = exp((1.-qFrac)*log(rateAlphaUpOmegaAv) + qFrac*log(rateAlphaUpQUpOmegaAv));
	}
      else  // use extrapolation
	{
	  slopeAlphaUp = (log(rateAlphaUp2QUpOmegaAv)-log(rateAlphaUpQUpOmegaAv))/qStep;
	  rateAlphaUpQAv = exp(log(rateAlphaUpQUpOmegaAv)-slopeAlphaUp*((1.-qFrac)*qStep));
	}
    }
  // interpolate linearly for small omega
  else
    {
      rateQAv = (1.-qFrac)*rateOmegaAv + qFrac*rateQUpOmegaAv;
      rateAlphaUpQAv = (1.-qFrac)*rateAlphaUpOmegaAv + qFrac*rateAlphaUpQUpOmegaAv;
    }

  result = (1.-alphaFrac)*rateQAv + alphaFrac*rateAlphaUpQAv;
  
  // the absolute normalization doesn't matter when sampling the shape
  // it only matters in "totalRate" etc.
  return result;
}


double LambertWM(double z)
{
  double w_new, w_old, ratio, e_old, tol;
  int n;

  tol = 1.0e-14;

  if(z <= -exp(-1.0))
    {
      fprintf(stderr, "LambertWM is not defined for z = %e\n", z);
      fprintf(stderr, "z needs to be bigger than %e\n", -exp(-1.0));
      exit(0);
    }

  if(z > 3.0)
    {
      w_old = log(z) - log(log(z));
    }
  else {w_old = 1.0;}
 
  w_new = 0.0;
  ratio = 1.0;
  n = 0;
  while(std::abs(ratio) > tol) 
    {
      e_old = exp(w_old);
      ratio = w_old*e_old - z;
      ratio /= ( e_old*(w_old + 1.0) - (w_old+2.0)*(w_old*e_old-z)/(2.0*w_old + 2.0) );
      w_new = w_old - ratio;
      w_old = w_new;
      n++;
      if(n > 99) 
	{
	  fprintf(stderr, "LambertWM is not converging after 100 iterations.\n");
	  fprintf(stderr, "LambertWM: z = %e\n", z);
	  fprintf(stderr, "LambertWM: w_old = %e\n", w_old);
	  fprintf(stderr, "LambertWM: w_new = %e\n", w_new);
	  fprintf(stderr, "LambertWM: ratio = %e\n", ratio);
	  exit(0);
	}
    }

  return w_new;
}// LambertW by Sangyong Jeon

