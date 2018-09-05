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

#include "Martini_inel.h"
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


Martini_inel::Martini_inel()
{
  SetId("Martini_inel");
  VERBOSE(8);

  //vectors for elastic rates:
  dGamma_qq = new vector<double>;
  dGamma_qg = new vector<double>;
  dGamma_qq_q = new vector<double>;
  dGamma_qg_q = new vector<double>;
}

Martini_inel::~Martini_inel()
{
  VERBOSE(8);
}

void Martini_inel::Init()
{
  INFO<<"Intialize Martini_inel ...";

  // Redundant (get this from Base) quick fix here for now
  tinyxml2::XMLElement *eloss= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" );
  if ( !eloss )     throw std::runtime_error("Eloss not properly initialized in XML file ...");

  tinyxml2::XMLElement *martini=eloss->FirstChildElement("Martini_inel");
  // check that all is there
  if ( !martini )     throw std::runtime_error("Martini_inel not properly initialized in XML file ...");
  if ( !martini->FirstChildElement( "name" ) )     throw std::runtime_error("Martini_inel not properly initialized in XML file ...");
  if ( !martini->FirstChildElement( "Q0" ) )     throw std::runtime_error("Martini_inel not properly initialized in XML file ...");
  if ( !martini->FirstChildElement( "alpha_s" ) )     throw std::runtime_error("Martini_inel not properly initialized in XML file ...");
  if ( !martini->FirstChildElement( "pcut" ) )     throw std::runtime_error("Martini_inel not properly initialized in XML file ...");
  if ( !martini->FirstChildElement( "hydro_Tc" ) )     throw std::runtime_error("Martini_inel not properly initialized in XML file ...");

  string s = martini->FirstChildElement( "name" )->GetText();
  JSDEBUG << s << " to be initilizied ...";

  Q0 = 1.0;
  martini->FirstChildElement("Q0")->QueryDoubleText(&Q0);

  alpha_s = 0.3;
  martini->FirstChildElement("alpha_s")->QueryDoubleText(&alpha_s);
    
  pcut = 2.0;
  martini->FirstChildElement("pcut")->QueryDoubleText(&pcut);

  hydro_Tc = 0.16;
  martini->FirstChildElement("hydro_Tc")->QueryDoubleText(&hydro_Tc);

  g = sqrt(4.*M_PI*alpha_s);
  alpha_em = 1./137.;

  // Path to additional data
  if ( !martini->FirstChildElement( "path" ) )     throw std::runtime_error("Martini_inel not properly initialized in XML file ...");
  PathToTables=martini->FirstChildElement( "path" )->GetText();

  // Initialize random number distribution
  ZeroOneDistribution = uniform_real_distribution<double> { 0.0, 1.0 };
  
  readRadiativeRate(&dat, &Gam);
}

void Martini_inel::DoEnergyLoss(double deltaT, double Time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
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
  double omega, q;           // transferred energy/momentum for scattering
  double xx, yy, zz;         // position of initial parton (pIn)
  FourVector pVec, pVecNew;  // 4 vectors for momenta before & after process
  FourVector pVecRest;       // 4 vector in the rest frame of fluid cell
  FourVector pVecNewRest;
  FourVector kVec;           // 4 vector for momentum of radiated particle
  FourVector xVec;           // 4 vector for position (for next time step!)
  double velocity_jet[4];    // jet velocity for MATTER
  double eta;                // pseudo-rapidity
  
  // flow info
  double vx, vy, vz;         // 3 components of flow velocity
  double T;                  // Temperature of fluid cell
  double beta, gamma;        // flow velocity & gamma factor
  double cosPhi;             // angle between flow and particle
  double cosPhiRest;         // angle between flow and particle in rest frame
  double boostBack;          // factor for boosting back to lab frame
  
  for (int i=0;i<pIn.size();i++) {

    // Particle infomration
    Id = pIn[i].pid();

    px = pIn[i].px();
    py = pIn[i].py();
    pz = pIn[i].pz();

    // In MARTINI, particles are all massless and on-shell
    pAbs = sqrt(px*px+py*py+pz*pz);
    pVec = FourVector ( px, py, pz, pAbs );

    xx = pIn[i].x_in().x();
    yy = pIn[i].x_in().y();
    zz = pIn[i].x_in().z();

    eta = pIn[i].eta();

    // Extract fluid properties
    std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;
    GetHydroCellSignal(Time, xx, yy, zz, check_fluid_info_ptr);
    VERBOSE(8)<< MAGENTA<<"Temperature from Brick (Signal) = "
	      <<check_fluid_info_ptr->temperature;

    vx = check_fluid_info_ptr->vx;
    vy = check_fluid_info_ptr->vy;
    vz = check_fluid_info_ptr->vz;
    T = check_fluid_info_ptr->temperature;

    beta = sqrt( vx*vx + vy*vy + vz*vz );

    // Only accept low t particles
    if (pIn[i].t() > Q0*Q0 + rounding_error || T < hydro_Tc) continue;
    TakeResponsibilityFor ( pIn[i] ); // Generate error if another module already has responsibility.

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

    velocity_jet[0]=1.0;
    velocity_jet[1]=pIn[i].jet_v().x();
    velocity_jet[2]=pIn[i].jet_v().y();
    velocity_jet[3]=pIn[i].jet_v().z();

    int process = DetermineProcess(pRest, T, deltaT, Id);
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
	pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity

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
		pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity

	      }

	    if (kRest > pcut)
	      {
		k = kRest*boostBack;
		kVec.Set( (px/pAbs)*k, (py/pAbs)*k, (pz/pAbs)*k, k );
		pOut.push_back(Parton(0, Id, 0, kVec, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
                pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity
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
                pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity
	      }

	    // photon doesn't have energy threshold; No absorption into medium
	    // However, we only keep positive energy photons
	    if (kRest > 0.)
	      {
		k = kRest*boostBack;
		kVec.Set( (px/pAbs)*k, (py/pAbs)*k, (pz/pAbs)*k, k );
		pOut.push_back(Parton(0, Id, 0, kVec, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
                pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity
	      }

	    return;
	  }
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
                pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity
	      }

	    if (kRest > pcut)
	      {
		k = kRest*boostBack;
		kVec.Set( (px/pAbs)*k, (py/pAbs)*k, (pz/pAbs)*k, k );
		pOut.push_back(Parton(0, Id, 0, kVec, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
                pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity
	      }

	    return;
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
                pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity
	      }

	    if (kRest > pcut)
	      {
		k = kRest*boostBack;
		kVec.Set( (px/pAbs)*k, (py/pAbs)*k, (pz/pAbs)*k, k );
		pOut.push_back(Parton(0, Id, 0, kVec, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
                pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity
	      }

	    return;
	  }
      } // Id==21
  } // particle loop
}

int Martini_inel::DetermineProcess(double pRest, double T, double deltaT, int Id)
{

  double dT = deltaT/hbarc;   // time step in [GeV^(-1)]

  // get the rates for each process
  // total Probability = dT*Rate
  RateRadiative rateRad;
  rateRad = getRateRadTotal(pRest, T);

  // evolution for quark (u, d, s)
  if (std::abs(Id) == 1 || std::abs(Id) == 2 || std::abs(Id) == 3)
    {
      // multiplying by (ef/e)^2
      if (abs(Id) == 1)
	{
	  rateRad.qqgamma *= 4./9.;
	}
      else
	{
	  rateRad.qqgamma *= 1./9.;
	}

      double totalQuarkProb = 0.;

      if (pRest > pcut) totalQuarkProb += rateRad.qqg*dT;

      // warn if total probability exceeds 1
      if (totalQuarkProb > 1.){
	WARN << " : Total Probability for quark processes exceeds 1 ("
	     << totalQuarkProb << "). "
	     << " : Most likely this means you should choose a smaller deltaT in the xml (e.g. 0.01).";
	//throw std::runtime_error ("Martini_inel probability problem.");
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

      // warn if total probability exceeds 1
      if (totalGluonProb > 1.){
	WARN << " : Total Probability for gluon processes exceeds 1 ("
	     << totalGluonProb << "). "
	     << " : Most likely this means you should choose a smaller deltaT in the xml (e.g. 0.01).";
	//throw std::runtime_error ("Martini_inel probability problem.");
      }

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

RateRadiative Martini_inel::getRateRadTotal(double pRest, double T)
{
  RateRadiative rate;

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

  return rate;
}

RateRadiative Martini_inel::getRateRadPos(double u, double T)
{
  RateRadiative rate;

  rate.qqg = (0.5322 - 3.1037/(u*u) + 2.0139/u - 0.9417/sqrt(u))*pow(g, 4.)*T;
  rate.ggg = (1.1923 - 11.5250/(u*u*u) + 3.3010/u - 1.9049/sqrt(u))*pow(g, 4.)*T;
  rate.gqq = (0.0004656 - 0.04621/(u*u) + 0.0999/u - 0.08171/pow(u,0.8) 
              + 0.008090/pow(u,0.2) - 1.2525*pow(10.,-8.)*u)*pow(g, 4.)*T*nf; 
  rate.qqgamma = 0.;

  return rate;
}

RateRadiative Martini_inel::getRateRadNeg(double u, double T)
{
  RateRadiative rate;

  rate.qqg = (0.3292 - 0.6759/(u*u) + 0.4871/pow(u,1.5) - 0.05393/u + 0.007878/sqrt(u))
    *pow(g, 4.)*T;
  rate.ggg = (0.7409 + 1.8608/(u*u*u) - 0.1353/(u*u) + 0.1401/u)*pow(g, 4.)*T;
  rate.gqq = (0.03215/(u*u*u) + 0.01419/(u*u) + 0.004338/u - 0.00001246/sqrt(u))
    *pow(g, 4.)*T*nf;
  rate.qqgamma = 0.;

  return rate;
}

double Martini_inel::getNewMomentumRad(double pRest, double T, int process)
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
	      y = 2.5/(LambertW3(2.59235*pow(10.,23.)*exp(-100.*randA)));

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
	      y = 5./(LambertW3(2.68812*pow(10., 45.)*exp(-50.*randA)));

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
double Martini_inel::area(double y, double u, int posNegSwitch, int process)
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

double Martini_inel::function(double u, double y, int process)
{
  if (process == 1)      return getRate_qqg(u, y);
  else if (process == 2) return getRate_qqgamma(u, y);
  else if (process == 3) return getRate_ggg(u, y);
  else if (process == 4) return getRate_gqq(u, y);

  return 0.;
}


// Reads in the binary stored file of dGamma values
void Martini_inel::readRadiativeRate(Gamma_info *dat, dGammas *Gam)
{
  FILE *rfile;
  string filename;
  filename = PathToTables+"radgamma";

  INFO << "Reading rates of inelastic collisions from file ";
  INFO << filename.c_str() << " ... ";
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



double Martini_inel::getRate_qqg(double p, double k)
{
  return use_table(p, k, Gam.qqg, 0);
}

double Martini_inel::getRate_gqq(double p, double k)
{
  if (k < p/2.) return use_table(p, k, Gam.gqq, 1);
  else return 0.;
}

double Martini_inel::getRate_ggg(double p, double k)
{
  if (k < p/2.) return use_table(p, k, Gam.ggg, 2);
  else return 0.;
}

double Martini_inel::getRate_qqgamma(double p, double k)
{
  return use_table(p, k, Gam.qqgamma, 3);
}

double Martini_inel::use_table(double p, double k, double dGamma[NP][NK], int which_kind)
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


double LambertW3(double z)
{
  double w_new, w_old, ratio, e_old, tol;
  int n;

  tol = 1.0e-14;

  if(z <= -exp(-1.0))
    {
      WARN << "LambertW is not defined for z = " << z;
      WARN << "z needs to be bigger than " << -exp(-1.0);
      throw std::runtime_error("LambertW small z problem");
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
          WARN << "LambertW is not converging after 100 iterations.";
          WARN << "LambertW: z = " << z;
          WARN << "LambertW: w_old = " << w_old;
          WARN << "LambertW: w_new = " << w_new;
          WARN << "LambertW: ratio = " << ratio;
          throw std::runtime_error("LambertW not conversing");
	}
    }

  return w_new;
}// LambertW by Sangyong Jeon
