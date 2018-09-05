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

#include "Tequila.h"
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


Tequila::Tequila()
{
  SetId("Tequila");
  VERBOSE(8);

}

Tequila::~Tequila()
{
  VERBOSE(8);
}

void Tequila::Init()
{
  INFO<<"Intialize Tequila ...";

  // Redundant (get this from Base) quick fix here for now
  tinyxml2::XMLElement *eloss= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" );
  if ( !eloss )     throw std::runtime_error("Eloss not properly initialized in XML file ...");

  tinyxml2::XMLElement *martini=eloss->FirstChildElement("Tequila");
  // check that all is there
  if ( !martini )     throw std::runtime_error("Tequila not properly initialized in XML file ...");
  if ( !martini->FirstChildElement( "name" ) )     throw std::runtime_error("Tequila not properly initialized in XML file ...");
  if ( !martini->FirstChildElement( "Q0" ) )     throw std::runtime_error("Tequila not properly initialized in XML file ...");
  if ( !martini->FirstChildElement( "alpha_s" ) )     throw std::runtime_error("Tequila not properly initialized in XML file ...");
  if ( !martini->FirstChildElement( "pcut" ) )     throw std::runtime_error("Tequila not properly initialized in XML file ...");
  if ( !martini->FirstChildElement( "hydro_Tc" ) )     throw std::runtime_error("Tequila not properly initialized in XML file ...");
  if ( !martini->FirstChildElement( "path" ) )     throw std::runtime_error("Tequila not properly initialized in XML file ...");
  if ( !martini->FirstChildElement( "omega_over_T_cutoff" ) )     throw std::runtime_error("Tequila not properly initialized in XML file ...");

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

  omega_over_T_cutoff = g;
  martini->FirstChildElement("omega_over_T_cutoff")->QueryDoubleText(&omega_over_T_cutoff);

  // Path to additional data
  path_to_tables=martini->FirstChildElement( "path" )->GetText();

  // Initialize random number distribution
  ZeroOneDistribution = uniform_real_distribution<double> { 0.0, 1.0 };

  allocate_memory_for_radiative_rate_table();

  // Load differential rate
  std::cout << "Loading differential collinear rate...\n";
  const double alpha_EM=1./137.; //Not currently a parameter
  const int Nf=3;
  //const std::string location_of_pretabulated_collinear_rates="./Tequila/";
  // Either compute or read from file the collinear differential rates and load them into member array "differential_rate_p_omega_qperp[][][]"
  load_differential_rate(alpha_s, alpha_EM, Nf, path_to_tables);

  // Compute total rate from differential rate and save in member array "rate_p[]"
  std::cout << "Computing integrated collinear rate...\n";
  evaluate_integrated_rate(omega_over_T_cutoff, differential_rate_qqg_p_omega_qperp,rate_qqg_p);
  evaluate_integrated_rate(omega_over_T_cutoff, differential_rate_ggg_p_omega_qperp,rate_ggg_p);
  evaluate_integrated_rate(omega_over_T_cutoff, differential_rate_gqq_p_omega_qperp,rate_gqq_p);
  std::cout << "Done computing integrated collinear rate.\n";

  
//  readRadiativeRate(&dat, &Gam);
}

void Tequila::DoEnergyLoss(double deltaT, double Time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
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

    enum radiative_process_type process = DetermineProcess(pRest, T, deltaT, Id);
    VERBOSE(8)<< MAGENTA
	      << "Time = " << Time << " Id = " << Id
	      << " process = " << process << " T = " << T
	      << " pAbs = " << pAbs << " " << px << " " << py << " " << pz 
	      << " | position = " << xx << " " << yy << " " << zz;


    // Do nothing for this parton at this timestep
    if (process == none) 
      {
	pOut.push_back(Parton(0, Id, 0, pVec, xVec));
	pOut[pOut.size()-1].set_form_time(0.);
	pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity

	return;
      }
    if (std::abs(Id) == 1 || std::abs(Id) == 2 || std::abs(Id) == 3)
      {
	// quark radiating gluon (q->qg)
	if (process == qqg)
	  {
	    if (pRest/T < AMYpCut) return;

	    // sample radiated parton's momentum
	    double omega, qperp;
	    sample_dgamma_dwdq(pRest, T,differential_rate_qqg_p_omega_qperp, omega, qperp);
	    kRest = omega; 
	    //kRest = getNewMomentumRad(pRest, T, process);
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
      }
    else if (Id == 21)
      {
	// gluon radiating gluon (g->gg)
	if (process == ggg)
	  {
	    if (pRest/T < AMYpCut) return;

	    // sample radiated parton's momentum
	    double omega, qperp;
	    sample_dgamma_dwdq(pRest, T,differential_rate_ggg_p_omega_qperp, omega, qperp);
	    kRest = omega; 
	    //kRest = getNewMomentumRad(pRest, T, process);
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
	if (process == gqq)
	  {
	    if (pRest/T < AMYpCut) return;

	    // sample radiated parton's momentum
	    double omega, qperp;
	    sample_dgamma_dwdq(pRest, T,differential_rate_gqq_p_omega_qperp, omega, qperp);
	    kRest = omega; 
	    //kRest = getNewMomentumRad(pRest, T, process);
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

enum radiative_process_type Tequila::DetermineProcess(double pRest, double T, double deltaT, int Id)
{

 double dT = deltaT/hbarc;   // time step in [GeV^(-1)]

// // get the rates for each process
// // total Probability = dT*Rate
// RateRadiative rateRad;
// rateRad = getRateRadTotal(pRest, T);

  // evolution for quark (u, d, s)
  if (std::abs(Id) == 1 || std::abs(Id) == 2 || std::abs(Id) == 3)
    {

      double totalQuarkProb = 0.;

      const double rate_qqg=rate(pRest, T, rate_qqg_p);

//      if (pRest > pcut) totalQuarkProb += rateRad.qqg*dT;
      if (pRest > pcut) totalQuarkProb += rate_qqg*dT;

      // warn if total probability exceeds 1
      if (totalQuarkProb > 1.){
	WARN << " : Total Probability for quark processes exceeds 1 ("
	     << totalQuarkProb << "). "
	     << " : Most likely this means you should choose a smaller deltaT in the xml (e.g. 0.01).";
	//throw std::runtime_error ("Tequila probability problem.");
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
	      //Prob = rateRad.qqg*dT/totalQuarkProb;
	      Prob = rate_qqg*dT/totalQuarkProb;
	      if (accumProb <= randProb && randProb < (accumProb + Prob))
		return qqg;

	      //accumProb += Prob;
	      //Prob = rateRad.qqgamma*dT/totalQuarkProb;
	      //if (accumProb <= randProb && randProb < (accumProb + Prob))
	      //  return 2;
	    }

	}
      else
	{
	  // nothing happens to quark
	  return none;
	}
    }
  // evolution for gluon
  else if (Id == 21)
    {
      double totalGluonProb = 0.;

      const double rate_gqq=rate(pRest, T, rate_gqq_p);
      const double rate_ggg=rate(pRest, T, rate_ggg_p);

      //if (pRest > pcut) totalGluonProb += (rateRad.ggg + rateRad.gqq)*dT;
      if (pRest > pcut) totalGluonProb += (rate_gqq+rate_ggg)*dT;

      // warn if total probability exceeds 1
      if (totalGluonProb > 1.){
	WARN << " : Total Probability for gluon processes exceeds 1 ("
	     << totalGluonProb << "). "
	     << " : Most likely this means you should choose a smaller deltaT in the xml (e.g. 0.01).";
	//throw std::runtime_error ("Tequila probability problem.");
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
	      Prob = rate_ggg*dT/totalGluonProb;
	      if (accumProb <= randProb && randProb < (accumProb + Prob))
		return ggg;

	      accumProb += Prob;
	      Prob = rate_gqq*dT/totalGluonProb;
	      if (accumProb <= randProb && randProb < (accumProb + Prob))
		return gqq;
	    }

	}
      else
	{
	  // nothing happens to gluon
	  return none;
	}
    }

  // if parton is other than u,d,s,g, do nothing
  return none;
}

void Tequila::allocate_memory_for_radiative_rate_table() {
	// Allocate memory for differential and integrated rate
	rate_ggg_p = new double [nb_points_in_p];
	rate_gqq_p = new double [nb_points_in_p];
	rate_qqg_p = new double [nb_points_in_p];
	//maximum_differential_rate = new double [nb_points_in_p];
	differential_rate_ggg_p_omega_qperp = new double ** [nb_points_in_p];
	differential_rate_gqq_p_omega_qperp = new double ** [nb_points_in_p];
	differential_rate_qqg_p_omega_qperp = new double ** [nb_points_in_p];
	for(int ip=0;ip<nb_points_in_p; ip++) {
		differential_rate_ggg_p_omega_qperp[ip]=new double * [nb_points_in_omega];
		differential_rate_gqq_p_omega_qperp[ip]=new double * [nb_points_in_omega];
		differential_rate_qqg_p_omega_qperp[ip]=new double * [nb_points_in_omega];
		for(int iomega=0;iomega<nb_points_in_omega; iomega++) {
			 differential_rate_ggg_p_omega_qperp[ip][iomega]=new double [nb_points_in_qperp];
			 differential_rate_gqq_p_omega_qperp[ip][iomega]=new double [nb_points_in_qperp];
			 differential_rate_qqg_p_omega_qperp[ip][iomega]=new double [nb_points_in_qperp];
		}
	}

}

//void Tequila::load_radiative_rate_tables() {
//	//const struct ERateParam &param, const ERateCase &kcs) {
////
////	//Isn't this check redundant? This function is clearly just for collinear rates...
////	if (kcs != kGColl) {
////		printf("Bad rate selection in switch in %s:%d, received %d\n", __FILE__, __LINE__, (int)kcs) ;
////	}
////	fCase = kcs ;
////
////	//Gluon in, gluon out
////	fId0 = 21 ; 
////	fIdOut0 = 21 ;
////	//   fM2Case = PSTabulator::kGToGPureglue ;
////
////	// Allocate memory for differential and integrated rate
////	rate_p = new double [nb_points_in_p];
////	//maximum_differential_rate = new double [nb_points_in_p];
////	differential_rate_p_omega_qperp = new double ** [nb_points_in_p];
////	for(int ip=0;ip<nb_points_in_p; ip++) {
////		differential_rate_p_omega_qperp[ip]=new double * [nb_points_in_omega];
////		for(int iomega=0;iomega<nb_points_in_omega; iomega++) {
////			 differential_rate_p_omega_qperp[ip][iomega]=new double [nb_points_in_qperp];
////		}
////	}
////
//	// Load differential rate
//	std::cout << "Loading differential collinear rate...\n";
//	const double alpha_EM=1./137.; //Not currently a parameter
//	//Either compute or read from file the collinear differential rates and load them into member array "differential_rate_p_omega_qperp[][][]"
//	load_differential_rate(param.alphas0(), alpha_EM, param.Nf(), param.fLocation_of_pretabulated_collinear_rates);
////	for(int i=0;i<nb_points_in_qperp;i++) std::cout << "dG/dp domega dq(q=" << get_qperp_over_T_from_index(i) << ")=" << differential_rate_p_omega_qperp[0][0][i] << "\n";
//
//	//Compute total rate from differential rate and save in member array "rate_p[]"
//	std::cout << "Computing integrated collinear rate...\n";
//	evaluate_integrated_rate(param.fCollinear_omega_over_T_cut);
//	std::cout << "Done computing integrated collinear rate.\n";
//	//for(int i=0;i<nb_points_in_p;i++) std::cout << "dG/dp(p=" << get_p_over_T_from_index(i) << ")=" << rate_p[i] << "\n"; // " vs " << omega_over_T_max(get_p_over_T_from_index(i))-fOmegaMinCutoff << "\n";
//	//for(int i=0;i<nb_points_in_p;i++) std::cout << "dG/dp(p=" << get_p_over_T_from_index(i) << ")=" << rate_p[i] << "\n";
//
//}

//Load rate differential in incoming parton energy p and radiated parton energy omega (not differential in transverse momentum q)
//The table contains dGamma/(dx domega)= gs^4*use_table(), where use_table() is the function from Moore's code that tabulates the collinear rate
void Tequila::load_differential_rate(const double alpha_s, const double alpha_EM, const int Nf, const std::string location_of_collinear_rates) {

	// Check if the tabulated rate is already available for these parameters
	// The format is "dGamma_dp_domega_NfX_alphasYYY" where "X" is the number of flavours and "YYY" is the value of alpha_s
	std::stringstream filename;
	filename << location_of_collinear_rates << "dGamma_dp_domega_Nf" << Nf << "_alphaEM" << alpha_EM << "_alphaS" << alpha_s;
	//Open file
	std::ifstream rate_file;
	rate_file.open(filename.str().c_str(),std::fstream::in);

	if (!rate_file.is_open()) {
		std::cout << "Can't open rate table file, aborting...\n";
		exit(1);
	}

	//Struct [from Moore's code] to store rates
	dGammas Moore_rate_arrays;
	Gamma_info Moore_rate_info;
	// Save file with proper name
	sprintf(Moore_rate_info.in_fname, "%s", filename.str().c_str());

	// If not, run Guy's program given parameters
	if (!rate_file.is_open()) {

		//
		std::cout << "Pre-tabulated collinear rate not available. Tabulating it now...\n";
	
		//Pre-initialized in data: Nc, Nf, alpha_s, alpha
		Moore_rate_info.Nf=Nf;
		Moore_rate_info.Nc=3;
		Moore_rate_info.alpha_s=alpha_s;
		Moore_rate_info.alpha=alpha_EM;
		build_table(&Moore_rate_info , &Moore_rate_arrays);
		
		std::cout << "... and writing it into file \"" << filename.str().c_str() << "\" for future use.\n";
		write_table(&Moore_rate_info , &Moore_rate_arrays);

	}
	else {
		std::cout << "Collinear rate available for value of alpha_s, alpha_EM and Nf. Reading from file.\n";

		rate_file.close();
		// Read rate file in temporary array differential in p and omega using Moore's functions
		read_table(&Moore_rate_info,&Moore_rate_arrays);
		std::cout << "Collinear rates read.\n";
	}

	const double gs4=alpha_s*alpha_s*(16*M_PI*M_PI);


	// Populate member array "differential_rate_p_omega_qperp[]" with p and omega rate from Moore's code and Gaussian distribution for qperp
	// Some useful parameters for the (temporary) qperp distrib
	// The (temporary) Gaussian envelope used in the q_perp direction (since the q_perp distribution is not currently known)
	const double envelope_sigma_sqr= 4.*M_PI*alpha_s/3.*(3 + Nf*0.5); //Debye mass
	for(int ip=0;ip<nb_points_in_p; ip++) {
		const double p_over_T_val=get_p_over_T_from_index(ip);
		for(int iomega=0;iomega<nb_points_in_omega; iomega++) {
			const double omega_over_T_val=get_omega_over_T_from_index(ip,iomega);
			//Object "Moore_rate_arrays" does not actually contains the rate, but rather the rate stripped of various factors 
			//Moore's function "use_table()" returns the actual rate, after multiplications by the proper factors
			//Using "use_table()" to do this is kind-of sub-optimal, but speed shouldn't be too much of an issue here
			//Need to multiply by g_s^4 since the factor is stripped from the rate in Moore's program
			const double tmp_rate_ggg=gs4*use_table(p_over_T_val , omega_over_T_val , Moore_rate_arrays.dGamma_ggg , 2 );
			const double tmp_rate_gqq=gs4*use_table(p_over_T_val , omega_over_T_val , Moore_rate_arrays.dGamma_gqq , 1 );
			const double tmp_rate_qqg=gs4*use_table(p_over_T_val , omega_over_T_val , Moore_rate_arrays.dGamma, 0 );

			for(int iq=0; iq<nb_points_in_qperp;iq++) {	
				// q_perp envelope function: Gaussian of width sigma \propto m_D?
				// Integrate[ (2 Pi r) Exp[-r^2/sigma^2]/Sqrt[Pi sigma^2]^2, {r, 0, Infinity}]=1
				//jacobian between E dG/dq^3 and dG/(dqperp domega)???
				const double qperp_over_T=get_qperp_over_T_from_index(iq);
				//const double envelope=exp(-qperp_over_T*qperp_over_T/envelope_sigma_sqr)/(M_PI*envelope_sigma_sqr);
				// Hack to make the integral over q_perp unity
				//const double jacobian=2*M_PI*qperp_over_T;
				const double envelope=1./(qperp_over_T_max()); 
				const double rate_ggg=tmp_rate_ggg*envelope;
				const double rate_gqq=tmp_rate_gqq*envelope;
				const double rate_qqg=tmp_rate_qqg*envelope;
				//if (ip ==0 && iomega == 0) std::cout << "qperp_over_T=" << qperp_over_T << " & rate=" << rate << "\n";
				//assign rate
				differential_rate_ggg_p_omega_qperp[ip][iomega][iq]=rate_ggg;
				differential_rate_gqq_p_omega_qperp[ip][iomega][iq]=rate_gqq;
				differential_rate_qqg_p_omega_qperp[ip][iomega][iq]=rate_qqg;
				//differential_rate_p_omega_qperp[ip][iomega][iq]=envelope;
			}
		}
	}

}


//Find position in array corresponding to value of p_over_T and omega_over_T
//Matches Moore's code (it must!)
double Tequila::get_index_from_omega_over_T(const double p_over_T, const double omega_over_T) {

	double b;
	
	if ( omega_over_T < 2 )
	{
		if ( omega_over_T < -1 )
		{
			if ( omega_over_T < -2 )
				b = 60 + 5*omega_over_T;
			else
				b = 70+10*omega_over_T;
		}
		else
		{
			if ( omega_over_T < 1 )
				b = 80 + 20*omega_over_T;
			else
				b = 90 + 10*omega_over_T;
		}
	}
	else if ( omega_over_T < p_over_T-2 )
	{ /* This is that tricky middle ground. */
		b = 190 - 10*log ( 1.000670700260932956l / 
				( 0.0003353501304664781l + (omega_over_T-2) / (p_over_T-4) ) - 1 );
	}
	else
	{
		if ( omega_over_T < p_over_T+1 )
		{
			if ( omega_over_T < p_over_T-1 )
				b = 290 + 10*(omega_over_T-p_over_T);
			else
				b = 300 + 20*(omega_over_T-p_over_T);
		}
		else
		{
			if ( omega_over_T < p_over_T+2 )
				b = 310 + 10*(omega_over_T-p_over_T);
			else
				b = 320 + 5*(omega_over_T-p_over_T);
		}
	}

	return b;

}


//Find the value of p_over_T and omega_over_T for a given position in array
//Matches Moore's code (it must!)
double Tequila::get_omega_over_T_from_index(const int index_p_over_T, const int index_omega_over_T) {

	//Need p_over_T to get omega_over_T
	double p_over_T=get_p_over_T_from_index(index_p_over_T);
	double omega_over_T;

	if ( index_omega_over_T < 50 )        /* spaced by 0.2  from -12 to -2 */
		omega_over_T = -12 + index_omega_over_T * 0.2;
	else if ( index_omega_over_T < 60 )   /* spaced by 0.1  from -2  to -1 */
		omega_over_T = -2 + (index_omega_over_T-50) * 0.1;
	else if ( index_omega_over_T < 100 )  /* spaced by 0.05 from -1  to +1 */
		omega_over_T = -1 + (index_omega_over_T-60) * 0.05;
	else if ( index_omega_over_T < 110 )  /* spaced by 0.1  from +1  to +2 */
		omega_over_T = 1 + (index_omega_over_T-100) * 0.1;
	else if ( index_omega_over_T < 270 )  /* spaced complicated, +2 to p_over_T-2 */
	{
		omega_over_T = 0.1 * (index_omega_over_T-190);
		omega_over_T = 2 + (p_over_T-4) * ( -0.0003353501304664781l
				+ 1.000670700260932956l / (1+exp(-omega_over_T)) );
	}
	else if ( index_omega_over_T < 280 )  /* spaced by 0.1  from p_over_T-2 to p_over_T-1 */
		omega_over_T = p_over_T - 2 + 0.1 * (index_omega_over_T-270);
	else if ( index_omega_over_T < 320 )  /* spaced by 0.05 from p_over_T-1 to p_over_T+1 */
		omega_over_T = p_over_T + 0.05 * (index_omega_over_T - 300);
	else if ( index_omega_over_T < 330 )  /* spaced by 0.1  from p_over_T+1 to p_over_T+2 */
		omega_over_T = p_over_T + 0.1 * (index_omega_over_T - 310);
	else                   /* spaced by 0.2  from p_over_T+2 to p_over_T+12 */
		omega_over_T = p_over_T + 0.2 * (index_omega_over_T - 320);

	return omega_over_T;

}

//Evaluated integrated rate, with cut-off "omega_over_T_cut" on omega, from differential rate stored in member object "differential_rate_p_omega_qperp[][]"
//Also save the maximum value of the rate
void Tequila::evaluate_integrated_rate(double omega_over_T_cut, double *** differential_rate, double * integrated_rate) {

	//current omega/T integration not very good if omega_over_T_cut<0.1
	if (omega_over_T_cut<0.1) std::cout << "Warning: omega/T integration is not very good for omega/T cut-off smaller than 0.1\n";

	//loop over all values of "p"
	for(int ip=0;ip<nb_points_in_p; ip++) {

		const double p_over_T=get_p_over_T_from_index(ip);

		// integrate omega/T from -infinity to omega_over_T_cut, and omega_over_T_cut to p/2T
		// the discretization is not uniform in omega/T
		// let's just use the GSL interpolation routine to integrate, since performance is not an issue here
		double integral_omega=0.0;

		// Maximum point in omega/T necessary
		int pOver2T_cut_pos_position_int=ceil(get_index_from_omega_over_T(p_over_T,p_over_T/2.0))+2;

		//Arrays to store the omega/T rates
		double * omega_rate_array = new double [pOver2T_cut_pos_position_int];
		double * omega_position_array = new double [pOver2T_cut_pos_position_int];

		// Fill an array with the qperp-integrated values for each omega points
		for(int iomega=0;iomega<pOver2T_cut_pos_position_int; iomega++) {

			omega_position_array[iomega]=get_omega_over_T_from_index(ip,iomega);

			// integrate in q_perp
			double integral_qperp=0.0;

			// track the position of various timesteps to compute the weight of each step properly with a minimum number of call to function "get_qperp_over_T_from_index()"
			// trapeze integration not very good for radial coordinate. might require improvements.
			double prev_qperp=0.0;
			double curr_qperp=get_qperp_over_T_from_index(0);
			double next_qperp;
			for(int iq=0; iq<nb_points_in_qperp;iq++) {	
				next_qperp=get_qperp_over_T_from_index(iq+1);
				// get weight
				const double weight_qperp=non_uniform_trapeze_weights(iq,nb_points_in_qperp,prev_qperp,curr_qperp,next_qperp);
				// disable jacobian for now
				// jacobian for radial integration
				//const double jacobian=2*M_PI*curr_qperp;
				const double jacobian=1.0;

				// add contribution from this timestep
				//integral_qperp+=weight_qperp*jacobian*differential_rate_p_omega_qperp[ip][iomega][iq];
				integral_qperp+=weight_qperp*jacobian*differential_rate[ip][iomega][iq];
				//update positions for q_perp
				prev_qperp=curr_qperp;
				curr_qperp=next_qperp;
			}

			//std::cout << "ip=" << ip << ", iomega=" << iomega << ", rate=" << integral_qperp << "\n";

			omega_rate_array[iomega]=integral_qperp;

		}

		// initialise GSL interpolator
		gsl_interp * interp = gsl_interp_alloc( gsl_interp_akima, pOver2T_cut_pos_position_int);
	//	gsl_interp * interp = gsl_interp_alloc( gsl_interp_linear, pOver2T_cut_pos_position_int);
	//	gsl_interp * interp = gsl_interp_alloc(gsl_interp_cspline, pOver2T_cut_pos_position_int);
		gsl_interp_accel *acc = gsl_interp_accel_alloc ();
		gsl_interp_init(interp,  omega_position_array, omega_rate_array, pOver2T_cut_pos_position_int); 

		// integral in omega from -infinity to -omega_over_T_cut
		integral_omega+=gsl_interp_eval_integ(interp, omega_position_array, omega_rate_array, omega_over_T_min(p_over_T), -1*omega_over_T_cut, acc);
		// integral in omega from omega_over_T_cut to p_over_T/2
		integral_omega+=gsl_interp_eval_integ(interp, omega_position_array, omega_rate_array , omega_over_T_cut, p_over_T/2., acc);

		// free memory
		gsl_interp_free(interp);
		gsl_interp_accel_free(acc);
		delete [] omega_position_array;
		delete [] omega_rate_array;

		//std::cout << ip << " "  << integral_omega << "\n";

		//rate_p[ip]=integral_omega;
		integrated_rate[ip]=integral_omega;
	}



}

// weights for the trapezoid rule on a non-uniform grid
// 1/2 Sum[f[getQfromIndex[position]]Which[0==position,(getQfromIndex[position+1]-getQfromIndex[position]),size_of_array-1==position,(getQfromIndex[position]-getQfromIndex[position-1]),True,(getQfromIndex[position+1]-getQfromIndex[position-1])],{position,0,Qnum-1}]
// 1/2 Sum[f[getQfromIndex[position]]Which[0==position,(next_x-curr_x),size_of_array-1==i,(curr_x-prev_x),True,(next_x-prev_x)],{i,0,size_of_array-1}]
double Tequila::non_uniform_trapeze_weights(int position, int size_of_array, double prev_x, double curr_x, double next_x) {

	double weight=0.5;
	
	if (0 == position) {
		weight*=(next_x-curr_x);	
	}
	else if (size_of_array -1 == position) {
		weight*=(curr_x-prev_x);	
	}
	else {
		weight*=(next_x-prev_x);
	}
	return weight;

}

//Returns T^2 dGamma/(dx domega d^2 qperp)
double Tequila::differential_rate(const double p_over_T, const double omega_over_T, const double qperp_over_T, double *** differential_rate_p_omega_qperp) {


	//tri-linear interpolation
	//somehow I thought this function would be shorter...

	if ((p_over_T<p_over_T_min())||(p_over_T>p_over_T_max())) return 0.0;
	if ((omega_over_T<omega_over_T_min(p_over_T))||(omega_over_T>omega_over_T_max(p_over_T))||(qperp_over_T>qperp_over_T_max())) return 0.0;

	//first, get position in grid of where rate is located
	const double tmp_p=get_index_from_p_over_T(p_over_T);
	const double tmp_omega=get_index_from_omega_over_T(p_over_T,omega_over_T);
	const double tmp_qperp=get_index_from_qperp_over_T(qperp_over_T);

	const int pos_array_p_over_T=floor(tmp_p);
	const int pos_array_omega_over_T=floor(tmp_omega);
	const int pos_array_qperp_over=floor(tmp_qperp);

	//actual positions of the grid points around the desired value
	const double p_over_T_low_val=get_p_over_T_from_index(pos_array_p_over_T);
	const double omega_over_T_low_val=get_omega_over_T_from_index(pos_array_p_over_T,pos_array_omega_over_T);
	const double qperp_over_T_low_val=get_qperp_over_T_from_index(pos_array_qperp_over);
	const double p_over_T_high_val=get_p_over_T_from_index(pos_array_p_over_T+1);
	const double omega_over_T_high_val=get_omega_over_T_from_index(pos_array_p_over_T,pos_array_omega_over_T+1);
	const double qperp_over_T_high_val=get_qperp_over_T_from_index(pos_array_qperp_over+1);

	//value of the rate at the above gridpoints
	const double v000=differential_rate_p_omega_qperp[pos_array_p_over_T][pos_array_omega_over_T][pos_array_qperp_over];
	const double v001=differential_rate_p_omega_qperp[pos_array_p_over_T][pos_array_omega_over_T][pos_array_qperp_over+1];
	const double v010=differential_rate_p_omega_qperp[pos_array_p_over_T][pos_array_omega_over_T+1][pos_array_qperp_over];
	const double v011=differential_rate_p_omega_qperp[pos_array_p_over_T][pos_array_omega_over_T+1][pos_array_qperp_over+1];
	const double v100=differential_rate_p_omega_qperp[pos_array_p_over_T+1][pos_array_omega_over_T][pos_array_qperp_over];
	const double v101=differential_rate_p_omega_qperp[pos_array_p_over_T+1][pos_array_omega_over_T][pos_array_qperp_over+1];
	const double v110=differential_rate_p_omega_qperp[pos_array_p_over_T+1][pos_array_omega_over_T+1][pos_array_qperp_over];
	const double v111=differential_rate_p_omega_qperp[pos_array_p_over_T+1][pos_array_omega_over_T+1][pos_array_qperp_over+1];

	//fraction of each corner to use
	const double frac_p_over_T=(p_over_T-p_over_T_low_val)/(p_over_T_high_val-p_over_T_low_val);
	const double frac_omega_over_T=(omega_over_T-omega_over_T_low_val)/(omega_over_T_high_val-omega_over_T_low_val);
	const double frac_qperp_over_T=(qperp_over_T-qperp_over_T_low_val)/(qperp_over_T_high_val-qperp_over_T_low_val);

	//get the value
	const double v00=v000*(1-frac_qperp_over_T)+v001*frac_qperp_over_T;
	const double v01=v010*(1-frac_qperp_over_T)+v011*frac_qperp_over_T;
	const double v10=v100*(1-frac_qperp_over_T)+v101*frac_qperp_over_T;
	const double v11=v110*(1-frac_qperp_over_T)+v111*frac_qperp_over_T;

	const double v0=v00*(1-frac_omega_over_T)+v01*frac_omega_over_T;
	const double v1=v10*(1-frac_omega_over_T)+v11*frac_omega_over_T;

	double res=v0*(1-frac_p_over_T)+v1*frac_p_over_T;

	return res;

}

////Rate
//double ERateColl::rate(const struct ERateParam &rate_params, const Pythia8::Vec4 &p0, const int &id0)
double Tequila::rate(double energy, double temp, double * rate_p)
{

//	const double temp=rate_params.T();
//
//	if (id0!=fId0) { 
//		// the incoming particle is not a gluon and this rate doesn't apply
//		return 0. ;
//	} 
//	if (temp/p0.e() >  ERateParam::MaxTOverP ) {
//		return 0. ;  // Don't evolve if the energy is too small
//	}

	//energy/temperature ratio of the incoming jet
	//const double energy_over_T=p0.e()/temp;
	const double energy_over_T=energy/temp;

	//if outside of tabulated range, set to 0 (that is, assume untabulated rate is tiny and can be ignored)
	if ((energy_over_T<p_over_T_min())||(energy_over_T>p_over_T_max())) return 0.0;
	//get the real-valued index
	double a = get_index_from_p_over_T(energy_over_T);
	//get the actual index of the array
	int n_p = int(a);
	//get the remainder
	a -= n_p;
	//get rate for energy/T
	double result = (1-a) * rate_p[n_p] + a * rate_p[n_p+1];

	//a factor of temperature is missing from the integrated rate
	result*=temp;


	return result;

}

//Given a value of parton momentum p, sample the rate in omega and qperp
void Tequila::sample_dgamma_dwdq(double p, double T, double *** differential_rate_p_omega_qperp, double &w, double &q) {
	//   double lam = pp.lambda(p0) ;
	int ntry = 0  ;
	const int ntry_max = 10000;   
	
	//helper variables
	//const double qperp_over_T_val_min=0.0;
	//const double qperp_over_T_val_max=qperp_over_T_max();
	const double p_over_T=p/T;
	const double omega_over_T_neg_min=omega_over_T_min(p_over_T);
	const double omega_over_T_neg_max=-omega_over_T_cutoff;
	const double omega_over_T_pos_min=omega_over_T_cutoff;
	const double omega_over_T_pos_max=p_over_T/2.0;
	const double omega_over_T_neg_span=(omega_over_T_neg_max-omega_over_T_neg_min);
	const double omega_over_T_pos_span=(omega_over_T_pos_max-omega_over_T_pos_min);

	const double max_rate=maximum_rate_p(p_over_T, differential_rate_p_omega_qperp);

	while (ntry < ntry_max) {
		//      fRateTable.sample_in_x1range(-GSL_DBL_MAX,lam/(2.*T), rng, params, w, q) ;
		//double r = max_rate*rng(params) ;
		double r = max_rate*ZeroOneDistribution(*GetMt19937Generator()) ;
		//      double v = (1. -  T*w/p0.e())/rratio(w) ;
		//const double qx_over_T_tmp=qperp_over_T_val_max*rng(params);
		//const double qy_over_T_tmp=qperp_over_T_val_max*rng(params);
		//const double q_over_T_test=sqrt(qx_over_T_tmp*qx_over_T_tmp+qy_over_T_tmp*qy_over_T_tmp);

		// sample omega_over_T
		// a bit more complicated due to the omega_over_T_cutoff around 0
		// With[{tmp = (rnd*(xnegspan + xposspan))},
		// If[tmp > xnegspan, xposmin + (tmp - xnegspan), xnegmin + tmp]
		//  ]
		//const double rnd2=rng(params);
		const double rnd2=ZeroOneDistribution(*GetMt19937Generator());
		const double tmp=rnd2*(omega_over_T_neg_span+omega_over_T_pos_span);
		double omega_over_T_test;
		if (tmp < omega_over_T_neg_span) omega_over_T_test=omega_over_T_neg_min+tmp;
		else omega_over_T_test=omega_over_T_pos_min+(tmp-omega_over_T_neg_span);

		// sample q_over_T, remembering it is a radial variable
		// also, q_over_T must be smaller than abs(omega_over_T_test)
		//const double q_over_T_sqr=omega_over_T_test*omega_over_T_test*rng(params);
		//const double q_over_T_test=sqrt(q_over_T_sqr);
		const double q_over_T_test=0.0;  //let's use qperp/T=0 for now

		double v = differential_rate(p_over_T, omega_over_T_test, q_over_T_test, differential_rate_p_omega_qperp);

		//safety check: did we pick the maximum rate correctly??
		if (v > max_rate) {
			std::cout << "Function \"maximum_rate_p()\" apparently does not return the maximum of the rate... This is bad. Fix it.\n";
			std::cout << "current guess for maximum=" << max_rate << " & sampled rate at omega/T="<< omega_over_T_test << "&q_perp/T=" << q_over_T_test << " is " << v << "\n";
			exit(1);
			//assert(false);
		}

		if (r < v) {
			w=omega_over_T_test*T;
			q=q_over_T_test*T;
			//std::cout << p_over_T << " " << w << " " << q << "\n";
			return ;
		}
		else {
			ntry++ ;
		}
	}
	w=omega_over_T_cutoff*T;
	q=0.0;
	std::cout << "*** ERateColl::sample_dgamma_dwdq *** Failed to find a sample "
		"after "<< ntry  << " iterations! Returning with w,q = " << w<<","<< q << std::endl;
}

//RateRadiative Tequila::getRateRadTotal(double pRest, double T)
//{
//  RateRadiative rate;
//
//  double u = pRest/T;  // making arguments in log to be dimensionless
//
//  rate.qqg = (0.8616 - 3.2913/(u*u) + 2.1102/u - 0.9485/sqrt(u))*pow(g, 4.)*T;
//  rate.ggg = (1.9463 + 61.7856/(u*u*u) - 30.7877/(u*u) + 8.0409/u - 2.6249/sqrt(u))
//    *pow(g, 4.)*T;
//  rate.gqq = (2.5830/(u*u*u) - 1.7010/(u*u) + 1.4977/u - 1.1961/pow(u,0.8) + 0.1807/sqrt(u))
//    *pow(g, 4.)*T*nf; 
//
//  double runningFactor = log(g*T*pow(10., 0.25)/.175)/log(g*T*pow(u, 0.25)/.175);
//  if (runningFactor < 1.)
//    {
//      rate.qqg *= runningFactor;
//      rate.gqq *= runningFactor;
//      rate.ggg *= runningFactor;
//    }
//
//  return rate;
//}
//
//RateRadiative Tequila::getRateRadPos(double u, double T)
//{
//  RateRadiative rate;
//
//  rate.qqg = (0.5322 - 3.1037/(u*u) + 2.0139/u - 0.9417/sqrt(u))*pow(g, 4.)*T;
//  rate.ggg = (1.1923 - 11.5250/(u*u*u) + 3.3010/u - 1.9049/sqrt(u))*pow(g, 4.)*T;
//  rate.gqq = (0.0004656 - 0.04621/(u*u) + 0.0999/u - 0.08171/pow(u,0.8) 
//              + 0.008090/pow(u,0.2) - 1.2525*pow(10.,-8.)*u)*pow(g, 4.)*T*nf; 
//
//  return rate;
//}
//
//RateRadiative Tequila::getRateRadNeg(double u, double T)
//{
//  RateRadiative rate;
//
//  rate.qqg = (0.3292 - 0.6759/(u*u) + 0.4871/pow(u,1.5) - 0.05393/u + 0.007878/sqrt(u))
//    *pow(g, 4.)*T;
//  rate.ggg = (0.7409 + 1.8608/(u*u*u) - 0.1353/(u*u) + 0.1401/u)*pow(g, 4.)*T;
//  rate.gqq = (0.03215/(u*u*u) + 0.01419/(u*u) + 0.004338/u - 0.00001246/sqrt(u))
//    *pow(g, 4.)*T*nf;
//
//  return rate;
//}
//
//double Tequila::getNewMomentumRad(double pRest, double T, int process)
//{
//  double newp = 0.;
//  double randA;
//  double x, y;
//  double fy, fyAct;
//
//  RateRadiative Pos, Neg;
//  double u = pRest/T;  // making arguments in log to be dimensionless
//
//  Pos = getRateRadPos(u, T);
//  Neg = getRateRadNeg(u, T);
//
//  // this switch will hold the decision whether k is positive or negative:
//  // 0 : negative, 1 : positive
//  int posNegSwitch = 1; 
//
//  /* process == 1 : quark radiating gluon
//     process == 2 : quark radiating photon
//     process == 3 : gluon radiating gluon
//     process == 4 : gluon split into quark-antiquark pair */
//
//  if (process == 1)
//    {
//      // decide whether k shall be positive or negative 
//      // if x (uniform on [0,1]) < area(k<0)/area(all k) then k < 0
//      if (ZeroOneDistribution(*GetMt19937Generator()) < Neg.qqg/(Neg.qqg+Pos.qqg))
//	posNegSwitch = 0;
//
//      if (posNegSwitch == 1) // if k > 0
//	{
//	  do
//	    {
//	      //randA is a uniform random number on [0, Area under the envelop function]
//	      randA = ZeroOneDistribution(*GetMt19937Generator())*area(u+12., u, posNegSwitch, 1);
//	      y = 2.5/(LambertW2(2.59235*pow(10.,23.)*exp(-100.*randA)));
//
//	      fy = 0.025/(y*y)+0.01/y;             // total area under the envelop function
//	      fyAct = function(u, y, process);     // actual rate
//
//	      x = ZeroOneDistribution(*GetMt19937Generator());       // random number, uniform on [0,1]
//
//	    } while (x > fyAct/fy); 
//	  // reject if x is larger than the ratio fyAct/fy
//	  newp = y;
//	}
//      else // if k < 0
//	{
//	  do
//	    {
//	      //randA is a uniform random number on [0, Area under the envelop function]
//	      randA = ZeroOneDistribution(*GetMt19937Generator())*area(-0.05, u, posNegSwitch, 1);
//	      y = -12./(1.+480.*randA);
//
//	      fy = 0.025/(y*y);                    // total area under the envelop function
//	      fyAct = function(u, y, process);     // actual rate
//
//	      x = ZeroOneDistribution(*GetMt19937Generator());       // random number, uniform on [0,1]
//
//	    } while (x > fyAct/fy); 
//	  // reject if x is larger than the ratio fyAct/fy
//	  newp = y;
//	}
//    }
//  else if (process == 2)
//    {
//      do
//	{
//	  randA = ZeroOneDistribution(*GetMt19937Generator())*area(1.15*u, u, posNegSwitch, 2);
//	  y = 83895.3*pow(pow(u, 0.5)*randA, 10./3.);
//
//	  fy = (0.01/(pow(y, 0.7)))/pow(u, 0.5); // total area under the envelop function
//	  fyAct = function(u, y, process);       // actual rate
//
//	  x = ZeroOneDistribution(*GetMt19937Generator());         // random number, uniform on [0,1]
//
//	} while (x > fyAct/fy); 
//      // reject if x is larger than the ratio fyAct/fy
//      newp = y;
//    }
//  else if (process == 3)
//    {
//      // decide whether k shall be positive or negative 
//      // if x (uniform on [0,1]) < area(k<0)/area(all k) then k < 0
//      if (ZeroOneDistribution(*GetMt19937Generator()) < Neg.ggg/(Neg.ggg+Pos.ggg))
//	posNegSwitch = 0;
//
//      if( posNegSwitch == 1 ) // if k > 0
//	{
//	  do
//	    {
//	      //randA is a uniform random number on [0, Area under the envelop function]
//	      randA = ZeroOneDistribution(*GetMt19937Generator())*area(u/2., u, posNegSwitch, 3);
//	      y = 5./(LambertW2(2.68812*pow(10., 45.)*exp(-50.*randA)));
//
//	      fy = 0.1/(y*y)+0.02/y;               // total area under the envelop function
//	      fyAct = function(u, y, process);     // actual rate
//
//	      x = ZeroOneDistribution(*GetMt19937Generator());       // random number, uniform on [0,1]
//
//	    } while (x > fyAct/fy); 
//	  // reject if x is larger than the ratio fyAct/fy
//	  newp = y;
//	}
//      else // if k < 0
//	{
//	  do
//	    {
//	      //randA is a uniform random number on [0, Area under the envelop function]
//	      randA = ZeroOneDistribution(*GetMt19937Generator())*area(-0.05, u, posNegSwitch, 3);
//	      y = -12./(1. + 120.*randA);
//
//	      fy = 0.1/(y*y);                      // total area under the envelop function
//	      fyAct = function(u, y, process);     // actual rate
//
//	      x = ZeroOneDistribution(*GetMt19937Generator());       // random number, uniform on [0,1]
//
//	    } while(x > fyAct/fy); 
//	  // reject if x is larger than the ratio fyAct/fy
//	  newp = y;
//	}
//    }
//  else if (process == 4)
//    {
//      // decide whether k shall be positive or negative 
//      // if x (uniform on [0,1]) < area(k<0)/area(all k) then k < 0
//      if (ZeroOneDistribution(*GetMt19937Generator()) < Neg.gqq/(Neg.gqq+Pos.gqq))
//	posNegSwitch = 0;
//
//      if( posNegSwitch == 1 ) // if k > 0
//	{
//	  do
//	    {
//	      //randA is a uniform random number on [0, Area under the envelop function]
//	      randA = ZeroOneDistribution(*GetMt19937Generator())*area(u/2., u, posNegSwitch, 4);
//	      y = 0.83333*(0.06*function(u, 0.05, process)+randA)/function(u, 0.05, process);
//
//	      fy = 1.2*function(u, 0.05, process); // total area under the envelop function
//	      fyAct = function(u, y, process);     // actual rate
//
//	      x = ZeroOneDistribution(*GetMt19937Generator());       // random number, uniform on [0,1]
//
//	    } while (x > fyAct/fy); 
//	  // reject if x is larger than the ratio fyAct/fy
//	  newp = y;
//	}
//      else // if k < 0
//	{
//	  do
//	    {
//	      //randA is a uniform random number on [0, Area under the envelop function]
//	      randA = ZeroOneDistribution(*GetMt19937Generator())*area(-0.05, u, posNegSwitch, 4);
//	      y = (2.5-u*log(7.81082*pow(10., -6.)*exp(14.5/u)+(-115.883+113.566*u)*randA))/(1.-0.98*u);
//
//	      fy = 0.98*exp((1.-1./u)*(-2.5+y))/u; // total area under the envelop function
//	      fyAct = function(u, y, process);     // actual rate
//
//	      x = ZeroOneDistribution(*GetMt19937Generator());       // random number, uniform on [0,1]
//
//	    } while (x > fyAct/fy); 
//	  // reject if x is larger than the ratio fyAct/fy
//	  newp = y;
//	}
//    }
//  else
//    {
//      WARN << "Invalid process number (" << process << ")";
//    }
//
//  return newp;
//}
//
//// calculates the area under the envelope function when using the rejection method
//// integrals had been solved analytically before
//double Tequila::area(double y, double u, int posNegSwitch, int process)
//{
//  if (process == 1)
//    {
//      if (posNegSwitch == 1)
//	return (0.5299 - 0.025/y + 0.01*log(y));
//      else 
//	return (-0.002083-0.025/y);
//    }
//  else if (process == 2)
//    {
//      return ((0.03333*pow(y,0.3))/pow(u,0.5));
//    }
//  else if (process == 3)
//    {
//      if (posNegSwitch == 1)
//	return (2.05991 - 0.1/y + 0.02*log(y));
//      else 
//	return (-0.008333 - 0.1/y);
//    }      
//  else if (process == 4)
//    {
//      if (posNegSwitch == 1)
//	return (1.2*function(u, 0.05, process)*(y-0.05));
//      else 
//	return ((6.8778*pow(10., -8.)*exp(14.5/u)
//		 -0.008805*exp((2.5-y+0.98*u*y)/u))/(1.0204-u));
//    }      
//
//  return 0.;
//}
//
//double Tequila::function(double u, double y, int process)
//{
//  if (process == 1)      return getRate_qqg(u, y);
//  else if (process == 3) return getRate_ggg(u, y);
//  else if (process == 4) return getRate_gqq(u, y);
//
//  return 0.;
//}
//
//
//// Reads in the binary stored file of dGamma values
//void Tequila::readRadiativeRate(Gamma_info *dat, dGammas *Gam)
//{
//  FILE *rfile;
//  string filename;
//  filename = PathToTables+"radgamma";
//
//  INFO << "Reading rates of inelastic collisions from file ";
//  INFO << filename.c_str() << " ... ";
//  size_t bytes_read;
//
//  rfile = fopen(filename.c_str(), "rb"); 
//  bytes_read = fread((char *)(&dat->ddf), sizeof(double), 1, rfile);
//  bytes_read = fread((char *)(&dat->dda), sizeof(double), 1, rfile);
//  bytes_read = fread((char *)(&dat->dcf), sizeof(double), 1, rfile);
//  bytes_read = fread((char *)(&dat->dca), sizeof(double), 1, rfile);
//  bytes_read = fread((char *)(&dat->Nc), sizeof(int), 1, rfile);
//  bytes_read = fread((char *)(&dat->Nf), sizeof(int), 1, rfile);
//  bytes_read = fread((char *)(&dat->BetheHeitler),sizeof(int) , 1, rfile);
//  bytes_read = fread((char *)(&dat->BDMPS), sizeof(int), 1, rfile);
//  bytes_read = fread((char *)(&dat->include_gluons), sizeof(int), 1, rfile);
//  bytes_read = fread((char *)Gam->qqg, sizeof(double), NP*NK, rfile);
//  bytes_read = fread((char *)Gam->tau_qqg, sizeof(double), NP*NK, rfile);
//  bytes_read = fread((char *)Gam->gqq, sizeof(double), NP*NK , rfile);
//  bytes_read = fread((char *)Gam->tau_gqq, sizeof(double), NP*NK, rfile);
//  bytes_read = fread((char *)Gam->ggg, sizeof(double), NP*NK, rfile);
//  bytes_read = fread((char *)Gam->tau_ggg, sizeof(double), NP*NK, rfile);
//  bytes_read = fread((char *)Gam->qqgamma, sizeof(double), NP*NK, rfile);
//  bytes_read = fread((char *)Gam->tau_qqgamma, sizeof(double), NP*NK, rfile);
//  fclose (rfile);
//
//  dat->Nf = nf;
//  dat->dp = 0.05;
//  dat->p_max = 20;
//  dat->p_min = 0;                 // exp(LogEmin) set to zero because my array starts at zero!...
//  dat->n_p = static_cast<int>(1.001+dat->p_max/dat->dp);    // np = int(0.4 + 121 / 0.5) = 241
//  dat->p_max = dat->dp*dat->n_p;                            // p_max = 0.5*241 = 120.5
//  dat->n_pmin = static_cast<int>(1.001+dat->p_min/dat->dp); // np_min = int(0.4 + 3.3 / 0.5) = 7
//  dat->n_p -= dat->n_pmin-1;                                // n_p = 241 - (7 - 1) = 235
//  dat->p_min = dat->dp * dat->n_pmin;                       // p_min = 0.5 * 7 = 3.5
//  dat->n_kmin = 1+2*(static_cast<int>(2./dat->dp));
//  dat->k_min = -dat->dp*dat->n_kmin;
//  dat->n_k = static_cast<int>((8+dat->p_max)/(2*dat->dp));
//  dat->k_max = 2*dat->dp*(dat->n_k-1)+dat->k_min;
//}
//
//
//
//double Tequila::getRate_qqg(double p, double k)
//{
//  return use_table(p, k, Gam.qqg, 0);
//}
//
//double Tequila::getRate_gqq(double p, double k)
//{
//  if (k < p/2.) return use_table(p, k, Gam.gqq, 1);
//  else return 0.;
//}
//
//double Tequila::getRate_ggg(double p, double k)
//{
//  if (k < p/2.) return use_table(p, k, Gam.ggg, 2);
//  else return 0.;
//}
//
//double Tequila::use_table(double p, double k, double dGamma[NP][NK], int which_kind)
///* Uses the lookup table and simple interpolation to get the value
//   of dGamma/dkdx at some value of p,k.
//   This works by inverting the relations between (p,k) and (n_p,n_k)
//   used in building the table, to find out what continuous values
//   of n_p, n_k should be considered; then linearly interpolates.     */
//{
//  double a, b, result;     // fraction of way from corner of box
//  int    n_p, n_k;         // location of corner of box
//
//  // out of range
//  if ((p < 4.01) || (p > 46000.) || (k < -12.) || (k > p+12.))
//    return 0.;
//
//  if ((which_kind % 3) && (k > p/2))
//    k = p - k;  // Take advantage of symmetry in these cases
//
//  a = 24.7743737154026*log(p*0.2493765586034912718l);
//  n_p = (int)a;
//  a -= n_p;
//  if (k < 2.)
//    {
//      if (k < -1)
//	{
//	  if (k < -2) b = 60.+5.*k;
//	  else b = 70.+10.*k;
//	}
//      else
//	{
//	  if (k < 1.) b = 80. + 20.*k;
//	  else b = 90.+10.*k;
//	}
//    }
//  else if ( k < p-2. )
//    { /* This is that tricky middle ground. */
//      b = 190.-10.*log(1.000670700260932956l/ 
//		       (0.0003353501304664781l+(k-2.)/(p-4.))-1.);
//    }
//  else
//    {
//      if (k < p+1.)
//	{
//	  if (k < p-1.) b = 290. + 10.*(k-p);
//	  else  b = 300. + 20.*(k-p);
//	}
//      else
//	{
//	  if (k < p+2.) b = 310. + 10.*(k-p);
//	  else b = 320. + 5.*(k-p);
//	}
//    }
//
//  n_k = (int)b;
//  b -= n_k;
//  result = (1.-a)*((1.-b)*dGamma[n_p][n_k]+b*dGamma[n_p][n_k+1])
//    +a*((1.-b)*dGamma[n_p+1][n_k]+b*dGamma[n_p+1][n_k+1]);
//
//  if (std::abs(k) > 0.001) // Avoid division by 0, should never get asked for
//    {
//      switch (which_kind)
//	{
//	case 0:
//	  result /= k;
//	  if (k < 20.)
//	    result /= 1.-exp(-k);
//	  if (k > p-20.)
//	    result /= 1. + exp(k-p);
//	  break;
//	case 1:
//	  result /= p;
//	  if (k < 20.)
//	    result /= 1 + exp(-k);
//	  if (k > p-20.)
//	    result /= 1. + exp(k-p);
//	  break;
//	case 2:
//	  result /= k*(p-k)/p;
//	  if (k < 20.)
//	    result /= 1.-exp(-k);
//	  if (k > p-20.)
//	    result /= 1.-exp(k-p);
//	  break;
//	case 3:
//	  result /= k;
//	  if (k < 0) result = 0.;
//	  if (k > p-20.)
//	    result /= 1. + exp(k-p);
//	  break;
//	}
//    }
//
//  return result;
//}
//
//
//double LambertW2(double z)
//{
//  double w_new, w_old, ratio, e_old, tol;
//  int n;
//
//  tol = 1.0e-14;
//
//  if(z <= -exp(-1.0))
//    {
//      WARN << "LambertW is not defined for z = " << z;
//      WARN << "z needs to be bigger than " << -exp(-1.0);
//      throw std::runtime_error("LambertW small z problem");
//    }
//
//  if(z > 3.0)
//    {
//      w_old = log(z) - log(log(z));
//    }
//  else {w_old = 1.0;}
// 
//  w_new = 0.0;
//  ratio = 1.0;
//  n = 0;
//  while(std::abs(ratio) > tol) 
//    {
//      e_old = exp(w_old);
//      ratio = w_old*e_old - z;
//      ratio /= ( e_old*(w_old + 1.0) - (w_old+2.0)*(w_old*e_old-z)/(2.0*w_old + 2.0) );
//      w_new = w_old - ratio;
//      w_old = w_new;
//      n++;
//      if(n > 99) 
//	{
//          WARN << "LambertW is not converging after 100 iterations.";
//          WARN << "LambertW: z = " << z;
//          WARN << "LambertW: w_old = " << w_old;
//          WARN << "LambertW: w_new = " << w_new;
//          WARN << "LambertW: ratio = " << ratio;
//          throw std::runtime_error("LambertW not conversing");
//	}
//    }
//
//  return w_new;
//}// LambertW by Sangyong Jeon

