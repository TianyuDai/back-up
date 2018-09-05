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

#ifndef MARTINI_H
#define MARTINI_H

#include <fstream>
#include <math.h>
#include "JetEnergyLossModule.h"
#include "JetScapeConstants.h"
#include "Tequila/evolve9_text_noExp.h"
#include <gsl//gsl_spline.h>

using namespace Jetscape;

//Basic.h//
struct RateRadiative
{
  double qqg;
  double ggg;
  double gqq;
};

enum radiative_process_type { qqg, ggg, gqq, none};

class Tequila : public JetEnergyLossModule<Tequila> //, public std::enable_shared_from_this<Tequila>
{  
 private:

  // AMY rates are calculated in p/T > AMYpCut
  static constexpr double AMYpCut = 4.01;

  double Q0;
  double alpha_s;
  double alpha_em;
  double g;
  double pcut;        // below this scale, no further Eloss
  double hydro_Tc;
  double omega_over_T_cutoff;


  // size of the array containing the differential rate in p, omega and q_\perp
  // those number in p and omega are, and must be, the same used in Guy Moore's code that is used to evaluate the collinear rates
  // the number of point in q_\perp is arbitrary at this point
  static const int nb_points_in_p=NP, nb_points_in_omega=NK;
  // Addionnal grid in q_perp, not defined in Guy Moore's code
  static const int nb_points_in_qperp=4;
  // the minimum and maximum parton energy p and parton energy loss omega are set in Guy Moore's code as well and must be the same here
  double p_over_T_min() { return 4.01; };
  double p_over_T_max() { return (p_over_T_min()*pow(1.04119,nb_points_in_p-1)); };
  double omega_over_T_min(const double p_over_T) { return -12.0; };
  double omega_over_T_max(const double p_over_T) { return p_over_T + 0.2 * (nb_points_in_omega -1 - 320); };
  // for future use
  //double qperp_over_T_min(const double p_over_T, const double omega_over_T) { return -5.0; };
  //double qperp_over_T_min() { return -10.0; };
  double qperp_over_T_max() { return 10.0; };
  //given the above tabulation of p/T, this function returns the index for a given p/T
  //the index is a real number, such that "int(index)" gives the position in the array and "index-int(index)" gives the remainder
  double get_index_from_p_over_T(const double p_over_T) { return 24.7743737154026 * log( p_over_T * .2493765586034912718l ); };
  double get_index_from_omega_over_T(const double p_over_T, const double omega_over_T); //Big function, defined elsewhere instead
  double get_index_from_qperp_over_T(const double qperp_over_T) { 
	  const double qmin=0.0;
	  const double qmax=qperp_over_T_max();
	  return (((nb_points_in_qperp-1)*(qperp_over_T-qmin))/(qmax-qmin));
  }
  double get_p_over_T_from_index(const int index_p_over_T) { return (p_over_T_min()*pow(1.04119,index_p_over_T)); };
  double get_omega_over_T_from_index(const int index_p_over_T, const int index_omega_over_T); //Big function, defined elsewhere instead
  // Assume uniform grid in q_perp for now
  double get_qperp_over_T_from_index(const int index_qperp_over_T) { 
	  const double qmin=0.0;
	  const double qmax=qperp_over_T_max();
	  return qmin+(qmax-qmin)/(nb_points_in_qperp-1)*index_qperp_over_T;
  };
  // arrays containing the differentail and integrated rates
  //double differential_rate_p_omega_qperp[nb_points_in_p][nb_points_in_omega][nb_points_in_qperp];
  //double rate_p[nb_points_in_p];
  double *** differential_rate_qqg_p_omega_qperp, *** differential_rate_ggg_p_omega_qperp, *** differential_rate_gqq_p_omega_qperp;
  double *rate_qqg_p, * rate_ggg_p, *rate_gqq_p;
  //maximum of the differential rate, used to sample the rate

  void allocate_memory_for_radiative_rate_table();

  double non_uniform_trapeze_weights(int position, int size_of_array, double prev_x, double curr_x, double next_x);

  // load differential rate from file into member array "differential_rate_p_omega_q"
  void load_differential_rate(const double alpha_s, const double alpha_EM, const int Nf, const std::string location_of_collinear_rates);
  // fill member array "rate_p" by integrating the content of "differential_rate_p_omega_q"
  //void evaluate_integrated_rate(double omegaMin);
  void evaluate_integrated_rate(double omega_over_T_cutoff, double *** differential_rate_gqq_p_omega_qperp, double * rate_gqq_p);

//  double differential_rate(const double p_over_T, const double omega_over_T, const double qperp_over_T);
  double differential_rate(const double p_over_T, const double omega_over_T, const double qperp_over_T, double *** differential_rate_p_omega_qperp);
  double rate(double energy, double temp, double * rate_p);
  void sample_dgamma_dwdq(double p, double T, double *** differential_rate_p_omega_qperp, double &w, double &q);
//  void sample_dgamma_dwdq(const struct ERateParam &pp, const Pythia8::Vec4 &p0, double (*rng)(void*params), void *params,  double &w, double &q);
//
  //educated guess on the highest value of the rate
  //double *maximum_differential_rate;
  double maximum_rate_p(double p_over_T, double *** differential_rate_p_omega_qperp) {
	  //the collinear rate should presumably be maximum at qperp/T=0 and either omega/T=+/-omega_over_T_cut 
	  const double val1=differential_rate(p_over_T,omega_over_T_cutoff,0,differential_rate_p_omega_qperp);
	  const double val2=differential_rate(p_over_T,-1*omega_over_T_cutoff,0,differential_rate_p_omega_qperp);
	  return val1 > val2 ? val1 : val2;
  };


  //Import.h//
//  static const int NP = 230;
//  static const int NK = 381;

//  static const int Nalphas = 11;
//  static const int Nomega = 120;
//  static const int Nq = 60;
//
//  static constexpr double omegaStep = 0.2;
//  static constexpr double qStep = 0.2;
//  static constexpr double alphaMin = 0.15;
//  static constexpr double alphaStep = 0.03;

//  typedef struct
//  {
//    double ddf;
//    double dda;
//    double dcf;
//    double dca;
//    int    include_gluons;
//    int Nc;
//    int Nf;
//    int BetheHeitler;
//    int BDMPS;
//    double alpha;
//    double delta_x;
//    double dx_max;
//    double dp;
//    double p_min;
//    double p_max;
//    long   n_p;
//    long   n_pmin;
//    double k_min;
//    double k_max;
//    long   n_k;
//    long   n_kmin;
//  } Gamma_info;

//  // Structure to store information about splitting functions
//  typedef struct
//  {
//    double qqg[NP][NK];
//    double gqq[NP][NK];
//    double ggg[NP][NK];
//    double qqgamma[NP][NK];
//
//    double tau_qqg[NP][NK];
//    double tau_gqq[NP][NK];
//    double tau_ggg[NP][NK];
//    double tau_qqgamma[NP][NK];
//  } dGammas;

//  Gamma_info dat;
//  dGammas    Gam;

//  vector<double> *dGamma_qq;
//  vector<double> *dGamma_qg;
//  vector<double> *dGamma_qq_q;
//  vector<double> *dGamma_qg_q;

 public:
  
  Tequila();
  virtual ~Tequila();

  //main//
  void Init();
  void DoEnergyLoss(double deltaT, double Time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
  //int DetermineProcess(double p, double T, double deltaT, int id);
  enum radiative_process_type DetermineProcess(double pRest, double T, double deltaT, int Id);
  void WriteTask(weak_ptr<JetScapeWriter> w) {};
  
//  //Radiative.h//
//  RateRadiative getRateRadTotal(double p, double T);
//  RateRadiative getRateRadPos(double u, double T);
//  RateRadiative getRateRadNeg(double u, double T);

//  double getNewMomentumRad(double p, double T, int process);
//  double area(double y, double u, int posNegSwitch, int process);
//  double function(double u, double y, int process);

//  double areaOmega(double u, int posNegSwitch, int process);
//  double areaQ(double u, double omega, int process);
//  Jetscape::FourVector getNewMomentumElas(Jetscape::FourVector pVec, double omega, double q);

//  //Import.h//
//  void readRadiativeRate(Gamma_info *dat, dGammas *Gam);

//  double getRate_qqg(double p, double k);
//  double getRate_gqq(double p, double k);
//  double getRate_ggg(double p, double k);
//  double getRate_qqgamma(double p, double k);
//  double use_table(double p, double k, double dGamma[NP][NK], int which_kind);


 protected:
  uniform_real_distribution<double> ZeroOneDistribution;

  std::string path_to_tables;
};

//double LambertW2(double z);

#endif // MARTINI_H

