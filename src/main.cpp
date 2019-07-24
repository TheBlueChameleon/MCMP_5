/* main.cpp
 * runs the HMC model simulation.
 */

// ========================================================================= //
// dependencies

#include <iostream>
#include <limits>

#include "globals.hpp"
#include "vectorStatistics.hpp"

#include "grid.hpp"
#include "HMC.hpp"

// ========================================================================= //
// dependencies

int main () {
  // ----------------------------------------------------------------------- //
  // initialize
  
  DEBUG_LAUNCH;
  RNG_init();
  
  // ----------------------------------------------------------------------- //
  // setup
  const int L    = 16;
  const int N_MC = 10000;
  
  const std::string fn_EMCX = "./out/EMCX.dat";
  std::ofstream     fh_EMCX(fn_EMCX);
  
  const std::string fn_cal = "./out/cal.dat";
  std::ofstream     fh_cal(fn_cal);
  
  fh_EMCX
      << "# "
      << "temperature T\t"
      << "heat density e\terror on e\t" 
      << "absolute magnetization density |m|\terror on m\t"
      << "specific heat density c\terror on c\t"
      << "magnetic susceptibility chi\terror on chi"
      << std::endl;
  
  fh_cal
      << "# n/dt calibration. First block is lo-T, second is hi-T\n"
      << "# n\tdt\t2 *n * tau\tn\tdt\t2 *n * tau"
      << std::endl;
  
  double dT = outerT_dT;
  
  std::cout << SEPARATOR;
  std::cout << "Setup with dimension L=" << L << std::endl;
  HMC model(L, Startcondition::COLDSTART);
  
  
  // ----------------------------------------------------------------------- //
  // calibration
  
  std::vector<int>    propN  = {   2,    4,   8,  16     };
  std::vector<double> propDT = {0.01, 0.05, 0.1, 0.2, 0.4};
  std::vector<double> valNTau( propN.size() * propDT.size() * 2 );    // holds all estimates for t_indep, for both temperatures -- hence the * 2.
  auto records = valNTau.size() / 2;
  
  auto lambda_recorder = [
    &model,
    &valNTau,
    &propN,
    &propDT
  ] (int offset = 0) {
    int i = 0;
    
    for   (auto  n : propN ) {
      for (auto dt : propDT) {
        model.setSteps( n);
        model.setDT   (dt);
        
        model.run(N_MC);
        
        valNTau[i + offset] = n * model.getTauM();
        i++;
      }
    }
  };
  
  
  std::cout << SEPARATOR << std::endl;
  std::cout << "Running low T calibration at T = "  << outerT_lo << std::endl;
  model.setT(outerT_lo);
  lambda_recorder();
  
  
  std::cout << SEPARATOR << std::endl;
  std::cout << "Running high T calibration at T = " << outerT_hi << std::endl;
  model.setT(outerT_hi);
  lambda_recorder(records);
  
  
  // file output
  for (auto i=0u; i<records; i++) {
    // lo T section
    fh_cal << propN  [ i / propDT.size() ] << "\t";      // "x-value": n
    fh_cal << propDT [ i % propN .size() ] << "\t";      // "y-value": dt
    fh_cal << valNTau[ i                 ] << "\t";      // "z-value": n * tau
    
    // hi T section
    fh_cal << propN  [ i / propDT.size() ] << "\t";      // "x-value": n
    fh_cal << propDT [ i % propN .size() ] << "\t";      // "y-value": dt
    fh_cal << valNTau[i        + records ] << "\t";      // "z-value": n * tau
    
    fh_cal << std::endl;
  }
  
  // find minima
  int    iLo = 0                                      , iHi = 0;
  double vLo = std::numeric_limits<double>::infinity(), vHi = std::numeric_limits<double>::infinity();
  for (auto i=0u; i<records; i++) {
    if (valNTau[i          ] < vLo) {vLo = valNTau[i          ]; iLo = i;}
    if (valNTau[i + records] < vHi) {vHi = valNTau[i + records]; iHi = i;}
  }
  const double dtLo = propDT[ iLo % propDT.size() ];
  const double dtHi = propDT[ iHi % propDT.size() ];
  const double  NLo = propN [ iLo % propN .size() ];
  const double  NHi = propN [ iHi % propN .size() ];
  
  std::cout <<   "Optimum for low  T: (n, dt) = (" << NLo << ", " << dtLo << ")" << std::endl;
  std::cout <<   "Optimum for high T: (n, dt) = (" << NHi << ", " << dtHi << ")" << std::endl;
  
  fh_cal    << "# Optimum for low  T: (n, dt) = (" << NLo << ", " << dtLo << ")" << std::endl;
  fh_cal    << "# Optimum for high T: (n, dt) = (" << NHi << ", " << dtHi << ")" << std::endl;
  
  // ----------------------------------------------------------------------- //
  // run
  
  std::cout << SEPARATOR << std::endl;
  std::cout << "About to begin productive run." << std::endl;
  
  double dt, Ns;
  for (auto T = outerT_lo; T <= outerT_hi; T += dT) {
    if (between(T, innerT_lo, innerT_hi)) {dT = innerT_dT;}
    else                                  {dT = outerT_dT;}
    
    dt = dtLo + (T - outerT_lo) / (outerT_hi - outerT_lo) * (dtHi - dtLo);
    Ns =  NLo + (T - outerT_lo) / (outerT_hi - outerT_lo) * ( NHi -  NLo);
    
    model.setDT   (dt);
    model.setSteps(Ns);
    
    model.setT(T);
    model.run(N_MC);
    
    std::cout << "computing primary and secondary quantities and their errors..." << std::flush;
    
    // file out
    fh_EMCX << T << "\t";
    
    if (std::isnan(model.getTauE()) || 
        std::isnan(model.getTauM())
    ) {
      fh_EMCX << "# --- insufficient data to compute meaningfull data ---" << std::endl;
      goto skipPointEMCX;
    }
    
    fh_EMCX << model.getValE() << "\t";
    fh_EMCX << model.getErrE() << "\t";
    
    fh_EMCX << model.getValM() << "\t";
    fh_EMCX << model.getErrM() << "\t";
    
    fh_EMCX << model.getValC() << "\t";
    fh_EMCX << model.getErrC() << "\t";
    
    fh_EMCX << model.getValX() << "\t";
    fh_EMCX << model.getErrX() << "\t";
    
    fh_EMCX << std::endl;
    
    
    skipPointEMCX:
    std::cout << "done" << std::endl;
  }
  // ----------------------------------------------------------------------- //
  // tidy up
  
  fh_EMCX.close();
  fh_cal  .close();
  
  std::cout << "simulation done." << std::endl;
  DEBUG_END;
} 
