/* HMC.hpp
 * defines TODO
 */

#ifndef HMC_H
#define HMC_H

// ========================================================================= //
// dependencies

#include <vector>
#include <cmath>

#include "grid.hpp"

// ========================================================================= //
// class

class HMC {
private:
  Grid G = Grid(1);
  
  int V;
  
  double spinVariation = 0.0;
  double T             = 1.0;
  
  double dt    = 0.1;
  int    steps = 5;
  
  std::vector<double> historyE;
  std::vector<double> historyM;
  
  double tauE = NAN;
  double valE = NAN;
  double errE = NAN;
  double valC = NAN;
  double errC = NAN;
  
  double tauM = NAN;
  double valM = NAN;
  double errM = NAN;
  double valX = NAN;
  double errX = NAN;
  
  void init(const int L, const double spinVariation);
  
public:
  HMC(const int L);
  HMC(const int L, const Startcondition SC);
  HMC(const int L, const double spinVariation);
  
  // ....................................................................... //
  // Grid quick interface
  
  inline int    getNeighbourOf(const int i, const int j) const;
  
  inline double getSpin       (const int i) const;
  inline double getTorque     (const int i) const;
  
  // ....................................................................... //
  // thermodynamic properties from grid
  
  double getM() const;
  double getE() const;
  
  double getMDensity() const;
  double getEDensity() const;
  
  double getHamiltonian() const;
  
  // ....................................................................... //
  // runtime params
  
  double getSpinVariation() const;
  void   setSpinVariation(const double s);
  
  double getT() const;
  void   setT(const double T);
  
  double getDT() const;
  void   setDT(const double DT);
  
  int    getSteps() const;
  void   setSteps(const int n);
  
  // ....................................................................... //
  // runtime
  
  void run(int N_MC);
  void reset();
  
  // ....................................................................... //
  // read results
  
  const std::vector<double> & getHistoryE () const;
  const std::vector<double> & getHistoryM () const;
  
  double getTauE          ();
  double getTauM          ();
  
  double getValE          ();
  double getErrE          ();
  
  double getValM          ();
  double getErrM          ();
  
  double getValC          ();
  double getErrC          ();    // implements bootstrap
  
  double getValX          ();
  double getErrX          ();    // implements bootstrap
};

#endif//HMC_H
