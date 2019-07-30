/* HMC.cpp
 * defines TODO
 */

// ========================================================================= //
// dependencies

#include <iostream>
#include <complex>

#include <vector>

#include "vectorStatistics.hpp"
#include "grid.hpp"
#include "HMC.hpp"


// ========================================================================= //
// CTor

void HMC::init(const int L, const double spinVariation) {
  if (L < 1) {
    std::cerr << "Invalid grid size." << std::endl;
    return;
  }
  
  this->G = Grid(L, spinVariation);
  this->V = this->G.getV();
  
  this->spinVariation = spinVariation;
  
  this->reset();
}
// ----------------------------------------------------------------------- //
HMC::HMC(const int L)                             {this->init    (L, 0);}
// ....................................................................... //
HMC::HMC(const int L, const Startcondition SC)    {this->init    (L, StartconditionToVariation(SC));}
// ....................................................................... //
HMC::HMC(const int L, const double spinVariation) {this->init    (L, spinVariation);}


// ========================================================================= //
// Grid quick interface

inline int    HMC::getNeighbourOf(const int i, const int j) const {return this->G.getNeighboursOf(i)[j];}

inline double HMC::getSpin     (const int i) const {return this->G.getQ()[i];}
inline double HMC::getTorque   (const int i) const {return this->G.getP()[i];}


// ========================================================================= //
// thermodynamic properties from grid

double HMC::getM() const {
  // magnetization
  // add up spins with their x/y components.
  // use complex numbers formalism to save on function calls.
  
  std::complex<double> z(0, 0);
  
  for (auto phi : this->G.getQ()) {z += std::exp<double>(I * phi);}
  
  return std::abs(z);
}
// ....................................................................... //
double HMC::getE() const {
  // Energy / Hamiltonian
  // add up projections of nearest forward neighbours
  
  double reVal = 0;
  double phi[3];    // spin orientation of self, x-neighbour and y-neighbour
  
  for (auto i=0; i<this->V; i++) {
    phi[0] = this->getSpin(i);
    phi[1] = this->getSpin(this->getNeighbourOf(i, 0));
    phi[2] = this->getSpin(this->getNeighbourOf(i, 1));
    
    reVal -= std::cos(phi[1] - phi[0]);
    reVal -= std::cos(phi[2] - phi[0]);
  }
  
  return reVal;
}
// ----------------------------------------------------------------------- //
double HMC::getMDensity() const {return this->getM() / this->V;}
// ....................................................................... //
double HMC::getEDensity() const {return this->getE() / this->V;}
// ----------------------------------------------------------------------- //
double HMC::getHamiltonian() const {
  double pot = this->getE() / this->T;
  double kin = 0;
  
  for (auto p : this->G.getP()) {kin += p * p;}
  kin /= 2;
  
  return (pot + kin) / (this->V);
}


// ========================================================================= //
// runtime params

double HMC::getSpinVariation() const         {return this->spinVariation;}
// ....................................................................... //
void   HMC::setSpinVariation(const double s) {if (s > 0) {this->spinVariation = s;}}
// ----------------------------------------------------------------------- //
double HMC::getT() const                     {return this->T;}
// ....................................................................... //
void   HMC::setT(const double T)             {if (T > 0) {this->T = T;}}
// ----------------------------------------------------------------------- //
double HMC::getDT() const                    {return this->dt;}
// ....................................................................... //
void   HMC::setDT(const double DT)           {if (DT > 0) {this->dt = DT;}}
// ----------------------------------------------------------------------- //
int    HMC::getSteps() const                 {return this->steps;}
// ....................................................................... //
void   HMC::setSteps(const int n)            {if (n > 0) {this->steps = n;}}


// ========================================================================= //
// runtime

void HMC::run(int N_MC) {
  std::cout << MEDSEP;
  std::cout << "Attempting to run simulation." << std::endl;
  std::cout << "\tT  = " << to_string_with_precision(this-> T   , 2) << std::endl;
  std::cout << "\tdt = " << to_string_with_precision(this->dt   , 2) << std::endl;
  std::cout << "\tn  = " << to_string_with_precision(this->steps, 0) << std::endl;
  
  this->reset();
  
  std::cout << "Ready." << std::endl;
  
  double beta = 1.0 / this->T;
  bool   accept = false;
  
  auto & G = this->G;
  
  std::vector<double> thisQ;
  std::vector<double> thisP;
  
  double lastE = this->getEDensity(),
         lastM = this->getMDensity();
  double lastH = this->getHamiltonian(),
         thisH = lastH;
  
  int nAcc = 0, nRej = 0;
  
  historyE.push_back(lastE);
  historyM.push_back(lastM);
  
  // ....................................................................... //
  auto lambda_makeP = [
    this,
    &thisQ, &thisP, 
    &G, &beta
  ] (double dt_factor = 1) {
    /* p[i](+1/2) = p(0) - dS(0)/dq[i]
      * S "=" beta * this.getE(): sum over cos(delta phi)
      * where delta phi = phi_i - phi_j
      * and j is forward neighbour of i
      * 
      * this gives two contributions in the derivative
      * - where d/q[i] adresses the central node
      * - where this adresses a neighbour node
      * 
      * 
      *      a     d/dq may refer to spins a, b, c.
      *      |     S = beta [... - cos (a - b) - cos(c - b) - ...]
      *   d..b--c
      *      :     since S is only forward-defined, b appears in the terms with
      *      e     centers b (2x), d, e (once, each) 
      * 
      * cos(x) = cos(-x)
      * d/dq -cos(t - q) = d/dq -cos(q - t) = +sin(q - t)
      * 
      * since p -= dS/dq, the final contribution is negative again
      */
    
    double phi, delta;
    
    for (auto i=0u; i<thisQ.size(); i++) {       // loop i: lattice site
      delta = 0;
      for (auto n : G.getNeighboursOf(i)) {
        phi = thisQ[i] - thisQ[n];
        delta += std::sin(phi);
      }
      delta *= beta * this->dt * dt_factor;   // dt_factor: allow for full or half time step
      thisP[i] -= delta;
    }
  };
  // ....................................................................... //
  
  for (auto m=0; m<N_MC; m++) {
    G.makeRandomP();
    thisQ = G.getQ();
    thisP = G.getP();
    
    lambda_makeP(0.5);                                                          // initial half step for momenta
    
    // move forward through time: n-1 steps
    for (auto j = 0; j < this->steps - 1; j++) {
      for (auto i=0u; i<thisQ.size(); i++) {thisQ[i] += thisP[i] * this->dt;}   // update spins with most recent
      lambda_makeP(1);                                                          // update momenta: same as before, but with full time step
    }
    
    for (auto i=0u; i<thisQ.size(); i++) {thisQ[i] += thisP[i] * this->dt;}     // final half step for spins
    lambda_makeP(0.5);                                                          // final half step for momenta (half-time step again)
    
    
    // the Hamiltonian method expects these to be in-place. Copy back lastQ on reject
    G.swapQ(thisQ);
    G.swapP(thisP);
    
    // accept/reject
    thisH  = this->getHamiltonian();
    
    if      (thisH < lastH)                                  {accept = true ;}
    else if (gsl_rng_uniform(RNG) < std::exp(lastH - thisH)) {accept = true ;}     // exp(-delta H) = exp(-(this - last)) = exp(last - this)
    else                                                     {accept = false;}
    
    if (accept) {
      lastH = thisH;
      lastE = this->getEDensity();
      lastM = this->getMDensity();
      nAcc++;
      
    } else {
      G.swapQ(thisQ);
      thisH = lastH;
      nRej++;
    }
    
    historyE.push_back(lastE);
    historyM.push_back(lastM);
    
    // show progress
    if ( !((m+1) % (N_MC / 100)) ) {std::cout << "." << std::flush;}
  }
  
  
  std::cout << std::endl << "Run finished." << std::endl;
  std::cout << nAcc << " accepted pts, " << nRej << " rejected." << std::endl;
}
// ....................................................................... //
void HMC::reset() {
  this->historyE.clear();
  this->historyM.clear();
  
  this->valE = NAN;
  this->valC = NAN;
  this->tauE = NAN;
  this->errE = NAN;
  this->errC = NAN;
  
  this->valM = NAN;
  this->valX = NAN;
  this->tauM = NAN;
  this->errM = NAN;
  this->errX = NAN;
  
  this->G.makeRandomQ(this->spinVariation);
}


// ========================================================================= //
// read results

const std::vector<double> & HMC::getHistoryE () const {return this->historyE;}
const std::vector<double> & HMC::getHistoryM () const {return this->historyM;}
// ----------------------------------------------------------------------- //
double HMC::getTauE          () {
  if (std::isnan(this->tauE)) {
    this->tauE = autocorrelationTime(this->historyE);
    if (std::isnan(this->tauE)) {return NAN;}
    
    if (20 * this->tauE > this->historyE.size()) {
      std::cerr << "Warning: Insufficient data to discard thermalization phase" << std::endl;
      return this->tauE;
    }
    
    this->tauE = autocorrelationTime(this->historyE, 20 * this->tauE);
    if (20 * this->tauE > this->historyE.size()) {
      std::cerr << "Warning: Insufficient data to discard thermalization phase" << std::endl;
    }
  }
  
  return this->tauE;
}
// ....................................................................... //
double HMC::getTauM          () {
  if (std::isnan(this->tauM)) {
    this->tauM = autocorrelationTime(this->historyM);
    if (std::isnan(this->tauM)) {return NAN;}
    
    if (20 * this->tauM > this->historyM.size()) {
      std::cerr << "Warning: Insufficient data to discard thermalization phase" << std::endl;
      return this->tauM;
    }
    
    this->tauM = autocorrelationTime(this->historyM, 20 * this->tauM);
    if (20 * this->tauM > this->historyM.size()) {
      std::cerr << "Warning: Insufficient data to discard thermalization phase" << std::endl;
    }
  }
  
  return this->tauM;
}
// ----------------------------------------------------------------------- //
double HMC::getValE          () {
  if (std::isnan(this->valE)) {
    this->valE = avg(this->historyE, 20 * this->getTauE());
  }
  return this->valE;
}
// ....................................................................... //
double HMC::getErrE          () {
  if (std::isnan(this->errE)) {this->errE = stderror(historyE, this->getTauE());}
  return this->errE;
}
// ----------------------------------------------------------------------- //
double HMC::getValM          () {
  if (std::isnan(this->valM)) {
    this->valM = std::abs(avg(this->historyM, 20 * this->getTauM()));
  }
  return this->valM;
}
// ....................................................................... //
double HMC::getErrM          () {
  if (std::isnan(this->errM)) {this->errM = stderror(historyM, this->getTauM());}
  return this->errM;
}
// ----------------------------------------------------------------------- //
double HMC::getValC          () {
  if (std::isnan(this->valC)) {
    this->valC = (1/(this->T*this->T)) * variance(this->historyE, 20 * this->getTauE()) * V;
  }
  return this->valC;
}
// ....................................................................... //
double HMC::getErrC          () {
  if (std::isnan(this->errC)) {
    this->errC = bootstrap_variance_error  (this->historyE, V / (T*T), 20 * this->tauE, this->tauE);
  }
  return this->errC;
}
// ----------------------------------------------------------------------- //
double HMC::getValX          () {
  if (std::isnan(this->valX)) {
    this->valX = (1/(    this->T    )) * variance(this->historyM, 20 * this->getTauM()) * V;
  }
  return this->valX;
}
// ....................................................................... //
double HMC::getErrX          () {
  if (std::isnan(this->errX)) {
    this->errX = bootstrap_variance_error  (this->historyM, V / T, 20 * this->tauM, this->tauM);
  }
  return this->errX;
}
