/* grid.cpp
 * defines a grid of size LxL of spins with one continuous degree of freedom
 * (component q[i] in [0, 2pi]) as well as a prospective change of spin
 * (component p[i] in [0, 2pi]).
 */

// ========================================================================= //
// dependencies
 
#include <vector>
#include <array>

#include <cmath>
#include <iostream>

#include "globals.hpp"
#include "grid.hpp"

// ========================================================================= //
// CTor

void Grid::init(const int L, const double spinVariation) {
  if (L < 1) {
    std::cerr << "Invalid grid size." << std::endl;
    return;
  }
  
  this->L = L;
  this->V = L * L;
  
  this->Q         .resize(this->V, 0);
  this->P         .resize(this->V, 0);
  this->neighbours.resize(this->V);
  
  for (auto i=0; i<this->V; i++) {
    auto theNeighbours = this->neighboursOf(arrayCoordinateX(i), arrayCoordinateY(i));
    std::copy(theNeighbours.cbegin(), theNeighbours.cend(), std::back_inserter(this->neighbours[i]));
  }
  
  this->makeRandomQ(spinVariation);
}
// ----------------------------------------------------------------------- //
Grid::Grid(const int L)                             {this->init(L, StartconditionToVariation(Startcondition::COLDSTART));}
// ....................................................................... //
Grid::Grid(const int L, const Startcondition SC)    {this->init(L, StartconditionToVariation(SC));}
// ....................................................................... //
Grid::Grid(const int L, const double spinVariation) {this->init(L, spinVariation);}


// ========================================================================= //
// lattice geometry

int Grid::getL() const {return this->L;}
int Grid::getV() const {return this->V;}

// translate 1D array index into 2D coordinates
inline int Grid::arrayCoordinateX (const int idx) const {return idx % this->L;}
inline int Grid::arrayCoordinateY (const int idx) const {return idx / this->L;}

// translate 2D coordinates into 1D array index, incorporating periodic boundary conditions
inline int Grid::arrayIndex(const int x, const int y) const {
  return 
    (y >= static_cast<int>(this->L)  ?  y-this->L  :  (y < 0  ?  (this->L+y)  :  y)) * this->L + 
    (x >= static_cast<int>(this->L)  ?  x-this->L  :  (x < 0  ?  (this->L+x)  :  x));
}

// get automatic list of nearest neighbours
inline std::array<int, 4> Grid::neighboursOf (const int x, const int y) const {
  std::array<int, 4> reVal;
  
  reVal[0] = this->arrayIndex(x+1, y  );
  reVal[1] = this->arrayIndex(x  , y+1);
  reVal[2] = this->arrayIndex(x-1, y  );
  reVal[3] = this->arrayIndex(x  , y-1);
  
  return reVal;
}


// ========================================================================= //
// vector getters/setters

const std::vector<double> & Grid::getQ() const {return this->Q;}
const std::vector<double> & Grid::getP() const {return this->P;}
// ----------------------------------------------------------------------- //
void Grid::setP (const std::vector<double> &P) {
  if (static_cast<unsigned int>(this->V) != P.size()) {
    std::cerr << "Invalid size of momentum vector." << std::endl;
    return;
  }
  
  this->P = P;
}
// ....................................................................... //
void Grid::setQ (const std::vector<double> &Q) {
  if (static_cast<unsigned int>(this->V) != Q.size()) {
    std::cerr << "Invalid size of coordinate vector." << std::endl;
    return;
  }
  
  this->Q = Q;
}
// ----------------------------------------------------------------------- //
void Grid::swapQ(      std::vector<double> & Q) {this->Q.swap(Q);}
// ....................................................................... //
void Grid::swapP(      std::vector<double> & P) {this->P.swap(P);}
// ----------------------------------------------------------------------- //
void Grid::makeRandomP () {
  if (RNG_initialized) {for (auto & p : this->P) {p = gsl_ran_gaussian(RNG, 1.);}} 
  else                 {std::cerr << "RNG not initialized." << std::endl; return;}
}
// ....................................................................... //
void Grid::makeRandomQ (const double spinVariation) {
  if (RNG_initialized) {
    for (auto & q : this->Q) {
      q  = randBetween(-spinVariation, +spinVariation);
      q  = std::fmod(q, 2 * PI);
      q += (q < 0  ?  2 * PI  :  0);
    }
  } else {
    std::cerr << "Failed to (re)initialize grid: RNG not ready" << std::endl;
  }
}
// ----------------------------------------------------------------------- //
const std::vector<int> & Grid::getNeighboursOf(const int i) const {return this->neighbours[i];}
