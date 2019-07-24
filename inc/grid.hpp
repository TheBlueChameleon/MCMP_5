/* grid.hpp
 * defines a grid of size LxL of spins with one continuous degree of freedom
 * (component q[i] in [0, 2pi]) as well as a prospective change of spin
 * (component p[i] in [0, 2pi]).
 */

#ifndef GRID_H
#define GRID_H

// ========================================================================= //
// dependencies

#include <vector>
#include <array>

#include "globals.hpp"

// ========================================================================= //
// symbols and proc

enum class Startcondition {
  COLDSTART, HOTSTART
};

inline double StartconditionToVariation (Startcondition SC) {
  switch (SC) {
    case Startcondition::COLDSTART : return 0.0;
    case Startcondition::HOTSTART  : return PI;
  }
  
  // this should not be able to happen
  std::cerr << "Warning: Encountered an invalid Startcondition symbol" << std::endl;
  return NAN;
}

// ========================================================================= //
// class

class Grid {
private:
  int L = -1;
  int V = -1;
  
  std::vector<double>           Q;
  std::vector<double>           P;
  std::vector<std::vector<int>> neighbours;
  
  void init(const int L, const double spinVariation);
  
public:
  Grid(const int L);
  Grid(const int L, const Startcondition SC);
  Grid(const int L, const double spinVariation);
  
  // ....................................................................... //
  // lattice geometry
  
  int getL() const;
  int getV() const;
  
  int arrayCoordinateX (const int idx) const;
  int arrayCoordinateY (const int idx) const;
  
  int arrayIndex(const int x, const int y) const;
  
  std::array<int, 4> neighboursOf (const int x, const int y) const;
  
  // ....................................................................... //
  // vector getters / setters
  
  const std::vector<double> & getQ() const;
  const std::vector<double> & getP() const;
  
  void setQ       (const std::vector<double> & Q);
  void setP       (const std::vector<double> & P);
  
  void swapQ      (      std::vector<double> & Q);
  void swapP      (      std::vector<double> & P);
  
  void makeRandomP ();
  void makeRandomQ (const double spinVariation);
  
  const std::vector<int> & getNeighboursOf(const int i) const;
};

#endif//GRID_H
