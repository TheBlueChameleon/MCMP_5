/* grid.cpp
 * defines a grid of size LxL of spins with one continuous degree of freedom
 * (component q[i] in [0, 2pi]) as well as a prospective change of spin
 * (component p[i] in [0, 2pi]).
 */

// ========================================================================= //
// dependencies
 
#include <vector>
#include <cmath>
#include <iostream>

#include "globals.hpp"
#include "grid.hpp"

// ========================================================================= //
// enum translator

inline double StartconditionToQ (Startcondition SC) {
  switch (SC) {
    case Startcondition::COLDSTART : return 0.0;
    case Startcondition::HOTSTART  : return 2 * PI;
  }
  
  // this should not be able to happen
  std::cerr << "Warning: Encountered an invalid IsingStart symbol" << std::endl;
  return NAN;
}
