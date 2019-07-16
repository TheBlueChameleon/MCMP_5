/* TODO: Module description
 * 
 */

#ifndef GRID_H
#define GRID_H

// ========================================================================= //
// dependencies

#include <vector>

// ========================================================================= //
// symbols and proc

enum class Startcondition {
  COLDSTART, HOTSTART
};

inline double StartconditionToQ (Startcondition SC);

// ========================================================================= //
// class

class Grid {
private:
  int L = -1;
  int V = -1;
  
  std::vector<double>           q;
  std::vector<double>           p;
  std::vector<std::vector<int>> neighbours;
  
  void init(const int L, const double spinVariation);
  
public:
  Grid(const int L);
  Grid(const int L, const double spinVariation);
  
  void makeRandomP ();
  void updateQwithP();
  
  std::vector<double> & getSpinField    () const;
  double                getSpinMagnitude() const;
  
  void setSpinChanges(const std::vector<double> &P);
};

#endif//GRID_H
