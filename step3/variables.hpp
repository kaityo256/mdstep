#pragma once
#include <vector>
//------------------------------------------------------------------------
struct Atom {
  double qx, qy, qz;
  double px, py, pz;
};
//------------------------------------------------------------------------
class Variables {
public:
  std::vector<Atom> atoms;
  double time;
  Variables(void) {time = 0.0;}
  void add_atoms(double x, double y, double v);
  void export_cdview(void);
  int number_of_atoms(void) {return static_cast<int>(atoms.size());}
  void set_initial_velocity(const double);
};
//------------------------------------------------------------------------
