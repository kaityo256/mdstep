#pragma once
#include <vector>
//------------------------------------------------------------------------
struct Atom {
  double qx, qy, qz;
  double px, py, pz;
};
//------------------------------------------------------------------------
struct Pair {
  int i, j;
};
//------------------------------------------------------------------------
class Variables {
public:
  std::vector<Atom> atoms;
  std::vector<int> neighbor_list;
  std::vector<int> i_position;
  std::vector<int> j_count;
  double time;
  double zeta; //For Nose-Hoover
  Variables(void) {time = 0.0;zeta = 0.0;}
  void add_atoms(double x, double y, double v);
  void export_cdview(void);
  int number_of_atoms(void) {return static_cast<int>(atoms.size());}
  void set_initial_velocity(const double);
  void make_neighbor_list(std::vector<Pair> &pair);
};
//------------------------------------------------------------------------
