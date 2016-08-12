#pragma once
#include <vector>
//------------------------------------------------------------------------
struct Atom{
  double qx,qy,qz;
  double px,py,pz;
  Atom(void){
  }
  Atom(double x, double y, double z){
    qx = x;
    qy = y;
    qz = z;
    px = 0.0;
    py = 0.0;
    pz = 0.0;
  }
};
//------------------------------------------------------------------------
class Variables{
  public:
  std::vector<Atom> atoms;
  double time;
  Variables(void){time = 0.0;}
  void add_atoms(double x, double y, double v);
  void export_cdview(void);
  int number_of_atoms(void){return static_cast<int>(atoms.size());}
  void set_initial_velocity(const double);
};
//------------------------------------------------------------------------
