#pragma once
#include "variables.hpp"
//------------------------------------------------------------------------
class Observer {
public:
  double kinetic_energy(Variables *vars);
  double potential_energy(Variables *vars, std::vector<Pair> & pairs);
  double temperature(Variables *vars) {return kinetic_energy(vars) / 1.5;}
  double total_energy(Variables *vars, std::vector<Pair> &pairs) {return kinetic_energy(vars) + potential_energy(vars, pairs);}
};
//------------------------------------------------------------------------
