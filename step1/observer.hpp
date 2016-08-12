#pragma once
#include "variables.hpp"
//------------------------------------------------------------------------
class Observer{
  public:
    double kinetic_energy(Variables *vars); 
    double potential_energy(Variables *vars); 
    double temperature(Variables *vars){return kinetic_energy(vars)/1.5;}
    double total_energy(Variables *vars){return kinetic_energy(vars) + potential_energy(vars);}
};
//------------------------------------------------------------------------
