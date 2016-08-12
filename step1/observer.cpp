//------------------------------------------------------------------------
#include "systemparam.hpp"
#include "observer.hpp"
//------------------------------------------------------------------------
double
Observer::kinetic_energy(Variables *vars){
  double k = 0;
  for(auto &a: vars->atoms){
    k += a.px * a.px;
    k += a.py * a.py;
    k += a.pz * a.pz;
  }
  k/= static_cast<double>(vars->number_of_atoms());
  return k*0.5;
};
//------------------------------------------------------------------------
double
Observer::potential_energy(Variables *vars){
  double v = 0.0;
  const int pn = vars->number_of_atoms();
  Atom *atoms = &(vars->atoms[0]);
  for(int i=0;i<pn-1;i++){
    for(int j=i+1;j<pn;j++){
      double dx = atoms[j].qx - atoms[i].qx;
      double dy = atoms[j].qy - atoms[i].qy;
      double dz = atoms[j].qz - atoms[i].qz;
      adjust_periodic(dx,dy,dz);
      double r2 = (dx*dx+dy*dy+dz*dz);
      if(r2>CL2)continue;
      double r6 = r2*r2*r2;
      double r12 = r6*r6;
      v += 4.0*(1.0/r12 - 1.0/r6) + C0;
    }
  }
  v /= static_cast<double>(pn);
  return v;
}
//------------------------------------------------------------------------
