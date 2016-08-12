//------------------------------------------------------------------------
#include <iostream>
#include <assert.h>
#include <math.h>
#include "md.hpp"
#include "systemparam.hpp"
#include "observer.hpp"
#include "variables.hpp"
//------------------------------------------------------------------------
MD::MD(void) {
  vars = new Variables();
  obs = new Observer();
  margin_length = 0.0;
}
//------------------------------------------------------------------------
MD::~MD(void) {
  delete vars;
  delete obs;
}
//------------------------------------------------------------------------
void
MD::makeconf(void) {
  const double density = 0.50;
  const double s = 1.0 / pow(density * 0.25, 1.0 / 3.0);
  const double hs = s * 0.5;
  const int is = static_cast<int>(L / s);
  for (int iz = 0; iz < is; iz++) {
    for (int iy = 0; iy < is; iy++) {
      for (int ix = 0; ix < is; ix++) {
        vars->add_atoms(ix * s, iy * s, iz * s);
        vars->add_atoms(ix * s + hs, iy * s, iz * s);
        vars->add_atoms(ix * s, iy * s + hs, iz * s);
        vars->add_atoms(ix * s, iy * s, iz * s + hs);
      }
    }
  }
  vars->set_initial_velocity(1.0);
}
//------------------------------------------------------------------------
void
MD::make_pair(void){
  pairs.clear();
  const int pn = vars->number_of_atoms();
  Atom *atoms = vars->atoms.data();
  for (int i = 0; i < pn - 1; i++) {
    for (int j = i + 1; j < pn; j++) {
      double dx = atoms[j].qx - atoms[i].qx;
      double dy = atoms[j].qy - atoms[i].qy;
      double dz = atoms[j].qz - atoms[i].qz;
      adjust_periodic(dx, dy, dz);
      double r2 = (dx * dx + dy * dy + dz * dz);
      if (r2 > ML2)continue;
      Pair p;
      p.i = i;
      p.j = j;
      pairs.push_back(p);
    }
  }
}
//------------------------------------------------------------------------
void
MD::check_pairlist(void){
  double vmax2 = 0.0;
  for (auto &a : vars->atoms) {
    double v2 = a.px*a.px + a.py*a.py + a.pz*a.pz;
    if (vmax2 < v2) vmax2 = v2;
  }
  double vmax = sqrt(vmax2);
  margin_length -= vmax*2.0*dt;
  if(margin_length < 0.0){
    margin_length = MARGIN;
    make_pair();
  }
}
//------------------------------------------------------------------------
void
MD::update_position(void) {
  const double dt2 = dt * 0.5;
  for (auto &a : vars->atoms) {
    a.qx += a.px * dt2;
    a.qy += a.py * dt2;
    a.qz += a.pz * dt2;
  }
}
//------------------------------------------------------------------------
void
MD::calculate_force(void) {
  const int pn = vars->number_of_atoms();
  Atom *atoms = vars->atoms.data();
  for (int i = 0; i < pn - 1; i++) {
    for (int j = i + 1; j < pn; j++) {
      double dx = atoms[j].qx - atoms[i].qx;
      double dy = atoms[j].qy - atoms[i].qy;
      double dz = atoms[j].qz - atoms[i].qz;
      adjust_periodic(dx, dy, dz);
      double r2 = (dx * dx + dy * dy + dz * dz);
      if (r2 > CL2)continue;
      double r6 = r2 * r2 * r2;
      double df = (24.0 * r6 - 48.0) / (r6 * r6 * r2) * dt;
      atoms[i].px += df * dx;
      atoms[i].py += df * dy;
      atoms[i].pz += df * dz;
      atoms[j].px -= df * dx;
      atoms[j].py -= df * dy;
      atoms[j].pz -= df * dz;
    }
  }
}
//------------------------------------------------------------------------
void
MD::calculate_force_pair(void) {
  const int pp = pairs.size();
  Atom *atoms = vars->atoms.data();
  for(int k=0;k<pp;k++){
    const int i = pairs[k].i;
    const int j = pairs[k].j;
    double dx = atoms[j].qx - atoms[i].qx;
    double dy = atoms[j].qy - atoms[i].qy;
    double dz = atoms[j].qz - atoms[i].qz;
    adjust_periodic(dx, dy, dz);
    double r2 = (dx * dx + dy * dy + dz * dz);
    if (r2 > CL2)continue;
    double r6 = r2 * r2 * r2;
    double df = (24.0 * r6 - 48.0) / (r6 * r6 * r2) * dt;
    atoms[i].px += df * dx;
    atoms[i].py += df * dy;
    atoms[i].pz += df * dz;
    atoms[j].px -= df * dx;
    atoms[j].py -= df * dy;
    atoms[j].pz -= df * dz;
  }
}
//------------------------------------------------------------------------
void
MD::periodic(void) {
  for (auto &a : vars->atoms) {
    if (a.qx < 0.0) a.qx += L;
    if (a.qy < 0.0) a.qy += L;
    if (a.qz < 0.0) a.qz += L;
    if (a.qx > L) a.qx -= L;
    if (a.qy > L) a.qy -= L;
    if (a.qz > L) a.qz -= L;
    assert(a.qx < L);
    assert(a.qy < L);
    assert(a.qz < L);
  }
}
//------------------------------------------------------------------------
void
MD::calculate(void) {
  update_position();
  check_pairlist();
  //make_pair();
  calculate_force_pair();
  update_position();
  periodic();
  vars->time += dt;
}
//------------------------------------------------------------------------
void
MD::run(void) {
  makeconf();
  const int STEPS = 10000;
  const int OBSERVE = 100;
  for (int i = 0; i < STEPS; i++) {
    if ( (i % OBSERVE) == 0) {
      double k = obs->kinetic_energy(vars);
      double v = obs->potential_energy(vars);
      std::cout << vars->time << " ";
      std::cout << k << " ";
      std::cout << v << " ";
      std::cout << k + v << std::endl;
      vars->export_cdview();
    }
    calculate();
  }
}
//------------------------------------------------------------------------
