//------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <random>
#include "systemparam.hpp"
#include "variables.hpp"
//------------------------------------------------------------------------
void
Variables::add_atoms(double x, double y, double z) {
  Atom a;
  a.qx = x;
  a.qy = y;
  a.qz = z;
  a.px = 0.0;
  a.py = 0.0;
  a.pz = 0.0;
  atoms.push_back(a);
}
//------------------------------------------------------------------------
void
Variables::set_initial_velocity(const double V0) {
  std::mt19937 mt(2);
  std::uniform_real_distribution<double> ud(0.0, 1.0);
  double avx = 0.0;
  double avy = 0.0;
  double avz = 0.0;
  for (auto &a : atoms) {
    double z = ud(mt) * 2.0 - 1.0;
    double phi = 2.0 * ud(mt) * M_PI;
    double vx = V0 * sqrt(1 - z * z) * cos(phi);
    double vy = V0 * sqrt(1 - z * z) * sin(phi);
    double vz = V0 * z;
    a.px = vx;
    a.py = vy;
    a.pz = vz;
    avx += vx;
    avy += vy;
    avz += vz;
  }
  const int pn = atoms.size();
  avx /= static_cast<double>(pn);
  avy /= static_cast<double>(pn);
  avz /= static_cast<double>(pn);
  for (auto &a : atoms) {
    a.px -= avx;
    a.py -= avy;
    a.pz -= avz;
  }
}
//------------------------------------------------------------------------
void
Variables::export_cdview(void) {
  static int count = 0;
  char filename[256];
  sprintf(filename, "conf%03d.cdv", count);
  ++count;
  std::ofstream ofs(filename);
  int i = 0;
  for (auto &a : atoms) {
    ofs << i << " " << "0" << " ";
    ofs << a.qx << " ";
    ofs << a.qy << " ";
    ofs << a.qz << " ";
    ofs << std::endl;
    ++i;
  }
}
//------------------------------------------------------------------------
