//------------------------------------------------------------------------
#include <iostream>
#include <assert.h>
#include <algorithm>
#include "systemparam.hpp"
#include "meshlist.hpp"
//------------------------------------------------------------------------
MeshList::MeshList(void) {
  const double SL = CUTOFF + MARGIN;
  m = static_cast<int>(L / SL) - 1;
  assert(m > 2);
  mesh_size = static_cast<double>(L) / m;
  assert(mesh_size > SL);
  number_of_mesh = m * m * m;
  count.resize(number_of_mesh);
  indexes.resize(number_of_mesh);
}
//------------------------------------------------------------------------
void
MeshList::make_pair(Variables *vars, std::vector<Pair> &pairs) {
  pairs.clear();
  Atom *atoms = vars->atoms.data();
  const int pn = vars->number_of_atoms();
  std::vector<int> particle_position(pn);
  std::vector<int> pointer(pn);
  std::fill(particle_position.begin(), particle_position.end(), 0);
  std::fill(count.begin(), count.end(), 0);
  std::fill(pointer.begin(), pointer.end(), 0);

  double im = 1.0 / mesh_size;
  for (int i = 0; i < pn; i++) {
    int ix = static_cast<int>(atoms[i].qx * im);
    int iy = static_cast<int>(atoms[i].qy * im);
    int iz = static_cast<int>(atoms[i].qz * im);
    if (ix < 0) ix += m;
    if (ix >= m) ix -= m;
    if (iy < 0) iy += m;
    if (iy >= m) iy -= m;
    if (iz < 0) iz += m;
    if (iz >= m) iz -= m;

    int index = ix + iy * m + iz * m * m;
    assert(index >= 0);
    assert(index < number_of_mesh);
    count[index]++;
    particle_position[i] = index;
  }
  indexes[0] = 0;
  int sum = 0;
  for (int i = 0; i < number_of_mesh - 1; i++) {
    sum += count[i];
    indexes[i + 1] = sum;
  }
  for (int i = 0; i < pn; i++) {
    int pos = particle_position[i];
    int j = indexes[pos] + pointer[pos];
    sorted_buffer[j] = i;
    ++pointer[pos];
  }
  for (int i = 0; i < number_of_mesh; i++) {
    search(i, vars, pairs);
  }
  vars->make_neighbor_list(pairs);
}
//------------------------------------------------------------------------
void
MeshList::search_other(int id, int ix, int iy, int iz, Variables *vars, std::vector<Pair> &pairs) {
  if (ix < 0) ix += m;
  if (ix >= m ) ix -= m;
  if (iy < 0) iy += m;
  if (iy >= m ) iy -= m;
  if (iz < 0) iz += m;
  if (iz >= m ) iz -= m;
  int id2 = ix + iy * m + iz * m * m;
  Atom *atoms = vars->atoms.data();
  for (int k = indexes[id]; k < indexes[id] + count[id]; k++) {
    for (int m = indexes[id2]; m < indexes[id2] + count[id2]; m++) {
      int i = sorted_buffer[k];
      int j = sorted_buffer[m];
      double dx = atoms[j].qx - atoms[i].qx;
      double dy = atoms[j].qy - atoms[i].qy;
      double dz = atoms[j].qz - atoms[i].qz;
      adjust_periodic(dx, dy, dz);
      double r2 = (dx * dx + dy * dy + dz * dz);
      if (r2 > ML2)continue;
      Pair p;
      if (i < j) {
        p.i = i;
        p.j = j;
      } else {
        p.i = j;
        p.j = i;
      }
      pairs.push_back(p);
    }
  }
}
//------------------------------------------------------------------------
void
MeshList::search(int id, Variables *vars, std::vector<Pair> &pairs) {
  int ix = id % m;
  int iy = (id / m) % m;
  int iz = (id / m / m);
  search_other(id, ix + 1, iy, iz, vars, pairs);
  search_other(id, ix - 1, iy + 1, iz, vars, pairs);
  search_other(id, ix    , iy + 1, iz, vars, pairs);
  search_other(id, ix + 1, iy + 1, iz, vars, pairs);

  search_other(id, ix - 1, iy    , iz + 1, vars, pairs);
  search_other(id, ix    , iy    , iz + 1, vars, pairs);
  search_other(id, ix + 1, iy    , iz + 1, vars, pairs);

  search_other(id, ix - 1, iy - 1, iz + 1, vars, pairs);
  search_other(id, ix    , iy - 1, iz + 1, vars, pairs);
  search_other(id, ix + 1, iy - 1, iz + 1, vars, pairs);

  search_other(id, ix - 1, iy + 1, iz + 1, vars, pairs);
  search_other(id, ix    , iy + 1, iz + 1, vars, pairs);
  search_other(id, ix + 1, iy + 1, iz + 1, vars, pairs);
  // Registration of self box
  int si = indexes[id];
  int n = count[id];
  Atom *atoms = vars->atoms.data();
  for (int k = si; k < si + n - 1; k++) {
    for (int m = k + 1; m < si + n; m++) {
      int i = sorted_buffer[k];
      int j = sorted_buffer[m];
      double dx = atoms[j].qx - atoms[i].qx;
      double dy = atoms[j].qy - atoms[i].qy;
      double dz = atoms[j].qz - atoms[i].qz;
      adjust_periodic(dx, dy, dz);
      double r2 = (dx * dx + dy * dy + dz * dz);
      if (r2 > ML2)continue;
      Pair p;
      if (i < j) {
        p.i = i;
        p.j = j;
      } else {
        p.i = j;
        p.j = i;
      }
      pairs.push_back(p);
    }
  }
}
//------------------------------------------------------------------------
