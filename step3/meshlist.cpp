//------------------------------------------------------------------------
#include <iostream>
#include <assert.h>
#include <algorithm>
#include "systemparam.hpp"
#include "meshlist.hpp"
//------------------------------------------------------------------------
MeshList::MeshList(void){
  const double SL = CUTOFF+MARGIN;
  m = static_cast<int>(L/SL);
  mesh_size = static_cast<double>(L)/m;
  number_of_mesh = m*m*m;
  count.resize(number_of_mesh);
  index.resize(number_of_mesh);
}
//------------------------------------------------------------------------
void
MeshList::make_mesh(Variables *vars){
  Atom *atoms = vars->atoms.data();
  const int pn = vars->number_of_atoms();
  std::vector<int> particle_position(pn);
  std::vector<int> pointer(pn);
  std::fill(particle_position.begin(),particle_position.end(),0);
  std::fill(count.begin(),count.end(),0);
  std::fill(pointer.begin(),pointer.end(),0);

  double im = 1.0/mesh_size;
  for(int i=0;i<pn;i++){
    int ix = static_cast<int>(atoms[i].qx * im);
    int iy = static_cast<int>(atoms[i].qy * im);
    int iz = static_cast<int>(atoms[i].qz * im);
    int index = ix + iy*m + iz *m*m;
    assert(index >=0);
    assert(index < number_of_mesh);
    count[index]++;
    particle_position[i] = index;
  }
  index[0] = 0;
  int sum = 0;
  for(int i=0;i<number_of_mesh-1;i++){
    sum += count[i];
    index[i+1] = sum;
  }
  for(int i=0;i<pn;i++){
    int pos = particle_position[i];
    int j = index[pos] + pointer[pos];
    sorted_buffer[j] = i;
    ++pointer[pos];
  }
}
//------------------------------------------------------------------------
