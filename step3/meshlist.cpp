//------------------------------------------------------------------------
#include <iostream>
#include "systemparam.hpp"
#include "meshlist.hpp"
//------------------------------------------------------------------------
MeshList::MeshList(void){
  const double SL = CUTOFF+MARGIN;
  int m = static_cast<int>(L/SL);
  mesh_size = static_cast<double>(L)/m;
  std::cout << mesh_size << " " << m << std::endl;
  std::cout << m*mesh_size << " " << L << std::endl;
  const int m3 = m*m*m;
  count.resize(m3);
  index.resize(m3);
}
//------------------------------------------------------------------------
