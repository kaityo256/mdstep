#pragma once
//------------------------------------------------------------------------
#include <vector>
#include "variables.hpp"
//------------------------------------------------------------------------
class MeshList{
  private:
    double mesh_size;
    int m;
    int number_of_mesh;
    std::vector<int> count;
    std::vector<int> index;
  public:
    MeshList(void);
    void make_mesh(Variables *vars);
};
//------------------------------------------------------------------------
