#pragma once
//------------------------------------------------------------------------
#include <vector>
//------------------------------------------------------------------------
class MeshList{
  private:
    double mesh_size;
    std::vector<int> count;
    std::vector<int> index;
  public:
    MeshList(void);
};
//------------------------------------------------------------------------
