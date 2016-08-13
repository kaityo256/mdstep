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
    std::vector<int> indexes;
    std::vector<int> sorted_buffer;
    void search(int index, Variables *vars, std::vector<Pair> &pairs);
    void search_other(int id, int ix, int iy, int iz, Variables *vars, std::vector<Pair> &pairs);
  public:
    MeshList(void);
    void make_pair(Variables *vars, std::vector<Pair> &pairs);
    void set_number_of_atoms(int pn){
      sorted_buffer.resize(pn);
    }
};
//------------------------------------------------------------------------
