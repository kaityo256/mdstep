#pragma once
#include "variables.hpp"
#include "observer.hpp"
#include "meshlist.hpp"
//------------------------------------------------------------------------
class MD {
private:
  Variables *vars;
  Observer *obs;
  MeshList *mesh;
  std::vector<Pair> pairs;
  double margin_length;
  void makeconf(void);
  void update_position(void);
  void calculate_force(void);
  void calculate_force_pair(void);
  void calculate_force_list(void);
  void periodic(void);
  void calculate(void);
  void make_pair(void);
  void check_pairlist(void);
  void velocity_scaling(const double aimed_temperature);

public:
  MD(void);
  ~MD(void);
  void run(void);
};
//------------------------------------------------------------------------
