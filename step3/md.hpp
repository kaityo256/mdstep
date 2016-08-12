#pragma once
#include "variables.hpp"
#include "observer.hpp"
//------------------------------------------------------------------------
struct Pair{
  int i,j;
};
//------------------------------------------------------------------------
class MD {
private:
  Variables *vars;
  Observer *obs;
  std::vector<Pair> pairs;
  double margin_length;
  void makeconf(void);
  void update_position(void);
  void calculate_force(void);
  void calculate_force_pair(void);
  void periodic(void);
  void calculate(void);
  void make_pair(void);
  void check_pairlist(void);

public:
  MD(void);
  ~MD(void);
  void run(void);
};
//------------------------------------------------------------------------
