#pragma once
#include "variables.hpp"
#include "observer.hpp"

class MD{
  private:
    Variables *vars;
    Observer *obs;
    void makeconf(void);
    void update_position(void);
    void calculate_force(void);
    void periodic(void);
    void calculate(void);

  public:
    MD(void);
    ~MD(void);
    void run(void);
};
