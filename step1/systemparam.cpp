//------------------------------------------------------------------------
#include "systemparam.hpp"
//------------------------------------------------------------------------
void
adjust_periodic(double &dx, double &dy, double &dz) {
  const double LH = L * 0.5;
  if (dx < -LH)dx += L;
  if (dx > LH) dx -= L;
  if (dy < -LH)dy += L;
  if (dy > LH) dy -= L;
  if (dz < -LH)dz += L;
  if (dz > LH) dz -= L;
}
//------------------------------------------------------------------------
