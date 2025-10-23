////////////////////////////////////////
//// general functions               //
////                                //
//////////////////////////////////////
#include "../include/general.hpp"

double xlinearInter(
  const double x1,
  const double y1,
  const double x2,
  const double y2,
  const double y
)
{
  double x = 0.0;

  x = x1 + (y - y1)*(x2 - x1)/(y2 - y1);

  return x;
}
