////////////////////////////////////////
//// general functions               //
////                                //
//////////////////////////////////////
#include "../include/general.hpp"

float xlinearInter(
  const float x1,
  const float y1,
  const float x2,
  const float y2,
  const float y
)
{
  float x = 0.0;

  x = x1 + (y - y1)*(x2 - x1)/(y2 - y1);

  return x;
}
