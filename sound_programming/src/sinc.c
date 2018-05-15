/* code copied (and possibly modified) from http://floor13.sakura.ne.jp/book06/book06.html */

#include <math.h>

double sinc(double x)
{
  double y;
  
  if (x == 0.0)
  {
    y = 1.0;
  }
  else
  {
    y = sin(x) / x;
  }
  
  return y;
}
