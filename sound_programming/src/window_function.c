/* code copied (and possibly modified) from http://floor13.sakura.ne.jp/book06/book06.html */

#include <math.h>

void Hanning_window(double w[], int N)
{
  int n;
  
  if (N % 2 == 0) /* Nが偶数のとき */
  {
    for (n = 0; n < N; n++)
    {
      w[n] = 0.5 - 0.5 * cos(2.0 * M_PI * n / N);
    }
  }
  else /* Nが奇数のとき */
  {
    for (n = 0; n < N; n++)
    {
      w[n] = 0.5 - 0.5 * cos(2.0 * M_PI * (n + 0.5) / N);
    }
  }
}
