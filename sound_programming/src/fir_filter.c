/* code copied (and possibly modified) from http://floor13.sakura.ne.jp/book06/book06.html */

#include <math.h>

static double sinc(double x)
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


void FIR_LPF(double fe, int J, double b[], double w[])
{
  int m;
  int offset;
  
  offset = J / 2;
  for (m = -J / 2; m <= J / 2; m++)
  {
    b[offset + m] = 2.0 * fe * sinc(2.0 * M_PI * fe * m);
  }
  
  for (m = 0; m < J + 1; m++)
  {
    b[m] *= w[m];
  }
}

void FIR_HPF(double fe, int J, double b[], double w[])
{
  int m;
  int offset;
  
  offset = J / 2;
  for (m = -J / 2; m <= J / 2; m++)
  {
    b[offset + m] = sinc(M_PI * m) - 2.0 * fe * sinc(2.0 * M_PI * fe * m);
  }
  
  for (m = 0; m < J + 1; m++)
  {
    b[m] *= w[m];
  }
}

void FIR_BPF(double fe1, double fe2, int J, double b[], double w[])
{
  int m;
  int offset;
  
  offset = J / 2;
  for (m = -J / 2; m <= J / 2; m++)
  {
    b[offset + m] = 2.0 * fe2 * sinc(2.0 * M_PI * fe2 * m)
                  - 2.0 * fe1 * sinc(2.0 * M_PI * fe1 * m);
  }
  
  for (m = 0; m < J + 1; m++)
  {
    b[m] *= w[m];
  }
}

void FIR_BEF(double fe1, double fe2, int J, double b[], double w[])
{
  int m;
  int offset;
  
  offset = J / 2;
  for (m = -J / 2; m <= J / 2; m++)
  {
    b[offset + m] = sinc(M_PI * m)
                  - 2.0 * fe2 * sinc(2.0 * M_PI * fe2 * m)
                  + 2.0 * fe1 * sinc(2.0 * M_PI * fe1 * m);
  }
  
  for (m = 0; m < J + 1; m++)
  {
    b[m] *= w[m];
  }
}

void FIR_filtering(const double x[], double y[], int L, double b[], int J)
{
  int n, m;
  
  for (n = 0; n < L; n++)
  {
    for (m = 0; m <= J; m++)
    {
      if (n - m >= 0)
      {
        y[n] += b[m] * x[n - m];
      }
    }
  }
}
