/* code copied (and possibly modified) from http://floor13.sakura.ne.jp/book06/book06.html */

#include <math.h>

void IIR_LPF(double fc, double Q, double a[], double b[])
{
  fc = tan(M_PI * fc) / (2.0 * M_PI);
  
  a[0] = 1.0 + 2.0 * M_PI * fc / Q + 4.0 * M_PI * M_PI * fc * fc;
  a[1] = (8.0 * M_PI * M_PI * fc * fc - 2.0) / a[0];
  a[2] = (1.0 - 2.0 * M_PI * fc / Q + 4.0 * M_PI * M_PI * fc * fc) / a[0];
  b[0] = 4.0 * M_PI * M_PI * fc * fc / a[0];
  b[1] = 8.0 * M_PI * M_PI * fc * fc / a[0];
  b[2] = 4.0 * M_PI * M_PI * fc * fc / a[0];
  
  a[0] = 1.0;
}

void IIR_HPF(double fc, double Q, double a[], double b[])
{
  fc = tan(M_PI * fc) / (2.0 * M_PI);
  
  a[0] = 1.0 + 2.0 * M_PI * fc / Q + 4.0 * M_PI * M_PI * fc * fc;
  a[1] = (8.0 * M_PI * M_PI * fc * fc - 2.0) / a[0];
  a[2] = (1.0 - 2.0 * M_PI * fc / Q + 4.0 * M_PI * M_PI * fc * fc) / a[0];
  b[0] = 1.0 / a[0];
  b[1] = -2.0 / a[0];
  b[2] = 1.0 / a[0];
  
  a[0] = 1.0;
}

void IIR_BPF(double fc1, double fc2, double a[], double b[])
{
  fc1 = tan(M_PI * fc1) / (2.0 * M_PI);
  fc2 = tan(M_PI * fc2) / (2.0 * M_PI);
  
  a[0] = 1.0 + 2.0 * M_PI * (fc2 - fc1) + 4.0 * M_PI * M_PI * fc1 * fc2;
  a[1] = (8.0 * M_PI * M_PI * fc1 * fc2 - 2.0) / a[0];
  a[2] = (1.0 - 2.0 * M_PI * (fc2 - fc1) + 4.0 * M_PI * M_PI * fc1 * fc2) / a[0];
  b[0] = 2.0 * M_PI * (fc2 - fc1) / a[0];
  b[1] = 0.0;
  b[2] = -2.0 * M_PI * (fc2 - fc1) / a[0];
  
  a[0] = 1.0;
}

void IIR_BEF(double fc1, double fc2, double a[], double b[])
{
  fc1 = tan(M_PI * fc1) / (2.0 * M_PI);
  fc2 = tan(M_PI * fc2) / (2.0 * M_PI);
  
  a[0] = 1.0 + 2.0 * M_PI * (fc2 - fc1) + 4.0 * M_PI * M_PI * fc1 * fc2;
  a[1] = (8.0 * M_PI * M_PI * fc1 * fc2 - 2.0) / a[0];
  a[2] = (1.0 - 2.0 * M_PI * (fc2 - fc1) + 4.0 * M_PI * M_PI * fc1 * fc2) / a[0];
  b[0] = (4.0 * M_PI * M_PI * fc1 * fc2 + 1.0) / a[0];
  b[1] = (8.0 * M_PI * M_PI * fc1 * fc2 - 2.0) / a[0];
  b[2] = (4.0 * M_PI * M_PI * fc1 * fc2 + 1.0) / a[0];
  
  a[0] = 1.0;
}

void IIR_resonator(double fc, double Q, double a[], double b[])
{
  fc = tan(M_PI * fc) / (2.0 * M_PI);
  
  a[0] = 1.0 + 2.0 * M_PI * fc / Q + 4.0 * M_PI * M_PI * fc * fc;
  a[1] = (8.0 * M_PI * M_PI * fc * fc - 2.0) / a[0];
  a[2] = (1.0 - 2.0 * M_PI * fc / Q + 4.0 * M_PI * M_PI * fc * fc) / a[0];
  b[0] = 2.0 * M_PI * fc / Q / a[0];
  b[1] = 0.0;
  b[2] = -2.0 * M_PI * fc / Q / a[0];
  
  a[0] = 1.0;
}

void IIR_notch(double fc, double Q, double a[], double b[])
{
  fc = tan(M_PI * fc) / (2.0 * M_PI);
  
  a[0] = 1.0 + 2.0 * M_PI * fc / Q + 4.0 * M_PI * M_PI * fc * fc;
  a[1] = (8.0 * M_PI * M_PI * fc * fc - 2.0) / a[0];
  a[2] = (1.0 - 2.0 * M_PI * fc / Q + 4 * M_PI * M_PI * fc * fc) / a[0];
  b[0] = (4.0 * M_PI * M_PI * fc * fc + 1.0) / a[0];
  b[1] = (8.0 * M_PI * M_PI * fc * fc - 2.0) / a[0];
  b[2] = (4.0 * M_PI * M_PI * fc * fc + 1.0) / a[0];
  
  a[0] = 1.0;
}

void IIR_low_shelving(double fc, double Q, double g, double a[], double b[])
{
  fc = tan(M_PI * fc) / (2.0 * M_PI);
  
  a[0] = 1.0 + 2.0 * M_PI * fc / Q + 4.0 * M_PI * M_PI * fc * fc;
  a[1] = (8.0 * M_PI * M_PI * fc * fc - 2.0) / a[0];
  a[2] = (1.0 - 2.0 * M_PI * fc / Q + 4.0 * M_PI * M_PI * fc * fc) / a[0];
  b[0] = (1.0 + sqrt(1.0 + g) * 2.0 * M_PI * fc / Q + 4.0 * M_PI * M_PI * fc * fc * (1.0 + g)) / a[0];
  b[1] = (8.0 * M_PI * M_PI * fc * fc * (1.0 + g) - 2.0) / a[0];
  b[2] = (1.0 - sqrt(1.0 + g) * 2.0 * M_PI * fc / Q + 4.0 * M_PI * M_PI * fc * fc * (1.0 + g)) / a[0];
  
  a[0] = 1.0;
}

void IIR_high_shelving(double fc, double Q, double g, double a[], double b[])
{
  fc = tan(M_PI * fc) / (2.0 * M_PI);
  
  a[0] = 1.0 + 2.0 * M_PI * fc / Q + 4.0 * M_PI * M_PI * fc * fc;
  a[1] = (8.0 * M_PI * M_PI * fc * fc - 2.0) / a[0];
  a[2] = (1.0 - 2.0 * M_PI * fc / Q + 4.0 * M_PI * M_PI * fc * fc) / a[0];
  b[0] = ((1.0 + g) + sqrt(1.0 + g) * 2.0 * M_PI * fc / Q + 4.0 * M_PI * M_PI * fc * fc) / a[0];
  b[1] = (8.0 * M_PI * M_PI * fc * fc - 2.0 * (1.0 + g)) / a[0];
  b[2] = ((1.0 + g) - sqrt(1.0 + g) * 2.0 * M_PI * fc / Q + 4.0 * M_PI * M_PI * fc * fc) / a[0];
  
  a[0] = 1.0;
}

void IIR_peaking(double fc, double Q, double g, double a[], double b[])
{
  fc = tan(M_PI * fc) / (2.0 * M_PI);
  
  a[0] = 1.0 + 2.0 * M_PI * fc / Q + 4.0 * M_PI * M_PI * fc * fc;
  a[1] = (8.0 * M_PI * M_PI * fc * fc - 2.0) / a[0];
  a[2] = (1.0 - 2.0 * M_PI * fc / Q + 4.0 * M_PI * M_PI * fc * fc) / a[0];
  b[0] = (1.0 + 2.0 * M_PI * fc / Q * (1.0 + g) + 4.0 * M_PI * M_PI * fc * fc) / a[0];
  b[1] = (8.0 * M_PI * M_PI * fc * fc - 2.0) / a[0];
  b[2] = (1.0 - 2.0 * M_PI * fc / Q * (1.0 + g) + 4.0 * M_PI * M_PI * fc * fc) / a[0];
  
  a[0] = 1.0;
}

void IIR_filtering(const double x[], double y[], int L, const double a[], const double b[], int I, int J)
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
    for (m = 1; m <= I; m++)
    {
      if (n - m >= 0)
      {
        y[n] += -a[m] * y[n - m];
      }
    }
  }
}
