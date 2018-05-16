/* code copied (and possibly modified) from http://floor13.sakura.ne.jp/book06/book06.html */

#include <math.h>

#define PI M_PI

void IIR_HPF(double fc, double Q, double a[], double b[])
{
  fc = tan(PI * fc) / (2.0 * PI);
  
  a[0] = 1.0 + 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc;
  a[1] = (8.0 * PI * PI * fc * fc - 2.0) / a[0];
  a[2] = (1.0 - 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc) / a[0];
  b[0] = 1.0 / a[0];
  b[1] = -2.0 / a[0];
  b[2] = 1.0 / a[0];
  
  a[0] = 1.0;
}

void IIR_BPF(double fc1, double fc2, double a[], double b[])
{
  fc1 = tan(PI * fc1) / (2.0 * PI);
  fc2 = tan(PI * fc2) / (2.0 * PI);
  
  a[0] = 1.0 + 2.0 * PI * (fc2 - fc1) + 4.0 * PI * PI * fc1 * fc2;
  a[1] = (8.0 * PI * PI * fc1 * fc2 - 2.0) / a[0];
  a[2] = (1.0 - 2.0 * PI * (fc2 - fc1) + 4.0 * PI * PI * fc1 * fc2) / a[0];
  b[0] = 2.0 * PI * (fc2 - fc1) / a[0];
  b[1] = 0.0;
  b[2] = -2.0 * PI * (fc2 - fc1) / a[0];
  
  a[0] = 1.0;
}

void IIR_BEF(double fc1, double fc2, double a[], double b[])
{
  fc1 = tan(PI * fc1) / (2.0 * PI);
  fc2 = tan(PI * fc2) / (2.0 * PI);
  
  a[0] = 1.0 + 2.0 * PI * (fc2 - fc1) + 4.0 * PI * PI * fc1 * fc2;
  a[1] = (8.0 * PI * PI * fc1 * fc2 - 2.0) / a[0];
  a[2] = (1.0 - 2.0 * PI * (fc2 - fc1) + 4.0 * PI * PI * fc1 * fc2) / a[0];
  b[0] = (4.0 * PI * PI * fc1 * fc2 + 1.0) / a[0];
  b[1] = (8.0 * PI * PI * fc1 * fc2 - 2.0) / a[0];
  b[2] = (4.0 * PI * PI * fc1 * fc2 + 1.0) / a[0];
  
  a[0] = 1.0;
}



void IIR_notch(double fc, double Q, double a[], double b[])
{
  fc = tan(PI * fc) / (2.0 * PI);
  
  a[0] = 1.0 + 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc;
  a[1] = (8.0 * PI * PI * fc * fc - 2.0) / a[0];
  a[2] = (1.0 - 2.0 * PI * fc / Q + 4 * PI * PI * fc * fc) / a[0];
  b[0] = (4.0 * PI * PI * fc * fc + 1.0) / a[0];
  b[1] = (8.0 * PI * PI * fc * fc - 2.0) / a[0];
  b[2] = (4.0 * PI * PI * fc * fc + 1.0) / a[0];
  
  a[0] = 1.0;
}

void IIR_low_shelving(double fc, double Q, double g, double a[], double b[])
{
  fc = tan(PI * fc) / (2.0 * PI);
  
  a[0] = 1.0 + 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc;
  a[1] = (8.0 * PI * PI * fc * fc - 2.0) / a[0];
  a[2] = (1.0 - 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc) / a[0];
  b[0] = (1.0 + sqrt(1.0 + g) * 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc * (1.0 + g)) / a[0];
  b[1] = (8.0 * PI * PI * fc * fc * (1.0 + g) - 2.0) / a[0];
  b[2] = (1.0 - sqrt(1.0 + g) * 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc * (1.0 + g)) / a[0];
  
  a[0] = 1.0;
}

void IIR_high_shelving(double fc, double Q, double g, double a[], double b[])
{
  fc = tan(PI * fc) / (2.0 * PI);
  
  a[0] = 1.0 + 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc;
  a[1] = (8.0 * PI * PI * fc * fc - 2.0) / a[0];
  a[2] = (1.0 - 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc) / a[0];
  b[0] = ((1.0 + g) + sqrt(1.0 + g) * 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc) / a[0];
  b[1] = (8.0 * PI * PI * fc * fc - 2.0 * (1.0 + g)) / a[0];
  b[2] = ((1.0 + g) - sqrt(1.0 + g) * 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc) / a[0];
  
  a[0] = 1.0;
}

void IIR_peaking(double fc, double Q, double g, double a[], double b[])
{
  fc = tan(PI * fc) / (2.0 * PI);
  
  a[0] = 1.0 + 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc;
  a[1] = (8.0 * PI * PI * fc * fc - 2.0) / a[0];
  a[2] = (1.0 - 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc) / a[0];
  b[0] = (1.0 + 2.0 * PI * fc / Q * (1.0 + g) + 4.0 * PI * PI * fc * fc) / a[0];
  b[1] = (8.0 * PI * PI * fc * fc - 2.0) / a[0];
  b[2] = (1.0 - 2.0 * PI * fc / Q * (1.0 + g) + 4.0 * PI * PI * fc * fc) / a[0];
  
  a[0] = 1.0;
}


