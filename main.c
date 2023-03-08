
#include "cmath.h"

#include <stdio.h>
#include <math.h>


int function(int n, double t, double x[], double dxdt[])
{
  dxdt[0] = -71 * x[0] - 70 * x[1] + pow(2.72, 1 - t * t);
  dxdt[1] = x[0] + sin(1 - t);
  return 0;
}

void mine_runge(double * x, double t, double h) {
  double k1[2];
  function(2, t, x, k1);
  k1[0] *= h;
  k1[1] *= h;

  double k2[2];
  double x_tmp[2] = {x[0] + k1[0] / 3, x[1] + k1[0] / 3};
  function(2, t + h / 3, x_tmp, k2);
  k2[0] *= h;
  k2[1] *= h;

  double k3[2];
  x_tmp[0]=x[0] + k1[0] * 0.16 + 0.24 * k2[0];
  x_tmp[1]=x[1] + k1[1] * 0.16 + 0.24 * k2[1];
  function(2, t + 0.4 * h, x_tmp, k3);
  k3[0] *= h;
  k3[1] *= h;

  double k4[2];
  x_tmp[0]=x[0] + k1[0] * 0.25 - 3 * k2[0] + 3.75 * k3[0];
  x_tmp[1]=x[1] + k1[1] * 0.25 - 3 * k2[1] + 3.75 * k3[1];
  function(2, t + h, x_tmp, k4);
  k4[0] *= h;
  k4[1] *= h;

  double k5[2];
  x_tmp[0]=x[0] + (k1[0] * 6 + 90 * k2[0] - 50 * k3[0] + 8 * k4[0]) / 81;
  x_tmp[1]=x[1] + (k1[1] * 6 + 90 * k2[1] - 50 * k3[1] + 8 * k4[1]) / 81;
  function(2, t + 2 * h / 3, x_tmp, k5);
  k5[0] *= h;
  k5[1] *= h;

  double k6[2];
  x_tmp[0]=x[0] + (k1[0] * 6 + 36 * k2[0] + 10 * k3[0] + 8 * k4[0]) / 75;
  x_tmp[1]=x[1] + (k1[1] * 6 + 36 * k2[1] + 10 * k3[1] + 8 * k4[1]) / 75;
  function(2, t + 4 * h / 5, x_tmp, k6);
  k6[0] *= h;
  k6[1] *= h;

  x[0] = x[0] + (23 * k1[0] + 125 * k3[0] - 81 * k5[0] + 125 * k6[0]);
  x[1] = x[1] + (23 * k1[1] + 125 * k3[1] - 81 * k5[1] + 125 * k6[1]);
}




int main()
{
  int n = 2, fail = 0;
  double relerr = 0.0001;
  double abserr = 0.0001;
  double t = 0.0;
  double h;
  double h_print = 0.2;
  int flag = 1, nfe, maxfe = 1000000;
  rkfinit(n, &fail);
  printf("%s\n\n", cmathmsg(RKFINIT_C, fail));

  double x[2] = {0.0, 1.0};
  double dxdt[2];
  printf("standart rkf45 \n");
  if (fail == 0) {
    for (int step = 0; step <= 20; ++step) {
      double tout = step * h_print;
      t = tout - h_print;
      rkf45(function, n, x, dxdt, &t, tout, &relerr, abserr, &h, &nfe, maxfe, &flag);
      if (flag != 2) {
        printf("%s\n", cmathmsg(RKF45_C, flag));
        break;
      }
      printf("%lf %lf %lf \n", tout, x[0], x[1]);
    }
    rkfend();
    printf("-----------------------------------------\n");
  }

  printf("mine_runge\n");
  x[0] = 0.0;
  x[1] = 1.0;
  h = 0.000001;
  double tout = 0.0;
  for (int step = 0; step <= 4000000; ++step) {
    mine_runge(x, tout, h);
    if (step % 100000 == 0) {
      printf("%lf %lf %lf \n", tout, x[0], x[1]);
    }
    tout += h;
  }
}



