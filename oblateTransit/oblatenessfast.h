#ifndef OBLATENESS_H
#define OBLATENESS_H
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<gsl/gsl_sf_ellint.h>
#define pi 3.1415926535897932384626433832795028841971693993751
using namespace std;
class Oblateness{
  public:
    Oblateness(double Req, double Rpole, double alpha, double sma, double inc, double u1, double u2);
    ~Oblateness();
    void relativeFlux(double *phi, int np, double *deficitFlux, int nf);
  private:
    double Ran1_(long *idum);
    double LimbDarkening_(double r2);
    void FindRoot_(double x1, double y1, double x2, double y2, double xc, double yc, double *xval, double *yval, int *flag);
    void IntersectionPoints_(double *x, double *y, int *n, double a, double b, double xc, double yc);
    double TriangleArea_(double *x, double *y);
    double CircleChordArea_(double *x, double *y, double xc, double yc, double rc);
    double u1_,u2_,inc_,Req_,Rpole_,alpha_,sma_;
    static const int IA=16807;
    static const int IM=2147483647;
    double AM;
    double NDIV;
    double EPS;
    double RNMX;
    static const int IQ=127773;
    static const int IR=2836;
    static const int NTAB=32;

};
#endif
