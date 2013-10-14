#ifndef ELLIPTIC_H_
#define ELLIPTIC_H_
#include <math.h>
#include <iostream>
class Elliptic{
  public:
    Elliptic();
    ~Elliptic();
    void Calee(double *k, int lk, double *ee, int lee);
    void Calkk(double *k, int lk, double *kk, int lkk);
    void Ellpic_bulirsch(double *n, int ln, double *k, int lk, double *nk, int lnk);
  private:
     double Pi_,ea1_,ea2_,ea3_,ea4_,eb1_,eb2_,eb3_,eb4_;
     double ka1_,ka2_,ka3_,ka4_,kb1_,kb2_,kb3_,kb4_,ka0_,kb0_;

     };
#endif
