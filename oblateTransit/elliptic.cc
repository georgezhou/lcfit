#include"elliptic.h"
Elliptic::Elliptic(){
  Pi_=3.1415926535898;
  ea1_=0.44325141463;
  ea2_=0.06260601220;
  ea3_=0.04757383546;
  ea4_=0.01736506451;
  eb1_=0.24998368310;
  eb2_=0.09200180037;
  eb3_=0.04069697526;
  eb4_=0.00526449639;
 
  ka0_=1.38629436112;
  ka1_=0.09666344259;
  ka2_=0.03590092383;
  ka3_=0.03742563713;
  ka4_=0.01451196212;
  kb0_=0.5;
  kb1_=0.12498593597;
  kb2_=0.06880248576;
  kb3_=0.03328355346;
  kb4_=0.00441787012;

}
Elliptic::~Elliptic(){
}
void Elliptic::Calee(double *k, int lk, double *ee, int lee){
  double m1,logm1,ee1,ee2;
  for (int i=0; i<lk; i++){
    m1=1.-pow(k[i],2);
    logm1 = log(m1);
    ee1=1.+m1*(ea1_+m1*(ea2_+m1*(ea3_+m1*ea4_)));
    ee2=m1*(eb1_+m1*(eb2_+m1*(eb3_+m1*eb4_)))*(-logm1);
    ee[i] = ee1+ee2;
  }

}
void Elliptic::Calkk(double *k, int lk, double *kk, int lkk){
  double m1,logm1,ek1,ek2;
  for (int i=0; i<lk; i++){
    m1=1.-pow(k[i],2);
    logm1 = log(m1);
    ek1=ka0_+m1*(ka1_+m1*(ka2_+m1*(ka3_+m1*ka4_)));
    ek2=(kb0_+m1*(kb1_+m1*(kb2_+m1*(kb3_+m1*kb4_))))*logm1;
    kk[i] = ek1-ek2;
  }
}

void Elliptic::Ellpic_bulirsch(double *n, int ln, double *k, int lk, double *nk, int lnk){
  double *p = new double [ln];
  double *d = new double [ln];
  double *kc = new double [ln];
  double *e = new double [ln];
  double *m0 = new double [ln];
  double *c = new double [ln];
  double *f = new double [ln];
  double *g = new double [ln];
  double err;
  for (int i = 0; i<ln; i++){
    kc[i]=sqrt(1.-pow(k[i],2)); 
    if(n[i]+1<0) {
      std::cout << 'Negative p' << std::endl;
      exit(1);
    } else {
      p[i] = sqrt(n[i]+1);
      d[i] = 1./p[i];
    }
    e[i] = kc[i];
    m0[i] = 1.;
    c[i] = 1.;
    f[i] = 0.;
    g[i] = 0;
  }
  do{
    err = 0;
    for (int i=0; i<ln; i++){
      f[i] = c[i];
      c[i] += d[i]/p[i];
      g[i] = e[i]/p[i];
      d[i] = 2.*(f[i]*g[i]+d[i]);
      p[i] += g[i];
      g[i] = m0[i];
      m0[i] += kc[i];
      if(fabs(1.-kc[i]/g[i])>err) err = 1.-kc[i]/g[i];
    }
    if (err > 1.e-8){
      for (int i=0; i<ln; i++){
        kc[i] = 2*sqrt(e[i]); 
        e[i]=kc[i]*m0[i];
      }
    }
    else{
      for (int i=0; i<ln; i++){
        nk[i] = 0.5*Pi_*(c[i]*m0[i]+d[i])/(m0[i]*(m0[i]+p[i]));
      }
      break;
    }
  } while (1);
  delete [] p;
  delete [] d;
  delete [] kc;
  delete [] e;
  delete [] m0;
  delete [] c;
  delete [] f;
  delete [] g;
}
