#include"oblateness.h"
Oblateness::Oblateness(double Req, double Rpole, double alpha, double sma, double inc, double u1, double u2):Req_(Req),Rpole_(Rpole),alpha_(alpha),sma_(sma),inc_(inc),u1_(u1),u2_(u2){
  AM=1.0/IM;
  NDIV=(1+(IM-1)/NTAB);
  EPS=1.2e-7;
  RNMX= (1.0 - EPS);
}
Oblateness::~Oblateness(){
}

double Oblateness::Ran1_(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--){
      k = (*idum)/IQ;
      *idum = IA*(*idum - k*IQ)-IR*k;
      if (*idum<0) *idum += IM;
      if (j<NTAB) iv[j] = *idum;
    }
    iy = iv[0];
  }

  k= (*idum)/IQ;
  *idum = IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0 ) *idum += IM;
  j = iy/NDIV;
  iy = iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
static long seed=0;

/* limb darkening function */
double Oblateness::LimbDarkening_(double r2)
{
	double t = 1-r2;
	return((u1_+2*u2_)*sqrt(t)-u2_*t);
}

/* find the intersection point on the line (x1,y1) to (x2,y2) */
void Oblateness::FindRoot_(double x1, double y1, double x2, double y2, double xc, double yc, double *xval, double *yval, int *flag)
{
	double d1, d2, d;
	d1 = (x1-xc)*(x1-xc)+(y1-yc)*(y1-yc)-1;
	if(fabs(d1)<1e-12) { *xval = x1; *yval = y1; *flag = 1; return;}
	d2 = (x2-xc)*(x2-xc)+(y2-yc)*(y2-yc)-1;
	if(fabs(d2)<1e-12) { *flag = 0; return;} /* to avoid double counting */
	/* if d1*d2>0, no intersection point is there */
	if(d1*d2>0) { *flag = 0; return;}

	double x, y;
	while(fabs(x1-x2)>1e-6 || fabs(y1-y2)>1e-6)
	{
		x = (x1+x2)/2.0; y = (y1+y2)/2.0;
		d = (x-xc)*(x-xc)+(y-yc)*(y-yc)-1;
		if(fabs(d)<1e-12) { 
			*xval = x; *yval = y; 
			*flag = 1; 
			return;
		}
		d1 = (x1-xc)*(x1-xc)+(y1-yc)*(y1-yc)-1.0;
		if(d*d1>0) {x1 = x; y1 = y;}
		else {x2 =x; y2 = y;}
	}
	d1 = (x1-xc)*(x1-xc)+(y1-yc)*(y1-yc)-1.0;
	if(d1>0) {
		*xval = x2; *yval = y2; *flag = 1;
	}
	else {
		*xval = x1; *yval = y1; *flag = 1;
	}
}

/* find the intersection points on the ellipse */
void Oblateness::IntersectionPoints_(double *x, double *y, int *n, double a, double b, double xc, double yc)
{
	int N=500, i, count=0;
	double phi1, phi2, dphi=2*pi/N;
	double xval, yval;
	double x1, y1, x2, y2;
	int flag;
	for(i=0; i<N; i++)
	{
		phi1 = i*dphi; phi2 = (i+1)*dphi;
		x1 = a*cos(phi1); x2 = a*cos(phi2);
		y1 = b*sin(phi1); y2 = b*sin(phi2);
		FindRoot_(x1, y1, x2, y2, xc, yc, &xval, &yval, &flag);
		if(flag == 0) continue;
		if(count > 1) {printf("More intersection points than expected!\n"); exit(1);}
		x[count] = xval; y[count] = yval; count++;
		/* there is one intersection point between (x1,y1) and (x2,y2) */
	}
	*n = count;
	return;
}

/* area calculation */
double Oblateness::TriangleArea_(double *x, double *y)
{
	double dx = x[1]-x[0];
	double dy = y[1]-y[0];
	double d1 = sqrt(dx*dx+dy*dy);
	double k, h;
	if(fabs(dx)<1e-6) h = fabs(x[1]);
	else {
		k = dy/dx;
		h = fabs(y[0]-k*x[0])/sqrt(1+k*k);
	}
	double area = 0.5*d1*h;
	return(area);
}

double Oblateness::CircleChordArea_(double *x, double *y, double xc, double yc, double rc)
{
	double x1=x[0]-xc, y1=y[0]-yc;
	double x2=x[1]-xc, y2=y[1]-yc;
	double dbeta;
	double l12;
	l12 = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
	double l1, l2;
	l1 = rc*rc; l2 = l1;
	dbeta = acos((l1+l2-l12)/(2*sqrt(l1*l2)));
	double area;
	area = 0.5*rc*rc*(dbeta-sin(dbeta));
	return(area);
}

/* calculate the deficite flux of the star at phase phi */
void Oblateness::relativeFlux(double *phi, int np, double *deficitFlux, int nf)
//void relativeFlux(double *variables, int n, double phi, double *deficitFlux, double *circleAnalogy)
{
	/* variables */
/*  variables[0] = Req;
	variables[1] = Rpol;
	varibales[2] = alpha;
	variables[3] = sma;
	variables[4] = period;
	variables[5] = inclination;
	variables[6] = u1; limbdarkening coefficient
	variables[7] = u2; limbdarkening coefficient
	variables[8] = n; observation times
	variables[9] = percentage; phase interval
*/
  double *d = new double [np];
  double *index = new double[np];
  double clc1=0.0,clc2=0.0,clc3=0.0;
  clock_t init,final;
  double b0=sma_*cos(inc_);
  double epsilon = Req_/Rpole_;
  double epsilon2 = epsilon*epsilon;
  double *xc = new double[np], *yc=new double[np], *thetaa = new double[np],*thetab = new double[np];
  double beta = 0.0;
  double eta1, eta2, theta, R, x, y;
  double contribution=0.0, separation1, separation2;
  int Nrays = 3000; //hard wire in sample parameters
  double *etaa=new double[Nrays], *etab=new double[Nrays];
  //printf("%x %x %x %x %x %x\n",d,index,xc,yc,etaa,etab);
  for (int i=0;i<np;i++){
    /* distance between centers of the planet and the star */
    init = clock();
	  if(cos(phi[i]*2*pi) <0){
      //when planet is at the back
      deficitFlux[i] = 0.0;
      index[i] = 0.0;
      continue;
    }
    d[i]=sma_*sqrt(sin(phi[i]*2*pi)*sin(phi[i]*2*pi)+cos(phi[i]*2*pi)*cos(phi[i]*2*pi)*cos(inc_)*cos(inc_));
	
	  /* the embedded sphere's contribution */
	  //*circleAnalogy = standardCurve(variables[1], d, variables[6], variables[7]);

	  /* if no transit happens */
	  if(d[i] >= (1.0+Req_)){ 
      deficitFlux[i] = 0.0;
      index[i] = 0.0;
      continue;
    }
	
	  /* angle between the major axis and the line connecting two centers */
	  if(phi[i]<=0) beta = alpha_+asin(b0/d[i]); /* ingress */
	  if(phi[i]>0)  beta = alpha_-pi-asin(b0/d[i]); /* egress */

	    /* star's position and intersection points */
	  xc[i] = d[i]*cos(beta);
	  yc[i] = d[i]*sin(beta);
	  double xinter[2], yinter[2];
	  int npoints;
    final = clock()-init;
    clc1+=((double)final/((double)CLOCKS_PER_SEC));
    init = clock();
	  IntersectionPoints_(xinter, yinter, &npoints, Req_, Rpole_, xc[i], yc[i]);
    final = clock()-init;
    clc2+=((double)final/((double)CLOCKS_PER_SEC));
	
	  if(d[i]>1.0 && npoints<2) { /* the planet is outside the stellar plane, and there is no intersectiosn */
		  deficitFlux[i] = 0.0;
      index[i] = 0.0;
      continue;
	  }

	  /* calculate the limitation region */
	  double theta1=0, theta2=0, t;
	  double F00=0; /* F(x;a,b,alpha,0,0)-F(x;b,0,0) */
	  if(d[i]<1.0 && npoints<2) { /* the planet is inside the stellar plane and there is no intersections */
		  theta1 = 0.0;
	  	theta2 = 2*pi;
		  F00 = pi*Rpole_*(Req_-Rpole_);
	    index[i] = 1.0;
    } else {
      index[i] = 0.0;
    }
	  if(npoints==2) { /* two intersections */
		  if(yinter[0] >= 0)
	  		theta1 = acos(xinter[0]/Req_);
		  else
	  		theta1 = 2*pi-acos(xinter[0]/Req_);
		  if(yinter[1] >= 0)
	  		theta2 = acos(xinter[1]/Req_);
		  else
	  		theta2 = 2*pi-acos(xinter[1]/Req_);
		  if(theta1 > theta2){
			  t = theta1;
			  theta1 = theta2;
			  theta2 = t;
		  }
	  }

	  double theta0 = (theta1+theta2)*0.5;
	  double xtest, ytest;
	  xtest = Req_*cos(theta0);
	  ytest = Rpole_*sin(theta0);
	  double dtest;
	  dtest = (xtest-xc[i])*(xtest-xc[i])+(ytest-yc[i])*(ytest-yc[i])-1.0;
	  if(dtest>0) {
		  t = theta2-2*pi;
		  theta2 = theta1;
		  theta1 = t;
	  }
    thetaa[i] = theta1;
    thetab[i] = theta2;
	  double chordArea1, chordArea2;
	  double dtheta=theta2-theta1;
	  double ksi, psi, cc, p2;
	  if(npoints == 2) {
		  if(dtheta<=pi){
		    chordArea1 = dtheta*Req_*Rpole_*0.5-TriangleArea_(xinter, yinter);
      }else{
		    chordArea1 = dtheta*Req_*Rpole_*0.5+TriangleArea_(xinter, yinter);
      }
		  chordArea2 = CircleChordArea_(xinter, yinter, xc[i], yc[i], 1.0);
		  if(d[i]>=(Rpole_+1.0)) /* the two circles are separated */
      {
		    F00 = chordArea1+chordArea2;
      } else if(d[i]>(1.0-Rpole_)) { /* the two circles have intersection points */
			  ksi = acos((Rpole_*Rpole_+d[i]*d[i]-1.0)/(2*Rpole_*d[i]));
			  psi = acos((1.0+d[i]*d[i]-Rpole_*Rpole_)/(2.0*d[i]));
			  p2 = Rpole_*Rpole_;
			  cc = ksi*p2+psi-sqrt(4*d[i]*d[i]-(1+d[i]*d[i]-p2)*(1+d[i]*d[i]-p2))*0.5;
			  F00 = chordArea1+chordArea2-cc;
		  } else /* the small circle is included in the larger one */{
		    F00 = chordArea1+chordArea2-pi*Rpole_*Rpole_;
      }
	  }

    /*Deal with the ingress and egress*/
    if (index[i]==0.0){    
	    contribution = 0.0;
            init = clock();
	    for(int j=0; j<Nrays; j++)
	    { 
		    eta1 = Ran1_(&seed);
		    eta2 = Ran1_(&seed);
		    theta = (1.0-eta2)*theta1+eta2*theta2;
		    R = sqrt(eta1+(1-eta1)/epsilon2);
		    x = Req_*R*cos(theta);
		    y = Rpole_*R*sin(theta);
		    /* distance to the center of the star */
		    separation1 = (x-xc[i])*(x-xc[i])+(y-yc[i])*(y-yc[i]);
        
		    /* if the ray locates outside the stellar plane, or inside the circular plane of planet, ignore it */
		    if(separation1 > 1.0) continue;
		    /* distance to the center of the planet */
		    separation2 = x*x+y*y;
		    if(separation2 < Rpole_*Rpole_) continue;
		    /* otherwise, add its contribution on */
		    contribution += LimbDarkening_(separation1);
	    }
	
      final = clock()-init;
      clc3+=((double)final/((double)CLOCKS_PER_SEC));
	    double totalArea, Istar;
	    totalArea = dtheta*(Req_*Req_-Rpole_*Rpole_)*0.5/epsilon;
	    Istar = contribution*totalArea/Nrays;
	    deficitFlux[i] = (1.0-u1_-u2_)*F00+Istar;
      //printf("%d %f %f %f %f\n", i, F00,contribution,xc[i],yc[i]);
    } else {
      deficitFlux[i] = (1.0-u1_-u2_)*F00;
    }
  }
  /********************************************************/
  /*Deal with planet completely inside the star differently*/ 
  int flag = 0;
  init = clock();
  //double contribution=0.0, separation1, separation2;
  for (int i=0;i<np;i++){
    if(index[i]==1.0 and flag==0){
      contribution = 0.0;
      for(int j=0; j<Nrays;j++)
	    {
        //printf("%d %f %f %f %f %f\n",i,phi[i],d0,correction,xc[i],yc[i]);
		    eta1 = Ran1_(&seed);
		    eta2 = Ran1_(&seed);
		    theta = (1.0-eta2)*thetaa[i]+eta2*thetab[i];
		    R = sqrt(eta1+(1-eta1)/epsilon2);
		    x = Req_*R*cos(theta);
		    y = Rpole_*R*sin(theta);
		    etaa[j]=-1; etab[j]=-1; 
        /* distance to the center of the star */
		    separation1 = (x-xc[i])*(x-xc[i])+(y-yc[i])*(y-yc[i]);
		    /* if the ray locates outside the stellar plane, or inside the circular plane of planet, ignore it */
		    if(separation1 > 1.0) continue;
		    /* distance to the center of the planet */
		    separation2 = x*x+y*y;
		    if(separation2 < Rpole_*Rpole_) continue;
		    /* otherwise, add its contribution on */
		    contribution += LimbDarkening_(separation1);
        etaa[j]=eta1; etab[j]=eta2;
        flag = 1;
    }
      double totalArea, Istar;
	    totalArea = (thetab[i]-thetaa[i])*(Req_*Req_-Rpole_*Rpole_)*0.5/epsilon;
	    Istar = contribution*totalArea/Nrays;
	    deficitFlux[i] += Istar;

    } else { 
      if (index[i]==1.0 and flag ==1){	
        contribution = 0.0;
        for (int j=0;j<Nrays;j++){
          if(etaa[j]==-1) continue;
          
		      theta = (1.0-etab[j])*thetaa[i]+etab[j]*thetab[i];
		      R = sqrt(etaa[j]+(1-etaa[j])/epsilon2);
		      x = Req_*R*cos(theta);
		      y = Rpole_*R*sin(theta);
		      separation1 = (x-xc[i])*(x-xc[i])+(y-yc[i])*(y-yc[i]);
		      contribution += LimbDarkening_(separation1);
        }
      double totalArea, Istar;
	    totalArea = (thetab[i]-thetaa[i])*(Req_*Req_-Rpole_*Rpole_)*0.5/epsilon;
	    Istar = contribution*totalArea/Nrays;
	    deficitFlux[i] += Istar;

      }
    }
  }
  final = clock()-init;
  clc3+=((double)final/((double)CLOCKS_PER_SEC));

  //printf("three zones: %f %f %f\n",clc1,clc2,clc3);
  delete [] d;
  delete [] index;
  delete [] xc;
  delete [] yc;
  delete [] thetaa; 
  delete [] thetab;
  delete [] etaa;
  delete [] etab;
}


