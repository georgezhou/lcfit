#include "oblateness.h"
int main()
{
  double variables[10];
	double rmean;
	double f;
	double alpha, sma, period, inc, u1, u2;
	int N, i, n=10;
	double percent;
	FILE *fp;
	if((fp=fopen("system_config", "r")) == NULL)
	{
		printf("Input file not found!\n");
		exit(1);
	}
	for(i=0; i<10; i++)
	{
		fscanf(fp, "%lf", variables+i);
		while(fgetc(fp) != '\n');
	}
	variables[2] = variables[2]/180.0*pi;
	variables[5] = variables[5]/180.0*pi;
	/* parameters */
	rmean = variables[0];
	f	  = variables[1];
	alpha = variables[2]; //angle between planet orbit and its major axis
	sma   = variables[3]; //semimajor axis normalized by radius
	period= variables[4];
	inc   = variables[5];
	u1    = variables[6];
	u2    = variables[7];
	N     = (int)variables[8];
	percent= variables[9]; //phase

	double req = variables[0]/sqrt(1-variables[1]); //Req
	double rpol = sqrt(1-variables[1])*variables[0]; //Rpole
	Oblateness obl = Oblateness(req,rpol,alpha,sma,inc,u1,u1);

	double b0 = sma*cos(variables[5]); //inpact parameter
	if(b0>(1+req))
	{
		printf("No trnasit happens!\n");
		exit(1);
	}
	double dphi = 2*percent/(variables[8]-1);
	double* phi = new double [N];
  double* deficitFlux = new double [N];
	double totalFlux,  amp, circleAnalogy;
	FILE *lc;

	if(fabs(f)<1e-3)
	  lc = fopen("limbDarkenedLC-mean.dat", "w");
	else
	  lc = fopen("limbDarkenedLC.dat", "w");
	totalFlux = pi*(1.0-variables[6]/3-variables[7]/6);
  for(int i =0;i<N;i++){
    phi[i]=-1*percent+dphi*i;
  }
		/* functions to calculate the relative flux at this moment */
	obl.relativeFlux(phi,N, deficitFlux, N);
		//amp = circleAnalogy-deficitFlux/totalFlux;
  for(int i =0;i<N;i++){
    amp = deficitFlux[i]/totalFlux;
		fprintf(lc, "%f %f\n", phi[i], amp);
	}
	fclose(fp);
	fclose(lc);
	return(1);
}
