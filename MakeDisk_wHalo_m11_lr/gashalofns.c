#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"




#define  N  800000   /* number of mass bins for halo */



/* a number of tables */


static double *hg_mt,*hg_m2r,*hg_r2m;
static double *hg_rhofinal,*hg_rho2r;
static double *hg_rt;




  /* ---------------------------------
       halogas profile normalization 
     --------------------------------- */
/*
  RHO_0= M_HALOGAS/(4*PI*qromb(halogas_profile_int,1e-6,Rvir));
  printf("\nhalogas normalization= %g\n",RHO_0);
*/




/* ----------------------------------------------------
    Halogas mass, rho, distribution function vs. r
   ---------------------------------------------------- */
double gashalo_mass(double r)
{
  double x;
  if(M_GASHALO <= 0) { return 0; }

  if(r>hg_rt[N])
    x=hg_mt[N];
  else
    splint(hg_rt,hg_mt,hg_m2r,N,r,&x);
  return x;
}


double gashalo_q_to_r(double q)
{
  double m,x;
  m=hg_mt[N]*q;

  splint(hg_mt,hg_rt,hg_r2m,N,m,&x);
  //printf("hg_mt[N]=%g , q=%g, m=%g , x=%g , N=%d , \n",hg_mt[N],q,m,x,N); fflush(stdout);
  return x;
}



double gashalo_rho(double r)
{
  double x;

  if(M_GASHALO <= 0) { return 0; }

  if(r>hg_rt[N])
    x=0;
  else
    splint(hg_rt,hg_rhofinal,hg_rho2r,N,r,&x);

  return x;
}




double halogas_profile_int(double r)
{
  /* unnormailized (so need 4 PI RHO_0 when called, but notice divide by RS in the
	  appropriate place here not before calling this function)
          and used for volume integration so have r^2 in there too */
  /*
  Beta Model - what do I use for rc? Makino rc=0.22rs  or use RS as above? */
  double RC;
  RC= GasHalo_Rc_over_Rs*RS;
  return r*r*pow(1+(r/RC)*(r/RC), -1.5*GasHalo_Beta);


  /* Hydrostatic eq. - Isothermal */
/*
  double b= 12.0;
  return r*r*pow(1+r/RS,b*RS/r);
*/

  /*
  Hydrostatic eq. - T function of r */
/*
  double alpha= 1.0;
  double x=r/RS;
  return r*r*exp(-(3.0/alpha)*x*(1-log(1+x)/x)/(log(1+x)-x/(1+x)));
*/
}


double mass_cumulative_gashalo(double R)
{
  return gashalo_mass(R);
}



void setup_halogas_massprofile(void)
{
  int i;
/*double q,s,qq,f,f_,ds,r;*/
  double dlogr;


  hg_mt=dvector(1,N);
  hg_r2m=dvector(1,N);
  hg_m2r=dvector(1,N);
  hg_rhofinal=dvector(1,N);
  hg_rho2r=dvector(1,N);

  hg_rt=dvector(1,N);
  dlogr= (log10(Rvir)-log10(1.0e-7))/N;
  dlogr= (log10(10.0*Rvir)-log10(1.0e-6))/N;
  for(i=1;i<=N;i++)
    {
	hg_rt[i]= pow(10, log10(1.0e-6) + i*dlogr);
    }

  for(i=2, hg_mt[1]=0 ;i<=N;i++)
    {
      hg_rhofinal[i]=RHO_0*halogas_profile_int(hg_rt[i])/(hg_rt[i]*hg_rt[i]);
      hg_mt[i]=4.*PI*RHO_0*qromb(halogas_profile_int,1.0e-6,hg_rt[i]);
      //printf("hg_rt=%g  hg_mt=%g RHO_0=%g hgprof=%g \n",hg_rt[i],hg_mt[i],RHO_0,halogas_profile_int(hg_rt[i]));fflush(stdout);
    }

  spline(hg_rt,hg_mt,N,1e40,1e40,hg_m2r);
  spline(hg_mt,hg_rt,N,1e40,1e40,hg_r2m);
  spline(hg_rt,hg_rhofinal,N,1e40,1e40,hg_rho2r);
}





