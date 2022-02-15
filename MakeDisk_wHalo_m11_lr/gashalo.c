#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"





void compute_velocity_dispersions_gashalo(void)
{
  int i,j;
  double z,R;
  double rho;
 

  printf("gashalo velocity dispersion field...\n"); fflush(stdout);  

  for(i=0;i<=RSIZE;i++)
    {
      printf("gashalo A, %d\n",i);
      
      for(j=0;j<=ZSIZE;j++)
	{
	  xl[j+1]=list_z[j];
	  yl[j+1]=Dphi_z[i][j] * comp_rho_gashalo(list_R[i],list_z[j]);
	}

      spline(xl,yl,ZSIZE+1,1e40,1e40,D2yl);

      /* Xiangcheng: This is equivalent to solving hydrostatic equilibrium equations,
	 as it actually integrates rho*dPhi/dz, which gives pressure. 
	 For boundary conditions, it is problematic to assume p=0. 
	 Below is motivated by Virial theorm, which reduces scatter. */
      VelDispRz_gashalo[i][ZSIZE] = (1./3.) *
		fabs(list_R[i]*Dphi_R[i][ZSIZE]+list_z[ZSIZE]*Dphi_z[i][ZSIZE]);
      VelDispRz_gashalo[i][ZSIZE] *= comp_rho_gashalo(list_R[i],list_z[ZSIZE]);

      for(j=ZSIZE-1;j>=0;j--)
	{

	  VelDispRz_gashalo[i][j] =  VelDispRz_gashalo[i][j+1];
	  if(fabs(yl[j+2])>1e-100 && fabs(yl[j+1])>1e-100)
	    VelDispRz_gashalo[i][j]+=
	      qromb(splint_xl_yl_D2yl,list_z[j],list_z[j+1]);
	}
    }

  for(i=0;i<=RSIZE;i++)
    {
      printf("gashalo B, %d\n",i);

      for(j=0;j<=ZSIZE;j++)
	{
	  xl[j+1]=list_z[j];
	  yl[j+1]=Dphi_z_dR[i][j] * comp_rho_gashalo(list_RplusdR[i],list_z[j]);

	}

      spline(xl,yl,ZSIZE+1,1e40,1e40,D2yl);
 
     /* Xiangcheng: the same as above, we even don't have to bother it is R+dR */
      VelDispRz_gashalo[i][ZSIZE] = (1./3.) *
		fabs(list_R[i]*Dphi_R[i][ZSIZE]+list_z[ZSIZE]*Dphi_z[i][ZSIZE]);
      VelDispRz_gashalo[i][ZSIZE] *= comp_rho_gashalo(list_R[i],list_z[ZSIZE]);

      for(j=ZSIZE - 1;j>=0;j--)
	{
	  VelDispRz_dR_gashalo[i][j] =  VelDispRz_dR_gashalo[i][j+1];
	  if(fabs(yl[j+2])>1e-100 && fabs(yl[j+1])>1e-100)
	    VelDispRz_dR_gashalo[i][j]+=
	      qromb(splint_xl_yl_D2yl,list_z[j],list_z[j+1]);
	}
      
    }

  
  for(i=0;i<=RSIZE;i++)
    {
      for(j=0;j<=ZSIZE;j++)
	{
	  R=list_R[i];
	  z=list_z[j];
	  
	  rho = comp_rho_gashalo(R,z);

	  if(rho>0)
	    {
	      if(i>0)
		VelDispPhi_gashalo[i][j]=R/rho * (VelDispRz_dR_gashalo[i][j]-VelDispRz_gashalo[i][j])/(list_RplusdR[i]-list_R[i]);
	      else
		VelDispPhi_gashalo[i][j]=0;

	      VelDispRz_gashalo[i][j]/=rho;
	    }
	  else
	    VelDispRz_gashalo[i][j]=VelDispPhi_gashalo[i][j]=0;
	
	  VelVc2_gashalo[i][j]=R*Dphi_R[i][j];
	  
	  VelDispPhi_gashalo[i][j] += VelVc2_gashalo[i][j]+VelDispRz_gashalo[i][j];
	  
	  /* TEST */
	  //VelDispPhi_gashalo[i][j] = R*Dphi_R[i][j] + z*Dphi_z[i][j];

	  VelStreamPhi_gashalo[i][j] = 3.*halo_spinfactor * sqrt(VelVc2_gashalo[i][j]);
  
	  VelDispPhi_gashalo[i][j]-=VelStreamPhi_gashalo[i][j]*VelStreamPhi_gashalo[i][j];


	  if(VelDispRz_gashalo[i][j]<0)
	    VelDispRz_gashalo[i][j]=0;

	  if(VelDispPhi_gashalo[i][j]<0)
	    VelDispPhi_gashalo[i][j]=0;
	}
    }

  
  printf("done.\n"); fflush(stdout);  
}




double comp_Dphi_z_gashalo(double R,double z)
{
  double r,M_r;

  r=sqrt(R*R+z*z);
  
  M_r=gashalo_mass(r);

  if((r>0)&&(M_r>0))
    return G*z/(r*r*r)*M_r;
  else
    return 0;
}


double comp_Dphi_R_gashalo(double R,double z)
{
  double r,M_r;

  r=sqrt(R*R+z*z);

  M_r=gashalo_mass(r);

  if((r>0)&&(M_r>0))
    return G*R/(r*r*r)*M_r;
  else
    return 0;
}


double comp_rho_gashalo(double R,double z)
{
  double r,x;
  double gashalo_rho(double);

  r=sqrt(R*R+z*z);
  x = gashalo_rho(r);

  return x;
}

/*
double gashalo_m_lessthan_r(double r)
{
   return RHO_0*qromb(gashalo_profile_int,0,r);
}
*/


double gashalo_profile(double r)
{
  /* Beta Model */
  double RC=GasHalo_Rc_over_Rs*RS;
  return RHO_0*pow(1+(r/RC)*(r/RC),-1.5*GasHalo_Beta);

  /* Hydrostatic eq. - Isothermal 
  return pow(1+r, 1/r);
  */

  /* 
  Hydrostatic eq. - T function of r
  double alpha= 1.0;
  return exp(-(3.0/alpha)*r*(1-log(1+r)/r)/(log(1+r)-r/(1+r)));
  */

}
