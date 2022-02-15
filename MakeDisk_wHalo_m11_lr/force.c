#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"


#include "prototypes.h"
#include "globvars.h"

 




compute_force_field()
{
  int i,j;
  double k2,dphi_R_dr;


  printf("Start computing force field.\n");

  /* --- Don't know if I want this, but it stops the core dump for now --- */
  for(j=0;j<=ZSIZE;j++) Dphi_R[0][j] = 0;

  for(i=1;i<=RSIZE;i++)
    {
      printf("force field %d(%d)\n",i,RSIZE);  fflush(stdout);

      for(j=0;j<=ZSIZE;j++)
	{
	  if(j==0)
	    {
	      Dphi_z[i][j] = 0;
	      Dphi_z_dR[i][j] =0;
	    }
	  else
	    {
	      Dphi_z[i][j]    = comp_Dphi_z( list_R[i], list_z[j] );
	      Dphi_z_dR[i][j] = comp_Dphi_z( list_RplusdR[i] , list_z[j] );
	    }
	  /* --- Never actually does this check, since i starts at 1 above --- */
	  if(i==0)
	    Dphi_R[i][j] = 0;
	  else
	    Dphi_R[i][j] = comp_Dphi_R( list_R[i], list_z[j]);
	}

    }


  
  for(j=0;j<=ZSIZE;j++)
    {
      Dphi_z[0][j]    = Dphi_z[1][j];
      Dphi_z_dR[0][j] = Dphi_z_dR[1][j];
    }


  for(i=1,epi_gamma2[0]=1;i<=RSIZE;i++)
    {
      dphi_R_dr = comp_Dphi_R( list_RplusdR[i], 0);
      
      k2=3/list_R[i] *  Dphi_R[i][0] + (dphi_R_dr - Dphi_R[i][0]) /(list_RplusdR[i]- list_R[i]);
      
      epi_gamma2[i]=4/list_R[i]*Dphi_R[i][0]/k2;
      epi_kappa2[i]=k2;
    }

  epi_kappa2[0]=epi_kappa2[1];

  printf("Force field finished.\n");
}






double comp_Dphi_z(double R,double z)
{
  /*
  printf("R=%g z=%g mgh=%g DPHI_DZ_HALO=%g, DPHI_DZ_BUL=%g, DPHI_DZ_DSK=%g DPHI_DZ_BH=%g DPHI_DZ_GAS=%g DPHI_DZ_GH=%g \n",
    R,z,gashalo_mass(R),comp_Dphi_z_halo(R,z),comp_Dphi_z_bulge(R,z),comp_Dphi_z_disk(R,z),
    comp_Dphi_z_BH(R,z),comp_Dphi_z_gas(R,z),comp_Dphi_z_gashalo(R,z)); fflush(stdout);
  */
  return comp_Dphi_z_constant(fabs(R),fabs(z))
		// +comp_Dphi_z_halo(fabs(R),fabs(z))
        +comp_Dphi_z_bulge(fabs(R),fabs(z)) 
        +comp_Dphi_z_disk(fabs(R),fabs(z))
        +comp_Dphi_z_BH(fabs(R),fabs(z))
	    +comp_Dphi_z_gas(fabs(R),fabs(z));
	    //+comp_Dphi_z_gashalo(fabs(R),fabs(z))

}



double comp_Dphi_R(double R,double z)
{
  /*
  printf("R=%g z=%g mgh=%g DPHI_DR_HALO=%g, DPHI_DR_BUL=%g, DPHI_DR_DSK=%g DPHI_DR_BH=%g DPHI_DR_GAS=%g DPHI_DR_GH=%g DPHI_DR_constant=%g\n",
    R,z,gashalo_mass(R),comp_Dphi_R_halo(R,z),comp_Dphi_R_bulge(R,z),comp_Dphi_R_disk(R,z),
    comp_Dphi_R_BH(R,z),comp_Dphi_R_gas(R,z),comp_Dphi_R_gashalo(R,z),comp_Dphi_R_constant(R,z)); fflush(stdout);
  */
  return comp_Dphi_R_constant(fabs(R),fabs(z))
//		+comp_Dphi_R_halo(fabs(R),fabs(z))
        +comp_Dphi_R_bulge(fabs(R),fabs(z)) 
        +comp_Dphi_R_disk(fabs(R),fabs(z))
        +comp_Dphi_R_BH(fabs(R),fabs(z))
	    +comp_Dphi_R_gas(fabs(R),fabs(z));
//	    +comp_Dphi_R_gashalo(fabs(R),fabs(z))

}


double comp_Dphi_R_constant(double RR,double zz)
{
  double r;

  r=sqrt(RR*RR+zz*zz);

  if (r>0)
    return VC*VC*RR/(r*r);
  else
    return 0;
}


double comp_Dphi_z_constant(double RR,double zz)
{
  double r;

  r=sqrt(RR*RR+zz*zz);

  if (r>0)
    return VC*VC*zz/(r*r);
  else
    return 0;
}
