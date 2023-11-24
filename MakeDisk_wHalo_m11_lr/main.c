


#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"
#include "globvars.h"
#include "prototypes.h"

#define RES_FACTOR 1  /* should equal 8e4/desired resolution. if disk/bulge/gas disk mass changes, change N_DISK / N_Gas / N_bulge accordingly  */
#define f_gas 0.2    /* should equal gas fraction in disk  */

int main(int argc,char *argv[])
{
  int i;
  char filename[100];
  FILE *fd;
  int write_header(FILE *fd);

  /*******************************************/

  VC		        =     200;                          /* for analytic isothermal profile */
  CC                    =     10.;                          /* not used */
  Mvir                  =     100;              	    /* not used */
  LAMBDA                =   0.033;       		    /* not used          */
  M_HALO                =      84;		            /* not used */
  M_DISK                =     6;   		            /* total disk stellar mass in units of 10^10 Msolar */
  M_GAS                 =     M_DISK*f_gas;   		    /* total disk gas mass in units of 10^10 Msolar */
  M_GASHALO             =     2.48;              	    /* not used */
  M_BULGE               =     0.02;   		            /* total bulge mass in units of 10^10 Msolar */
  M_BH                  =     0.00002;			    /* seed BH mass in units of 10^10 Msolar */

  DARKMASS_IN_ROPT      =    -10.0;	                    /* not used  */
  H                     =     3.75;              	    /* manually set radial disk scale length, this decouples
                           		                       the angular momentum of the disk and the halo and hence
                           		                       the spin doesn't set the size of H and Jd is not
                           		                       fixed.   */
  DiskHeight            =      0.1;    	                    /* thickness of disk in units of radial scale length */
  N_HALO                =        0;                         /* should be 0, analytic gravity is used instead */
  N_DISK                =     M_DISK*125000*RES_FACTOR;     /* desired number of collisionless particles in disk, 8e4 for RES_FACTOR=1 */
  N_GAS                 =     N_DISK*f_gas;                 /* number of gas particles in disk, 8e4 for RES_FACTOR=1 */
  N_BULGE               =     M_BULGE*125000*RES_FACTOR;    /* number of bulge particles, 8e4  for RES_FACTOR=1  */
  N_GASHALO             =        0;                         /* should be 0, added later */
  HUBBLE                =        1;			    /* Hubble parameter (1 means units of h-1)*/
  Z                     =        0;		   	    /* Redshift of Galaxy */
  GasDistribution       =        0;			    /* 0 = exp. (normal, same Rd as disk)
				   			       1 = exp. (with Rd -> Rd*Alpha)
				   			       2 = Power Law (with PowerLawGamma < 2) 
				   			       3 = exp. (with Rd -> Rd*Alpha) 
				   				         with cutoff inside HoleRadius */
  GasExpAlpha           =        1;                          /* gas is exp. with Rd*Alpha scale length */
  PowerLawGamma         =        1;               	    /* power-law index, sigma ~ r^-gamma   (gamma=1 is mestel) - must be < 2 */
  PowerLawCutOff        =        4;			    /* in units of Rd, when gas disk is terminated */
  HoleRadius            =     0.03;			    /* in units of Rd, when inner gas disk is cut off */
  HoleGamma             =      2.0;         	            /* speed with which that cutoff occurs (power-law slope of cutoff >0) */
  SetInitModeAmp        =      0.0;  	  	            /* amplitude (fractional) of 'seed' mode M in initial gas */
  SetInitModeCut        =      0.0;   		            /* cutoff radius (absolute) of 'seed' mode M in initial gas */
  SetInitModeM          =        1;			    /* modenumber m for 'seed' mode in initial gas */
  BulgeSize             =   0.05*2;     	            /* bulge scale length in ABSOLUTE UNITS (*NOT* units of R_d) */
  BulgeDistribution     =        0;    	                    /* 0 = Hernquist profile (BulgeSize sets a)
                                	                       1 = Spherical exp. (BulgeSize sets 3D Rd) */
  GasHalo_Rc_over_Rs    =        1;                         /* 0.5, ratio of gas halo core to halo scale radius */
  GasHalo_Beta          =      0.5;                         /* 1.33, beta value for the beta-profile */
  HI_GasMassFraction    =      0.0;	                    /* in terms of the total gas mass */
  HI_GasDiskScaleLength =     3.72;	                    /* in terms of scale length of the disk */ 
  Qstabilizefactor      =      1.0;			    /* Toomre Q for initializing the radial dispersion in 
  							       the stellar disk (gas disk pressure set by DiskHeight) */
  ColdFlowExtraMgas     =      0.0;		            /* extra gas mass to assign to a cold flow
  									infalling on the disk */
  ColdFlowMdot          =      2.0;		            /* inflow rate through cold flow [Msun/yr] */
  ColdFlowTheta         =     60.0;	       	            /* angle (theta) of flow direction */
  ColdFlowPhi           =      5.0;		            /* angle (phi) of flow direction */
  ColdFlowRinner        =     10.0;		            /* innermost radius of the cold flow (to start) */
  ColdFlowRadius        =      9.0;		            /* radius [kpc] of cold flow cross section */


  /**********************************************************/

  MD= (M_DISK+M_GAS)/Mvir;
  MB= M_BULGE/Mvir;
  GasFraction= M_GAS/(M_DISK+M_GAS);
  
  N_GAS_DISK = N_GAS;
  N_GAS_FLOW = (int)(ColdFlowExtraMgas/(M_GAS/((float)N_GAS_DISK)));
  N_GAS = N_GAS_DISK + N_GAS_FLOW;

  /**********************************************************/


  if(argc!=2)
    {
      fprintf(stderr,"\n\nwrong argument(s).  Specify an output filename.\n\n");
      exit(0);
    }
  strcpy(filename, argv[1]);


  init_units();        /* set system of units */
  structure();         /* determine structure of halo, disk, and bulge */
  init();              /* allocate arrays */

  set_halo_positions();
  set_disk_positions();
  set_bulge_positions();
  set_gas_positions();
  set_gashalo_positions();
  
  compute_force_field();

  compute_velocity_dispersions_halo();  
  compute_velocity_dispersions_disk();
  compute_velocity_dispersions_bulge();  
  compute_velocity_dispersions_gas();
  compute_velocity_dispersions_gashalo();
 
  compute_local_escape_speed();

  set_halo_velocities();
  set_disk_velocities();    
  set_bulge_velocities();
  set_gas_velocities();
  set_gashalo_velocities();


  
  save_particles(filename);


  strcat(filename, ".parameters");
  if(fd=fopen(filename,"w"))
    {
      printf("writing circular velocity curve + Toomre's Q\n");
      write_header(fd);
      plot_circular_speeds(fd);
      plot_toomre_stability(fd);
      fclose(fd);
    }
  else
    {
      fprintf(stderr,"Can't open file '%s'.\n",filename);
      exit(0);
    }
  printf("done.\n");


  printf("Disk scale length: %g\n",H);
  printf("Rvir: %g\n",Rvir);
}



int write_header(FILE *fd)
{
  fprintf(fd,"GAL1    \n");
  fprintf(fd,"c       \t%g\n",CC);
  fprintf(fd,"mvir    \t%g\n",Mvir);
  fprintf(fd,"lambda  \t%g\n",LAMBDA);
  fprintf(fd,"mdisk   \t%g\n",M_DISK);
  fprintf(fd,"mbulge  \t%g\n",M_BULGE);
  fprintf(fd,"mgas    \t%g\n",M_GAS);
  fprintf(fd,"mbh     \t%g\n",M_BH);
  fprintf(fd,"disklen \t%g\n",H);
  fprintf(fd,"jd      \t%g\n",JD);
  fprintf(fd,"f_gas   \t%g\n",GasFraction);
  fprintf(fd,"disk_ht \t%g\n",DiskHeight);
  fprintf(fd,"bulge_a \t%g\n",BulgeSize);
  fprintf(fd,"n_halo  \t%d\n",N_HALO);
  fprintf(fd,"n_disk  \t%d\n",N_DISK);
  fprintf(fd,"n_gas   \t%d\n",N_GAS);
  fprintf(fd,"n_gashal\t%d\n",N_GASHALO);
  fprintf(fd,"n_bulge \t%d\n",N_BULGE);
  fprintf(fd,"hubble  \t%g\n",HUBBLE);
  fprintf(fd,"redshft \t%g\n",Z);
  fprintf(fd,"GasDist \t%d\n",GasDistribution);
  fprintf(fd,"GasExpA \t%g\n",GasExpAlpha);
  fprintf(fd,"PLGamma \t%g\n",PowerLawGamma);
  fprintf(fd,"PLCutO  \t%g\n",PowerLawCutOff);
//  fprintf(fd,"HoleR   \t%g\n",HoleRadius);
//  fprintf(fd,"HoleGam \t%g\n",HoleGamma);
  fprintf(fd,"BDist   \t%d\n",BulgeDistribution);
//  fprintf(fd,"m_iAmp  \t%g\n",SetInitModeAmp);
//  fprintf(fd,"m_Rcut  \t%g\n",SetInitModeCut);
  fprintf(fd,"Rvir    \t%g\n",Rvir);
  fprintf(fd,"Vvir    \t%g\n",Vvir);
  fprintf(fd,"Rd      \t%g\n",H);
  fprintf(fd,"Q       \t%g\n",Qstabilizefactor);
  fprintf(fd,"HI_gmf  \t%g\n",HI_GasMassFraction);
  fprintf(fd,"HI_gds  \t%g\n",HI_GasDiskScaleLength);
  fprintf(fd,"junk\n");
}



