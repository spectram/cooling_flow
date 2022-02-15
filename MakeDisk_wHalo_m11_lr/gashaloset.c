#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"
#include "prototypes.h"
#include "globvars.h"





void set_gashalo_velocities(void)
{
  int i,iz,ir;
  double r,R,z,ur,uz,vx,vy,vz,vdisp_rz,vdisp_phi,vstream_phi,vr,vphi,vc2tot;

  if(N_GASHALO==0) return;
  printf("set gashalo velocities...\t"); fflush(stdout);  
  
  for(i=1;i<=N_GASHALO;i++)
    {
      R=sqrt(xp_gashalo[i]*xp_gashalo[i] + yp_gashalo[i]*yp_gashalo[i]);
      z=zp_gashalo[i];

      ir=(int)( log(R/LL*(pow(FR,RSIZE)-1)+1)/log(FR));
      ur=( log(R/LL*(pow(FR,RSIZE)-1)+1)/log(FR)) - ir;

      iz=(int)( log(fabs(z)/LL*(pow(FZ,ZSIZE)-1)+1)/log(FZ));
      uz=( log(fabs(z)/LL*(pow(FZ,ZSIZE)-1)+1)/log(FZ)) - iz;

      vdisp_rz= VelDispRz_gashalo[ir][iz]*(1-ur)*(1-uz)
               +VelDispRz_gashalo[ir+1][iz]*(ur)*(1-uz)
   	       +VelDispRz_gashalo[ir][iz+1]*(1-ur)*(uz) 
               +VelDispRz_gashalo[ir+1][iz+1]*(ur)*(uz);

      vdisp_phi=VelDispPhi_gashalo[ir][iz]*(1-ur)*(1-uz)
	       +VelDispPhi_gashalo[ir+1][iz]*(ur)*(1-uz)
	       +VelDispPhi_gashalo[ir][iz+1]*(1-ur)*(uz) 
               +VelDispPhi_gashalo[ir+1][iz+1]*(ur)*(uz);

      vstream_phi=VelStreamPhi_gashalo[ir][iz]*(1-ur)*(1-uz)
	       +VelStreamPhi_gashalo[ir+1][iz]*(ur)*(1-uz)
	       +VelStreamPhi_gashalo[ir][iz+1]*(1-ur)*(uz) 
               +VelStreamPhi_gashalo[ir+1][iz+1]*(ur)*(uz);
       
       double tmp_r = Dphi_R[ir][iz]*(1-ur)*(1-uz)
	       +Dphi_R[ir+1][iz]*(ur)*(1-uz)
	       +Dphi_R[ir][iz+1]*(1-ur)*(uz) 
               +Dphi_R[ir+1][iz+1]*(ur)*(uz);
       double tmp_z = Dphi_z[ir][iz]*(1-ur)*(1-uz)
	       +Dphi_z[ir+1][iz]*(ur)*(1-uz)
	       +Dphi_z[ir][iz+1]*(1-ur)*(uz) 
               +Dphi_z[ir+1][iz+1]*(ur)*(uz);
       
    vc2tot = R*fabs(tmp_r) + fabs(z)*fabs(tmp_z);
      
      if(vdisp_rz<0)
	{
	  printf("in gashalo: vdisp_rz:%g   %g %g %d %d \n",vdisp_rz,ur,uz,ir,iz);
	  vdisp_rz=-vdisp_rz;
	}
      if(vdisp_phi<0)
	{
	  printf("in gashalo: vdisp_phi:%g  %g %g %d %d\n",vdisp_phi,ur,uz,ir,iz);
	  
	  vdisp_phi=-vdisp_phi;
	}

      /* --------------------------------------------
	    Set Gas Velocity
	 -------------------------------------------- */

      vr=0.0;
      vz=0.0;
      vphi=vstream_phi;
      vx=vr*xp_gashalo[i]/R - vphi*yp_gashalo[i]/R;
      vy=vr*yp_gashalo[i]/R + vphi*xp_gashalo[i]/R;

      vxp_gashalo[i]=vx;
      vyp_gashalo[i]=vy;
      vzp_gashalo[i]=vz;

      /* --------------------------------------------
	    Set Gas Temperature
	 -------------------------------------------- */

      // Xiangcheng: gas energy follows hydrostatic equilibrium
      u_gashalo[i] = (vdisp_rz-vstream_phi*vstream_phi/3.) / (GAMMA-1.);
      if (u_gashalo[i]<0.) u_gashalo[i] = 0.;

      //u_gashalo[i] = vdisp_rz / (GAMMA-1.);

      if(u_gashalo[i] > 1.e10)
      {
      printf("WARNING: LL=%g FR/FZ=%g/%g RSIZE/ZSIZE=%d/%d u=%g vd_rz=%g vd_rp=%g vstrm=%g x/y/z=%g/%g/%g vx/vy/vz=%g/%g/%g ir/ur/iz/uz=%d/%g/%d/%g \n",
        LL,FR,FZ,RSIZE,ZSIZE,u_gashalo[i],vdisp_rz,vdisp_phi,vstream_phi,xp_gashalo[i],yp_gashalo[i],zp_gashalo[i],
        vx,vy,vz,ir,ur,iz,uz);
      }
      if(u_gashalo[i]<=1.0) u_gashalo[i]=1.0;

    }
  printf("done.\n"); fflush(stdout);
}




void set_gashalo_positions(void)
{
  int   i,countr,countz;
  double q,R,phi,theta;
  //double AVOIDANCE_ZONE = 1.0;
  double AVOIDANCE_ZONE = 0.1;

  if(N_GASHALO==0) return;

  srand48(GHRAND);

  printf("set gashalo positions...\t"); fflush(stdout);

  for(i=1,countr=countz=0;i<=N_GASHALO;)
    {
      do
	{
	  q=drand48();
	  R=gashalo_q_to_r(q);
	}
      while((R > Rvir)||(R>LL));

      phi=drand48()*PI*2;
      theta=acos(drand48()*2-1);
	  
      xp_gashalo[i]=R*sin(theta)*cos(phi);
      yp_gashalo[i]=R*sin(theta)*sin(phi);
      zp_gashalo[i]=R*cos(theta);

      /* avoid the central disk, since it screws things up and the contribution 
	   to the mass there is negligible anyways */
      if((R < AVOIDANCE_ZONE*H) && (fabs(zp_gashalo[i]) < AVOIDANCE_ZONE*Z0))
	{
	  printf("\nrejected i=%d R=%g z=%g H=%g Z0=%g AV=%g.... ",
	    i,R,zp_gashalo[i],H,Z0,AVOIDANCE_ZONE); fflush(stdout);
	  continue;
	}

      i++;
    }


  for(i=1;i<=N_GASHALO;i++)
    mp_gashalo[i]=M_GASHALO/N_GASHALO;

  printf("done\n"); fflush(stdout);

}




