

double set_halo_positions(void);
double set_disk_positions(void);
double set_bulge_positions(void);
double set_gas_positions(void);
void set_gashalo_positions(void);


double set_halo_velocities(void);
double set_disk_velocities(void);
double set_bulge_velocities(void);
double set_gas_velocities(void);
void set_gashalo_velocities(void);

void setup_massprofile(int mode);
double vc2_sph_function(double r);


void read_parameterfile(char *fname);
void save_particles(char *fname);

double mass_total_dark(double radius);
double mass_cumulative_bulge(double radius);


double comp_Dphi_R_constant(double R,double z);
double comp_Dphi_z_constant(double R,double z);


double comp_Dphi_z_BH(double R,double z);
double comp_Dphi_R_BH(double R,double z);


double comp_Dphi_z_halo(double R,double z);
double comp_Dphi_R_halo(double R,double z);
double comp_rho_halo(double R,double z);


double comp_Dphi_z_bulge(double R,double z);
double comp_Dphi_R_bulge(double R,double z);
double comp_rho_bulge(double R,double z);


double comp_Dphi_z_disk(double R,double z);
double comp_Dphi_R_disk(double R,double z);
double comp_rho_disk(double R,double z);
double comp_Dphi_R_disk_razorthin(double RR,double zz);


double comp_Dphi_z_gas(double R,double z);
double comp_Dphi_R_gas(double R,double z);
double comp_rho_gas(double R,double z);
double comp_Dphi_R_gas_razorthin(double RR,double zz);


double comp_Dphi_z_gashalo(double R,double z);
double comp_Dphi_R_gashalo(double R,double z);
double comp_rho_gashalo(double R,double z);
void compute_velocity_dispersions_gashalo(void);


double comp_Dphi_z(double R,double z);
double comp_Dphi_R(double R,double z);


double epicyclic_kappa2(double R);


double splint_zs_y_D2y(double z);
double splint_xl_yl_D2yl(double t);


double   drand48(void);

double halogas_profile_int(double r);
void setup_halogas_massprofile(void);
double gashalo_q_to_r(double q);

double mass_cumulative_disk(double);
double mass_cumulative_bulge(double);
double mass_cumulative_gas(double);
double mass_cumulative_gashalo(double R);
double gashalo_mass(double r);
double gashalo_profile(double r);



