#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sim5lib.h"
//#include "sim5lib.c"

static char model_name[] = "Novikov-Thorne";
static double bh_mass = 0.0;
static double bh_spin = 0.0;


int diskmodel_init(double M, double a, char *params)
{
    bh_mass = M;
    bh_spin = a;
    double alpha = 0.1;
    double mdot  = 0.1;
    char* value;

    // read alpha parameter or use default value
    if ((value=key_value_get(params,"alpha")) != NULL) alpha = atof(value);

    // read mdot parameter or use default value
    if ((value=key_value_get(params,"mdot")) != NULL) mdot = atof(value);

    // read mdot parameter or use default value
    if ((value=key_value_get(params,"mdot_abs")) != NULL) mdot = atof(value)/bh_mass/Mdot_Edd;

    // read luminosity parameter or use default value for mdot
    if ((value=key_value_get(params,"lumi")) != NULL) mdot = -atof(value);

    return disk_nt_setup(bh_mass, bh_spin, fabs(mdot), alpha, (mdot>0)?0:DISK_NT_OPTION_LUMINOSITY);
}



void diskmodel_done()
// model finalization
{
    disk_nt_done();
}


char* diskmodel_name()
// model name
{
    return model_name;
}


void diskmodel_help()
// model help
{
    fprintf(stderr, "disk model '%s' options (orderd with increasing priority):\n", model_name);
    fprintf(stderr, "alpha    ... alpha-viscosity value (defaults to 0.1)\n");
    fprintf(stderr, "mdot     ... mass accretion rate [Mdot_Edd] (defaults to 0.1)\n");
    fprintf(stderr, "mdot_abs ... mass accretion rate [g/s]\n");
    fprintf(stderr, "lumi     ... luminosity (instead of mdot) [L_Edd]\n");
}


double diskmodel_r_min()
{
    return disk_nt_r_min();
}


double diskmodel_flux(double r)
{
    return disk_nt_flux(r);
}


double diskmodel_lumi()
{
    return disk_nt_lumi();
}


double diskmodel_mdot()
{
    return disk_nt_mdot();
}



double diskmodel_sigma(double r)
{
    return disk_nt_sigma(r);
}



double diskmodel_l(double r)
// specific angular momentum
// result is dimensionless (for a unit mass)
{
    return disk_nt_ell(r);
}



double diskmodel_vr(double r)
// radial velocity
// result in units of sound speed
{
    return disk_nt_vr(r);
}



double diskmodel_h(double r)
// disk height
// result in GM/c2 units
{
    return disk_nt_h(r);
}



double diskmodel_dhdr(double r)
// slope of the disk surface
// result is dimensionless
{
    return disk_nt_dhdr(r);
}


double diskmodel_eval(double r, int quantity)
{
    return 0.0;
}


void diskmodel_params(FILE* output)
// model name
{
    fprintf(output, "# model = %s\n", model_name); 
    fprintf(output, "# M     = %.3f [M_sun]\n", bh_mass); 
    fprintf(output, "# a     = %.3f\n", bh_spin); 
    fprintf(output, "# r_bh  = %.3f [GM/c2]\n", r_bh(bh_spin)); 
    fprintf(output, "# r_ms  = %.3f [GM/c2]\n", r_ms(bh_spin)); 
    fprintf(output, "# rmin  = %.3f [GM/c2]\n", diskmodel_r_min()); 
    fprintf(output, "# lumi  = %.5e [Ledd]\n", diskmodel_lumi()); 
    fprintf(output, "# mdot  = %.5e [Mdot_Edd] (%.5e x 10^18 g/s)\n", diskmodel_mdot(), diskmodel_mdot()*bh_mass*Mdot_Edd/1e18); 
}


void diskmodel_dump()
{
    printf("# Disk model dump\n"); 
    printf("#-------------------------------------------\n"); 
    printf("# model = %s\n", model_name); 
    printf("# M     = %.3f [M_sun]\n", bh_mass); 
    printf("# a     = %.3f\n", bh_spin); 
    printf("# r_bh  = %.3f [GM/c2]\n", r_bh(bh_spin)); 
    printf("# r_ms  = %.3f [GM/c2]\n", r_ms(bh_spin)); 
    printf("# rmin  = %.3f [GM/c2]\n", diskmodel_r_min()); 
    printf("# rmax  = %.3f [GM/c2]\n", 2000.); 
    printf("# lumi  = %.5e [Ledd]\n", diskmodel_lumi()); 
    printf("# mdot  = %.5e [Mdot_Edd] (%.5e x 10^18 g/s)\n", diskmodel_mdot(), diskmodel_mdot()*bh_mass*Mdot_Edd/1e18); 
    printf("#-------------------------------------------\n"); 
    printf("# col1: radius [GM/c2]\n");
    printf("# col2: flux (one side) [erg s-1 cm-2]\n");
    printf("# col3: sigma [g cm-2]\n");
    printf("# col4: specific angular momentum [none]\n");
    printf("# col5: radial velocity [speed of light]\n");
    printf("# col6: disk height [GM/c2]\n");
    printf("# col7: disk slope (derivative dH/dR of height with respect to equatorial radius) [none]\n");
    printf("#-------------------------------------------\n"); 
    double r;
    for (r=r_bh(bh_spin); r<2000.; r*=1.05) {
        printf("%e  %e  %e  %e  %+e  %+e  %+e\n",
            r, 
            diskmodel_flux(r),
            diskmodel_sigma(r),
            diskmodel_l(r),
            diskmodel_vr(r),
            diskmodel_h(r),
            diskmodel_dhdr(r)
        );
    }
}




#ifdef DISK_NT_DEBUG    
int main(int argc, char *argv[])
{
    double mdot = 0.1;
    double spin = 0.0;
    
    if (argc > 1) spin = atof(argv[1]);
    if (argc > 2) mdot = atof(argv[2]);
    
    printf("# model: %s\n", model_name);
    printf("# mdot = %.4f\n", mdot);
    printf("# spin = %.4f\n", spin);
    
    char options[256];
    sprintf(options, "mdot=%e", mdot);
    
    diskmodel_init(10.0, spin, options);
    diskmodel_dump();
    return 0;
}
#endif


