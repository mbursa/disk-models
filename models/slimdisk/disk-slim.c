#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sim5lib.h"
#include "disk-slim-dsi.h"



static char model_name[] = "slim-disk";
static double bh_mass = 0.0;
static double bh_spin = 0.0;
static double disk_mdot = 0.0;
static double disk_alpha = 0.0;
static double r_h, r_min, r_max;



#define disk_alpha_min   0.005
#define disk_alpha_max   0.100

#define disk_mdot_min   1e-02
#define disk_mdot_max   5e+02


// fwd declarations
static double find_mdot_for_luminosity(char* path, double spin, double alpha, double L0);
static double qspline(double x1, double y1, double dydx1, double x2, double y2, double dydx2, double x);



int diskmodel_init(double M, double a, char *params)
{
    bh_mass = M;
    bh_spin = a;
    r_h = r_bh(a);
    disk_mdot  = 0.1;
    disk_alpha = 0.1;
    char* value;

    char path[256] = ".";

    // read path parameter or use default value
    if ((value=key_value_get(params,"path")) != NULL) {
        strncpy(path, value, sizeof(path)-1);
    }
    
    // read alpha parameter or use default value
    if ((value=key_value_get(params,"alpha")) != NULL) {
        disk_alpha = atof(value);

        if ((disk_alpha < disk_alpha_min) || (disk_alpha > disk_alpha_max)) {
            fprintf(stderr, "# ERROR (disk-sd::diskmodel_init): alpha parameter is out of range [%.2e..%.2e]\n", disk_alpha_min, disk_alpha_max);
            return -1;
        }
    }

    // read mdot parameter if provided
    if ((value=key_value_get(params,"mdot")) != NULL) {
        disk_mdot = atof(value);
    }

    // read luminosity parameter if provided and convert it to mdot
    if ((value=key_value_get(params,"lumi")) != NULL) {
        disk_mdot = find_mdot_for_luminosity(path, bh_spin, disk_alpha, atof(value));
    }

    // setup DSI model
    int dsi_setup_result = dsi_init(path, disk_mdot, bh_spin, disk_alpha, &r_min, &r_max);

    // check result
    if (dsi_setup_result != 0) {
        fprintf(stderr, "# ERROR (disk-sd::diskmodel_init): model could not initialize for given parameters (a=%.3e/mdot=%.3e/alpha=%.3e)\n", bh_spin, disk_mdot, disk_alpha);
        return -1;
    }

    return 0;
}



void diskmodel_done()
// model finalization
{
    dsi_done();
}


char* diskmodel_name()
// model name
{
    return model_name;
}


double diskmodel_r_min()
{
    return r_min;
}


double diskmodel_flux(double r)
// total outgoing flux
// note: r has the meaning of "equatorial-plane R", not the B-L coordinate
// result in [erg cm-2 s-1]
{
    if (r <= r_h+1e-3) return 0.0;
 
    double F;   
    if (r < r_max)
        F = DSI_MASS/bh_mass * 0.5*dsi_eval(r,DSI_COL_F);   // only one-side flux (factor 0.5)
    else {
        F = DSI_MASS/bh_mass * 0.5*dsi_eval(r_max,DSI_COL_F) * pow(r/(r_max), -3.);
    }

    return ((F>0.0)&&(!isnan(F))) ? F/gu_sb_sigma*sb_sigma : 0.0;
    // gu_sb_sigma*sb_sigma factor makes conversion from geometrical to cgs units
}


double diskmodel_lumi()
//! Gets total disk luminosity.
//! Luminosity is obtained by integrating local flux over the surface area _without_ relativistic corrections.
//! ```L = 2*2pi \int F r dr```
//! @result Total disk luminosity through both surfaces in Eddinton units.
{
    const float disk_rmax = 1e5;

    // integrate disk luminosity from r_ms to disk_rmax rg
    // - the integration uses 'logarithmic rule': L = \int f(x) dx \int f(x)*x d(log(x))
    // - it also makes conversion from local to coordinate flux
    double func_luminosity(double log_r)
    {
        double diskmodel_l(double);
        double r = exp(log_r);
        // calculate U_t
        double gtt = -1. + 2./r;
        double gtf = -2.*bh_spin/r;
        double gff = sqr(r) + sqr(bh_spin) + 2.*sqr(bh_spin)/r;
        double ell = diskmodel_l(r);
        double Omega = -(gtf + ell*gtt) / (gff + ell*gtf);
        double U_t = sqrt(-1.0/(gtt + 2.*Omega*gtf + sqr(Omega)*gff)) * (gtt + Omega*gtf);
        double F = diskmodel_flux(r);
        // dL = 2pi*r*F(r) dr, extra r comes from log integration
        return 2.*M_PI*r*2.0*(-U_t)*F * r;
    }

    double L = integrate_simpson(func_luminosity, log(r_bh(bh_spin)+1e-2), log(disk_rmax), 1e-5);

    // fix units to erg/s
    L *= sqr(bh_mass*grav_radius);

    return L/(L_Edd*bh_mass);
}


double diskmodel_mdot()
{
    return disk_mdot;
}



double diskmodel_temp(double r)
// temperature
// note: r is r is the "equatorial-plane r", not the B-L one
// result in kelvins
{
    return sqrt4(diskmodel_flux(r)/sb_sigma);
}



double diskmodel_sigma(double r)
// surface density
// note: r is r is the "equatorial-plane r", not the B-L one
// result in g/cm2
{
    // approx. mass scaling: sigma~1/v_radial, v_r~Omega*R => no M scaling
    if (r>0.9*r_max) r=0.9*r_max;
    double S = 0.5*dsi_eval(r,DSI_COL_SIGMA)/7.42474e-27;  // return only midplane column density (factor 0.5)

    if (disk_alpha < disk_alpha_min) S /= disk_alpha/disk_alpha_min;

    return ((S>0.0)&&(!isnan(S))) ? S : 0.0;
}



double diskmodel_l(double r)
// specific angular momentum U_f/U_t
// note: r is r is the "equatorial-plane r", not the B-L one
// result is dimensionless (for a unit mass)
{
// original radial model:
//    r = fmax(r, r_min);
//    if (r>r_max) return dsi_eval(r_max,DSI_COL_ELL)/(DSI_MASS*grav_radius_cgs)*pow(r/r_max,0.5);
//    return dsi_eval(r,DSI_COL_ELL)/(DSI_MASS*grav_radius_cgs);

// polytropic model
    r = fmax(r, r_min+1e-3); // for r<r_min assume that ell is constant
    if (r<r_max)
        return dsi_eval(r,DSI_COL_ELL);
    else
        return ellK(r,bh_spin);//dsi_eval(r_max,DSI_COL_ELL)*pow(r/r_max,0.5);
}



double diskmodel_vr(double r)
// radial velocity
// note 1: r is the "equatorial-plane r", not the B-L one;
// result in units of light speed
{
    // note: DSI returns a positive value, we must make it negative as it represents inward velocity

    // approx. mass scaling: v_r~Omega*R => no M scaling
    if (r<r_min) {
        double r1 = 1.0*r_min;
        double r2 = 1.1*r_min;
        double v1 = dsi_eval(r1,DSI_COL_VR);
        double v2 = dsi_eval(r2,DSI_COL_VR);
        return -fmin(0.9999, qspline(r_h, 0.99, 2.0, r1, v1, (v2-v1)/(r2-r1), r));
    }
    if (r>r_max) return -dsi_eval(r_max,DSI_COL_VR);  // keep the last value (should be sufficiently low)
    return -fmin(dsi_eval(r,DSI_COL_VR), 0.9999);
}



double diskmodel_h(double r)
// disk height
// note: r is r is the "equatorial-plane r", not the B-L one
// result in GM/c2 units
{
//    // original dsi
//    if (r<r_min) return r*dsi_eval(r_min,DSI_COL_HR);
//    if (r>.8*r_max) return r*dsi_eval(.8*r_max,DSI_COL_HR);
//    return r*dsi_eval(r, DSI_COL_HR);

    // dsi(polytropic) returns H in meters => rescaling by [DSI_MASS*si_grav_radius] is needed

    if (r<r_min) {
        //return dsi_eval(r_min,DSI_COL_HR)/(DSI_MASS*si_grav_radius) * r/r_min;
        double r1 = 1.00*r_min;
        double r2 = 1.02*r_min;
        double h1 = dsi_eval(r1,DSI_COL_HR)/(DSI_MASS*si_grav_radius);
        double h2 = dsi_eval(r2,DSI_COL_HR)/(DSI_MASS*si_grav_radius);
        return qspline(r_h, 2e-3, 0.0, r1, h1, (h2-h1)/(r2-r1), r);
    }

    if (r>0.8*r_max) return dsi_eval(0.8*r_max,DSI_COL_HR)/(DSI_MASS*si_grav_radius) * r/(0.8*r_max); // keep H/R
    return dsi_eval(r, DSI_COL_HR)/(DSI_MASS*si_grav_radius);
}



double diskmodel_dhdr(double r)
// slope of the disk surface
// note: r is r is the "equatorial-plane r", not the B-L one
// result is dimensionless
{
    double H1,H2;
    double r1,r2;

    if (r>0.80*r_max) {
        // for r>0.8*r_max H/R is kept constant => dH/dR = H/R
        return diskmodel_h(0.8*r_max)/(0.80*r_max);

        // old scheme: extrapolate, keeping the derivative
        //r1 = 0.78*r_max;
        //r2 = 0.80*r_max;
        //H1 = flux_model_eval_h(r1);
        //H2 = flux_model_eval_h(r2);
        //return (H2-H1)/(r2-r1);
    }

    r1 = fmax(0.98*r, r_h);
    r2 = 1.02*r;
    H1 = diskmodel_h(r1);
    H2 = diskmodel_h(r2);
    double D1 = (H2-H1)/(r2-r1);

    r1 = fmax(0.94*r, r_h);
    r2 = 1.06*r;
    H1 = diskmodel_h(r1);
    H2 = diskmodel_h(r2);
    double D2 = (H2-H1)/(r2-r1);

    r1 = fmax(0.90*r, r_h);
    r2 = 1.10*r;
    H1 = diskmodel_h(r1);
    H2 = diskmodel_h(r2);
    double D3 = (H2-H1)/(r2-r1);

    return sqrt3(D1*D2*D3); // geometrical mean
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
    printf("# rmax  = %.3f [GM/c2]\n", r_max);
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
    for (r=r_bh(bh_spin); r<2.*r_max; r*=1.005) {
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



static double find_mdot_for_luminosity(char* path, double spin, double alpha, double L0) {
    double mdot;

    //fprintf(stderr,"looking for L0=%.2e\n", L0);

    double fce(double xmdot) {
        double L;
        dsi_init(path, xmdot, spin, alpha, &r_min, &r_max);
        L = diskmodel_lumi();
        dsi_done();
        return L0-L;
    }

    int res = rtbis(disk_mdot_min, disk_mdot_max, 1e-6, fce, &mdot);
    return (res) ? mdot : 0.0;
}


static double qspline(double x1, double y1, double dydx1, double x2, double y2, double dydx2, double x) {
    double c = 0.5*(dydx2-dydx1)/(x2-x1);
    double b = (y2-y1)/(x2-x1) - c*(x2+x1);
    double a = y2 - b*x2 - c*x2*x2;
    return a + b*x + c*x*x;
}



#ifdef DISK_SD_DEBUG
int main(int argc, char *argv[])
{
    double mdot = 0.5;
    double spin = 0.5;
    
    if (argc > 1) spin = atof(argv[1]);
    if (argc > 2) mdot = atof(argv[2]);
    
    printf("# model: %s\n", model_name);
    printf("# mdot = %.4f\n", mdot);
    printf("# spin = %.4f\n", spin);
    
    char options[256];
    sprintf(options, "path=./data,mdot=%f", mdot);
    
    int res = diskmodel_init(10.0, spin, options);
    if (res != 0) return res;

    diskmodel_dump();

    return 0;
}
#endif


