#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sim5lib.h"
#include "disk-slim-dsi.h"



// debug print macro
#ifdef DISK_SD_DEBUG
#define debug(...) fprintf(stderr, __VA_ARGS__)
#else
#define debug(...) {}
#endif



#define DSI_EPS             1e-9        // small factor
#define DSI_COLS            12          // number of columns in solution files
#define DSI_NRAD            1000        // number of radial grid points
#define DSI_ALPHAS          6           // number of alpha values



struct solution {
    double alpha;
    double spin;
    double mdot;
    int l, m, n;
};


// static contained for the interpolated tables
static sim5interp dsi_interp[DSI_COLS] = {0};


// array of alpha values contained in slim disk tabular data
static double alphas[DSI_ALPHAS] = {0.005, 0.010, 0.025, 0.050, 0.075, 0.100};


// count lines in file (a privately used function)
static long dsi_count_lines(FILE* f) {
    long nlines = 0;
    char ch;
    fseek(f, 0, SEEK_SET);
    while(!feof(f)) if ((ch = fgetc(f)) == '\n') nlines++;
    fseek(f, 0, SEEK_SET);
    return nlines;
}


int dsi_init(char* path, double mdot, double spin, double alpha, double *rmin, double *rmax)
// Initialization routine for DSI slim disk data interpolator.
//
// The function reads the tabulated slim disk solutions data (http://users.camk.edu.pl/as/ss.tar.gz), 
// finds solutions closest to the required parameters and sets up interpolation tables.
// 
// On output, it returns status code and rmin/rmax values. 
// Status code is 0 in successful initialization and less than zero otherwise.
// rmin/rmax tells the range of radii at which the set of solutions that has been used to interpolate 
// requested model provides some solid data. Precisely, rmin is the smallest (not largest!!) common 
// radius of the selected set of tabular solutions (meaning some bit of extrapolation may be required 
// to obtain solutions) and rmax is the smallest common outer radius if the set not greated than 1000 rg.
// Note that the interpolator allows for extrapolation  (linear in log-log space), so that data can 
// be obtained also for radii that are smaller that rmin or larger than rmax with all the downsides that 
// such an extrapolation may have. It also has to be noted that usually all quantities have smooth 
// profiles at inner and (especially) outer boundary, so that in most cases this extrapolation 
// work well. Some exceptions include VR (radial velocity) extrapolation towards low radii where the 
// result may exceed one (=speed of light).
// 
// Params:
//      path: path to the directory with the solutions
//      mdot: mass accretion rate in Mdot_Edd (2.255e18 g/s M/Msun)
//      spin: BH spin parameter
//      alpha: alpha-viscosity parameter
//      rmin: inner radius for which selected model provides data 
//            (the smallest inner radius that is covered by tabular data)
//      rmin: outer radius for which selected model provides data 
//            (the lergest outer radius that is covered by tabular data)
// Returns:
//      0 on sucessful intialization, negative code on error
//      Inner/outer radius covered by data returned in rmin/rmax. 
{
    debug("====== dsi_init =======\n");
    debug("path: %s\n", path);
    debug("mdot: %.5f\n", mdot);
    debug("spin: %.5f\n", spin);
    debug("vics: %.5f\n", alpha);


    long i,j;
    int ialpha1, ialpha2;
    struct solution solindex[2][2][2];  // order: alpha, spin, mdot

    // determine alpha bracket
    ialpha1 = ialpha2 = 0;
    for (i=0; i<DSI_ALPHAS-1; i++) {
        if ((alphas[i]-DSI_EPS <= alpha) && (alphas[i+1]+DSI_EPS > alpha)) {
            ialpha1 = i;
            ialpha2 = i+1;
            break;
        }
    }

    // return error code if no alpha bracket can be found    
    if (ialpha1 == ialpha2) {
        fprintf(stderr, "dsi_setup: cannot bracket alpha (%.4e; min=%.4e; max=%.4e)\n", alpha, alphas[0], alphas[DSI_ALPHAS-1]);
        return -1;
    }

    debug("found alpha bracket: %d  %d  (%.5e / %.5e)\n", ialpha1, ialpha2, alphas[ialpha1], alphas[ialpha2]);


    // open index files for alpha1 and alpha2
    FILE* findex1;
    FILE* findex2;
    char indexpath1[256];
    char indexpath2[256];
    sprintf(indexpath1, "%s/alpha-%.3f/res.mamdot.dat", path, alphas[ialpha1]);
    sprintf(indexpath2, "%s/alpha-%.3f/res.mamdot.dat", path, alphas[ialpha2]);
    debug("opening index file1: %s\n", indexpath1);
    debug("opening index file2: %s\n", indexpath2);
    findex1 = fopen(indexpath1, "r");
    findex2 = fopen(indexpath2, "r");

    if (!findex1) {
        fprintf(stderr, "dsi_setup: cannot open index file %s\n", indexpath1);
        return -1;
    }

    if (!findex2) {
        fprintf(stderr, "dsi_setup: cannot open index file %s\n", indexpath2);
        return -1;
    }

    // initialize solindex
    solindex[0][0][0] = solindex[0][0][1] = (struct solution){0};
    solindex[0][1][0] = solindex[0][1][1] = (struct solution){0};
    solindex[1][0][0] = solindex[1][0][1] = (struct solution){0};
    solindex[1][1][0] = solindex[1][1][1] = (struct solution){0};

    // alpha1 bracket
    debug("== alpha1 bracket (a=%e mdot=%e)\n", spin, mdot);
    fseek(findex1, 0, SEEK_SET);
    for(;;) {
        double tmp_mdot1, tmp_mdot2, tmp_spin1, tmp_spin2;
        int tmp_a1, tmp_a2, tmp_m1, tmp_m2;

        if (fscanf(findex1, "%*f %*f %lf %lf %*f %*f %d %d\n", &tmp_mdot1, &tmp_spin1, &tmp_a1, &tmp_m1) == EOF) break;
        int tmppos = ftell(findex1);
        if (fscanf(findex1, "%*f %*f %lf %lf %*f %*f %d %d\n", &tmp_mdot2, &tmp_spin2, &tmp_a2, &tmp_m2) == EOF) break;
        
        if (tmp_a1 != tmp_a2) {
            fseek(findex1, tmppos, SEEK_SET);
            continue;
        }

        //debug("file1 scan: a1=%.4e mdot1=%.4e n1=%2d m1=%2d \n", tmp_spin1, tmp_mdot1, tmp_a1, tmp_m1);
        //debug("file1 scan: a2=%.4e mdot2=%.4e n2=%2d m2=%2d \n", tmp_spin2, tmp_mdot2, tmp_a2, tmp_m2);

        if (tmp_spin1+(spin>DSI_EPS?DSI_EPS:-DSI_EPS) <= spin) {
            if ((tmp_mdot1 <= mdot) && (mdot <= tmp_mdot2)) {
                debug("hit-lo: mdot1=%.4e mdot2=%.4e mdot=%.4e \n", tmp_mdot1, tmp_mdot2, mdot);
                solindex[0][0][0] = (struct solution){alphas[ialpha1], tmp_spin1, tmp_mdot1, ialpha1, tmp_a1, tmp_m1};
                solindex[0][0][1] = (struct solution){alphas[ialpha1], tmp_spin1, tmp_mdot2, ialpha1, tmp_a1, tmp_m2};
            }
        }

        if ((tmp_spin1 > solindex[0][0][0].spin) && (tmp_spin1+DSI_EPS >= spin)) {
            if ((tmp_mdot1 <= mdot) && (mdot <= tmp_mdot2)) {
                debug("hit-hi: mdot1=%.4e mdot2=%.4e mdot=%.4e \n", tmp_mdot1, tmp_mdot2, mdot);
                solindex[0][1][0] = (struct solution){alphas[ialpha1], tmp_spin1, tmp_mdot1, ialpha1, tmp_a1, tmp_m1};
                solindex[0][1][1] = (struct solution){alphas[ialpha1], tmp_spin1, tmp_mdot2, ialpha1, tmp_a1, tmp_m2};
                break;
            }
        }

        fseek(findex1, tmppos, SEEK_SET);
    }


    // alpha2 bracket
    debug("== alpha2 bracket (a=%e mdot=%e)\n", spin, mdot);
    fseek(findex2, 0, SEEK_SET);
    for(;;) {
        double tmp_mdot1, tmp_mdot2, tmp_spin1, tmp_spin2;
        int tmp_a1, tmp_a2, tmp_m1, tmp_m2;

        if (fscanf(findex2, "%*f %*f %lf %lf %*f %*f %d %d\n", &tmp_mdot1, &tmp_spin1, &tmp_a1, &tmp_m1) == EOF) break;
        int tmppos = ftell(findex2);
        if (fscanf(findex2, "%*f %*f %lf %lf %*f %*f %d %d\n", &tmp_mdot2, &tmp_spin2, &tmp_a2, &tmp_m2) == EOF) break;
        
        if (tmp_a1 != tmp_a2) {
            fseek(findex2, tmppos, SEEK_SET);
            continue;
        }

        //debug("file2 scan: a1=%.4e mdot1=%.4e n1=%2d m1=%2d \n", tmp_spin1, tmp_mdot1, tmp_a1, tmp_m1);
        //debug("file2 scan: a2=%.4e mdot2=%.4e n2=%2d m2=%2d \n", tmp_spin2, tmp_mdot2, tmp_a2, tmp_m2);

        if (tmp_spin1+(spin>DSI_EPS?DSI_EPS:-DSI_EPS) <= spin) {
            if ((tmp_mdot1 <= mdot) && (mdot <= tmp_mdot2)) {
                debug("hit-lo: mdot1=%.4e mdot2=%.4e mdot=%.4e \n", tmp_mdot1, tmp_mdot2, mdot);
                solindex[1][0][0] = (struct solution){alphas[ialpha2], tmp_spin1, tmp_mdot1, ialpha2, tmp_a1, tmp_m1};
                solindex[1][0][1] = (struct solution){alphas[ialpha2], tmp_spin1, tmp_mdot2, ialpha2, tmp_a1, tmp_m2};
            }
        }

        if ((tmp_spin1 > solindex[0][0][0].spin) && (tmp_spin1+DSI_EPS >= spin)) {
            if ((tmp_mdot1 <= mdot) && (mdot <= tmp_mdot2)) {
                debug("hit-hi: mdot1=%.4e mdot2=%.4e mdot=%.4e \n", tmp_mdot1, tmp_mdot2, mdot);
                solindex[1][1][0] = (struct solution){alphas[ialpha2], tmp_spin1, tmp_mdot1, ialpha2, tmp_a1, tmp_m1};
                solindex[1][1][1] = (struct solution){alphas[ialpha2], tmp_spin1, tmp_mdot2, ialpha2, tmp_a1, tmp_m2};
                break;
            }
        }

        fseek(findex2, tmppos, SEEK_SET);
    }

    // close index files
    fclose(findex1);
    fclose(findex2);

    debug("file1 bracket[0][0][0]: spin=%.4e  mdot=%.4e  m=%d  n=%d\n", solindex[0][0][0].spin, solindex[0][0][0].mdot, solindex[0][0][0].m, solindex[0][0][0].n);
    debug("file1 bracket[0][0][1]: spin=%.4e  mdot=%.4e  m=%d  n=%d\n", solindex[0][0][1].spin, solindex[0][0][1].mdot, solindex[0][0][1].m, solindex[0][0][1].n);
    debug("file1 bracket[0][1][0]: spin=%.4e  mdot=%.4e  m=%d  n=%d\n", solindex[0][1][0].spin, solindex[0][1][0].mdot, solindex[0][1][0].m, solindex[0][1][0].n);
    debug("file1 bracket[0][1][1]: spin=%.4e  mdot=%.4e  m=%d  n=%d\n", solindex[0][1][1].spin, solindex[0][1][1].mdot, solindex[0][1][1].m, solindex[0][1][1].n);

    debug("file2 bracket[1][0][0]: spin=%.4e  mdot=%.4e  m=%d  n=%d\n", solindex[1][0][0].spin, solindex[1][0][0].mdot, solindex[1][0][0].m, solindex[1][0][0].n);
    debug("file2 bracket[1][0][1]: spin=%.4e  mdot=%.4e  m=%d  n=%d\n", solindex[1][0][1].spin, solindex[1][0][1].mdot, solindex[1][0][1].m, solindex[1][0][1].n);
    debug("file2 bracket[1][1][0]: spin=%.4e  mdot=%.4e  m=%d  n=%d\n", solindex[1][1][0].spin, solindex[1][1][0].mdot, solindex[1][1][0].m, solindex[1][1][0].n);
    debug("file2 bracket[1][1][1]: spin=%.4e  mdot=%.4e  m=%d  n=%d\n", solindex[1][1][1].spin, solindex[1][1][1].mdot, solindex[1][1][1].m, solindex[1][1][1].n);

    // tests
    if (
        ((solindex[0][0][0].m == solindex[0][0][1].m) && (solindex[0][0][0].n == solindex[0][0][1].n)) ||
        ((solindex[0][1][0].m == solindex[0][1][1].m) && (solindex[0][1][0].n == solindex[0][1][1].n)) ||
        ((solindex[1][0][0].m == solindex[1][0][1].m) && (solindex[1][0][0].n == solindex[1][0][1].n)) ||
        ((solindex[1][1][0].m == solindex[1][1][1].m) && (solindex[1][1][0].n == solindex[1][1][1].n))
    ) {
        fprintf(stderr, "dsi_init: the requested solution lies outside of the grid\n");
        return -1;
    }



    // allocate dynamic array of interpolator objects with dimensions [DSI_COL][2][2][2]
    sim5interp interp[DSI_COLS][2][2][2];

    int iv, ia, im;
    *rmin = 1e9;
    *rmax = 1e3;

    for (iv=0; iv<2; iv++) {
        for (ia=0; ia<2; ia++) {
            for (im=0; im<2; im++) {
                char fn[256];
                FILE* fsol;
                sprintf(fn, "%s/alpha-%.3f/soltt.%d.%d.dat", path, solindex[iv][ia][im].alpha, solindex[iv][ia][im].m, solindex[iv][ia][im].n);
                debug(">>> openning %s \n", fn);
                fsol = fopen(fn, "r");  

                // count lines in solution file and allocate space for values arrays
                int nlines = dsi_count_lines(fsol);
                double* values[DSI_COLS];
                for (i=0; i<DSI_COLS; i++) values[i] = (double*) malloc(nlines*sizeof(double));
                debug("has %d lines\n", nlines);

                // read values 
                int line = 0;
                for(;;) {
                    int scanr = fscanf(
                        fsol, 
                        "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
                        &values[0][line], &values[1][line], &values[2][line], &values[3][line], &values[4][line], &values[5][line], 
                        &values[6][line], &values[7][line], &values[8][line], &values[9][line], &values[10][line], &values[11][line]
                    );
                    if (scanr == EOF) break;
                    line++;
                }
                
                // if radii go from high to low, we need to reverse order for all arrays
                if (values[0][0] > values[0][nlines-1]) {
                    debug("reversing arrays\n");
                    for (i=0; i<DSI_COLS; i++) {
                        int j1 = 0;             // first element
                        int j2 = nlines-1;      // last element
                        while (j1 < j2) {
                            double _temp = values[i][j1];
                            values[i][j1] = values[i][j2];
                            values[i][j2] = _temp;
                            j1++;
                            j2--;
                        }
                    }
                }

                // do some value adjustments
                for (i=0; i<nlines; i++) {
                    // convert radius to [rg] units
                    values[DSI_COL_R][i] = 2.*values[DSI_COL_R][i];

                    // make vr positive
                    values[DSI_COL_VR][i] = fabs(values[DSI_COL_VR][i]);
                }
                debug("rmin = %.3f  rmax = %.3f [rg]\n", values[0][0], values[0][nlines-1]);

                *rmin = fmin(*rmin, values[0][0]);
                *rmax = fmin(*rmax, values[0][nlines-1]-DSI_EPS);

                // setup interpolation objects for each quantity
                debug("building interpolators\n");
                for (i=1; i<DSI_COLS; i++) {
                    sim5_interp_init(&interp[i][iv][ia][im], values[0], values[i], nlines, INTERP_DATA_COPY, INTERP_TYPE_LOGLOG, INTERP_OPT_CAN_EXTRAPOLATE);
                }

                // free values arrays
                for (i=0; i<DSI_COLS; i++) free(values[i]);

                fclose(fsol);
            } // end of mdot cycle
        } // end of spin cycle
    } // end of alpha cycle

    debug("global rmin = %.4e\n", *rmin);
    debug("global rmax = %.4e\n", *rmax);

    // define temp interpolator    
    sim5interp tmp_interp;

    // setup target radial grid
    double* R = (double*) malloc(DSI_NRAD * sizeof(double));
    for (i=0; i<DSI_NRAD; i++) R[i] = exp(log(*rmin) + (log(*rmax+DSI_EPS)-log(*rmin))*(double)(i)/(double)(DSI_NRAD-1));
    
    // contract interpolators in mdot
    debug("contract interpolators in mdot\n");
    for (i=1; i<DSI_COLS; i++) {
        for (iv=0; iv<2; iv++) {
            for (ia=0; ia<2; ia++) {
                // calculate weight
                double w = (mdot - solindex[iv][ia][0].mdot)/(solindex[iv][ia][1].mdot - solindex[iv][ia][0].mdot);
                if (i==1) debug("v=%.4e  v1=%.4e  v2=%.4e  w=%.4f\n", mdot, solindex[iv][ia][0].mdot, solindex[iv][ia][1].mdot, w);
                // setup empty temp interpolator
                sim5_interp_init(&tmp_interp, NULL, NULL, DSI_NRAD, INTERP_DATA_BUILD, INTERP_TYPE_LOGLOG, INTERP_OPT_CAN_EXTRAPOLATE+INTERP_OPT_ACCEL);
                // fill the mdot-interpolated data
                for (j=0; j<DSI_NRAD; j++) {
                    double val = (1.-w)*sim5_interp_eval(&interp[i][iv][ia][0], R[j]) + w*sim5_interp_eval(&interp[i][iv][ia][1], R[j]);
                    sim5_interp_data_push(&tmp_interp, R[j], val);
                }
                // free the two original interpolators
                sim5_interp_done(&interp[i][iv][ia][0]);
                sim5_interp_done(&interp[i][iv][ia][1]);
                // and copy temp interpolator in place of the first original one
                interp[i][iv][ia][0] = tmp_interp;
            }
        }
    }

    // contract interpolators in spin
    debug("contract interpolators in spin\n");
    for (i=1; i<DSI_COLS; i++) {
        for (iv=0; iv<2; iv++) {
            // calculate weight
            double pw = 0.5;
            double w = (pow(r_ms(spin),pw) - pow(r_ms(solindex[iv][0][0].spin),pw))/(pow(r_ms(solindex[iv][1][0].spin),pw) - pow(r_ms(solindex[iv][0][0].spin),pw));
            if (i==2) debug("v=%.4e  v1=%.4e  v2=%.4e  w=%.4f\n", r_ms(spin), r_ms(solindex[iv][0][0].spin), r_ms(solindex[iv][1][0].spin), w);
            // setup empty temp interpolator
            sim5_interp_init(&tmp_interp, NULL, NULL, DSI_NRAD, INTERP_DATA_BUILD, INTERP_TYPE_LOGLOG, INTERP_OPT_CAN_EXTRAPOLATE+INTERP_OPT_ACCEL);
            // fill the mdot-interpolated data
            for (j=0; j<DSI_NRAD; j++) {
                // calculate weight in terms of rg instear of spin
                double val = (1.-w)*sim5_interp_eval(&interp[i][iv][0][0], R[j]) + w*sim5_interp_eval(&interp[i][iv][1][0], R[j]);
                sim5_interp_data_push(&tmp_interp, R[j], val);
            }
            // free the two mdot interpolators
            sim5_interp_done(&interp[i][iv][0][0]);
            sim5_interp_done(&interp[i][iv][1][0]);
            // and copy temp interpolator in place of the first original one
            interp[i][iv][0][0] = tmp_interp;
        }
    }

    // contract interpolators in alpha
    debug("contract interpolators in alpha\n");
    for (i=1; i<DSI_COLS; i++) {
        // calculate weight
        double w = (alpha - solindex[0][0][0].alpha)/(solindex[1][0][0].alpha - solindex[0][0][0].alpha);
        if (i==1) debug("v=%.4e  v1=%.4e  v2=%.4e  w=%.4f\n", alpha, solindex[0][0][0].alpha, solindex[1][0][0].alpha, w);

        // setup empty temp interpolator
        sim5_interp_init(&tmp_interp, NULL, NULL, DSI_NRAD, INTERP_DATA_BUILD, INTERP_TYPE_LOGLOG, INTERP_OPT_CAN_EXTRAPOLATE+INTERP_OPT_ACCEL);
        // fill the mdot-interpolated data
        for (j=0; j<DSI_NRAD; j++) {
            double val = (1.-w)*sim5_interp_eval(&interp[i][0][0][0], R[j]) + w*sim5_interp_eval(&interp[i][1][0][0], R[j]);
            if (i == DSI_COL_VR) val = fmin(val, 0.9999);
            sim5_interp_data_push(&tmp_interp, R[j], val);
        }
        // free the two mdot interpolators
        sim5_interp_done(&interp[i][0][0][0]);
        sim5_interp_done(&interp[i][1][1][0]);
        // and copy temp interpolator in place of the first original one
        interp[i][0][0][0] = tmp_interp;
    }
    
    // copy the contracted interpolators to static container
    debug("copy the contracted interpolators to static container\n");
    for (i=1; i<DSI_COLS; i++) dsi_interp[i] = interp[i][0][0][0];
    
    // free stuff
    free(R);
    
    return 0;
}




double dsi_eval(double R, int N)
// Evaluation of slim disk solutions.
//
// Makes a lookup in the interpolated table to find the value of the requested qunatity at given radius.
// It allows for extrapolation and does not complain if R is outside rmin/rmax range.
//
// Params:
//      R: equatorial radius (r at theta=pi/2)
//      N: quantity to obtain (see DSI_COL_xxx constants in header file)
//
//  Returns:
//      Value od the requested quantity at radius R.
{
    if ((N <= 0) || (N > 11)) {
        fprintf(stderr, "dsi_eval: invalid parameter index (%d)\n", N);
        return NAN;
    }

    return sim5_interp_eval(&dsi_interp[N], R);
}




void dsi_done()
// Frees memory used by interpolated tables.
{
    int i;
    for (i=1; i<DSI_COLS; i++) sim5_interp_done(&dsi_interp[i]);
}





/*
#ifdef DISK_SD_DEBUG
int main(int argc, char **argv)
{
    double mdot  = 0.1;
    double spin = 0.9;
    double alpha = 0.1;
    double rmin, rmax;
    double r;

    dsi_init("./data", mdot, spin, alpha, &rmin, &rmax);

    for (r = rmin; r<rmax; r*=1.01) printf(
        "%e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e\n", r, 
        dsi_eval(r,1), dsi_eval(r,2), dsi_eval(r,3), dsi_eval(r,4), dsi_eval(r,5), dsi_eval(r,6), dsi_eval(r,7), dsi_eval(r,8), dsi_eval(r,9), dsi_eval(r,10), dsi_eval(r,11)
    );

    dsi_done();
  
    return 0;
}
#endif
*/



