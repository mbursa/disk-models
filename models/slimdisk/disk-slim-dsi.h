#ifndef _DISK_SD_DSI_H_
#define _DISK_SD_DSI_H_

#define DSI_MASS            10.0        // canonical BH mass used in tables

#define DSI_COL_R           0           // radius
#define DSI_COL_VR          1           // radial velocity
#define DSI_COL_TC          2           // midplane temperature
#define DSI_COL_L           3           // angular momentum ((u_phi)
#define DSI_COL_SIGMA       4           // integrated (column) density
#define DSI_COL_HR          5           // H/R ratio
#define DSI_COL_P           6           // integrated pressure
#define DSI_COL_PRPG        7           // radiation to gass pressure ratio
#define DSI_COL_F           8           // flux (from both surfaces)
#define DSI_COL_ELL         10          // specific angular momentum (u_phi / u_t)


int dsi_init(char* path, double mdot, double spin, double alpha, double *rmin, double *rmax);

double dsi_eval(double R, int N);

void dsi_done();


#endif


