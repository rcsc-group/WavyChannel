#include "grid/quadtree.h"
#include "EBM_VOF/myembed.h"
#include "navier-stokes/centered.h"
#include "EBM_VOF/embed_contact.h"
#include "EBM_VOF/embed_two-phase.h"
#include "EBM_VOF/embed_tension.h"
#include "view.h"
#include "output-mpi.h"

////// Mesh size variables
#define MAX_LEVEL	11
#define MIN_LEVEL	7

////// Dimensional quanities
// #define rhoLeft		1000.		// Density of left  phase (kg/m^3) (water)
// #define rhoLeft		1100.		// Density of left  phase (kg/m^3) (resin)
#define rhoLeft		870.		// Density of left  phase (kg/m^3) (oil)

// #define muLeft		1.e-3		// Dynamic viscosity of left  phase (kg/m/s) (water)
// #define muLeft		1.0			// Dynamic viscosity of left  phase (kg/m/s) (resin)
#define muLeft		1.e-2			// Dynamic viscosity of left  phase (kg/m/s) (oil)

// #define rhoRight		1.2			// Density of right phase (kg/m^3) (air)
#define rhoRight	1000		// Density of right phase (kg/m^3) (water)

// #define muRight		1.8e-5		// Dynamic viscosity of right phase (kg/m/s) (air)
#define muRight		1.e-3		// Dynamic viscosity of right phase (kg/m/s) (water)

// #define sig			0.072		// Surface tension between phases (kg/s^2) (water-air)
// #define sig			0.045		// Surface tension between phases (kg/s^2) (resin-air)
// #define sig			0.035		// Surface tension between phases (kg/s^2) (oil-air)
#define sig			0.035		// Surface tension between phases (kg/s^2) (oil-water)

#define p_in		1.e5		// Atmospheric pressure at inlet
#define p_out		1.e5		// Atmospheric pressure at outlet


////// Dimensioned scales
#define L_ref		0.1			// Reference length of channel
#define R_ref		0.001		// Reference half-width of channel

#define U_ref		(sig*R_ref/(muLeft*L_ref)) 		// Reference velocity scale (From Ca = R/L)
#define T_ref		(R_ref/U_ref) 					// Reference time scale
#define P_ref		(sig/R_ref)						// Reference pressure scale


////// Corrensponding dimensionless quantities
#define eps			(R_ref/L_ref)					// Ratio of length scales
#define epsInv		(L_ref/R_ref)					// Inverse ratio of length scales
#define Re			(rhoLeft*U_ref*R_ref/muLeft)	// Reynolds number
#define Ca			(muLeft*U_ref/sig)				// Capillary number
#define rhoRatio	(rhoRight/rhoLeft)				// Ratio of densities
#define muRatio		(muRight/muLeft)				// Ratio of viscosities

double thetac		= 60.;							// Contact angle in degrees
#define thetar		(thetac*pi/180.)				// Contact angle in degrees


// Simulation time
double t_end = 100.0;

vector tmp_h[], o_interface[], ncc[], hnew1[];
double csTL = max(1.e-2, VFTL);

// Boundary conditions
u.t[embed] = dirichlet(0.); 
u.n[embed] = dirichlet(0.);

u.t[left] = neumann(0.); // Free influx
u.n[left] = neumann(0.);

u.t[right] = neumann(0.); // Free outflux
u.n[right] = neumann(0.);

u.t[bottom] = neumann(0.); // Symmetry
u.n[bottom] = dirichlet(0.);


p[left] = dirichlet(p_in/P_ref);
pf[left] = dirichlet_face(p_in/P_ref); // Unsure if dirichlet_face is required here, as opposed to dirichlet

p[right] = dirichlet(p_out/P_ref);
pf[right] = dirichlet_face(p_out/P_ref);

f[embed] = neumann(0.);
f[bottom] = neumann(0.);
f[left] = dirichlet(1.);
f[right] = dirichlet(0.);

cs[left] = neumann(0.);
cs[right] = neumann(0.);



int main() {

	size(epsInv);
	origin(0., 0.);
	
	// 1: left
	// 2: right
	rho1 = Re*Ca, rho2 = Re*Ca*rhoRatio;
	mu1 = Ca, mu2 = Ca*muRatio;

	f.sigma = 1.;

	tmp_c.height = tmp_h;
	tmp_c.hnew1 = hnew1;
	tmp_c.oxyi = o_interface;
	tmp_c.nc = ncc;

	// Parameters for Poisson solver
	// NITERMIN = 4;
	// TOLERANCE = 1.e-4;

	fprintf(stdout, "\n\nMax level: %2d\nMin level: %2d\n", MAX_LEVEL, MIN_LEVEL);
	fprintf(stdout, "\nReynolds: %0.9g\n", Re);
	fprintf(stdout, "epsReynolds: %0.9g\n", eps*Re);
	fprintf(stdout, "Capillary: %0.9g\n", Ca);
	fprintf(stdout, "Density ratio: %0.9g\n", rhoRatio);
	fprintf(stdout, "Viscosity ratio: %0.9g\n", muRatio);
	fprintf(stdout, "Velocity scale: %0.9g\n", U_ref);
	fprintf(stdout, "Time scale: %0.9g\n", T_ref);
	fprintf(stdout, "Pressure scale: %0.9g\n", P_ref);

	fprintf(stderr, "\n%6s %8s %9s ;; %9s %9s ;; %5s %5s %7s\n", "i", "t", "L^2 norm", "Wall time", "CPU time", "mgu.i", "mgp.i", "N");
	init_grid(1 << MAX_LEVEL);
	run();
	
}





event init(t=0) {
	solid(cs, fs, (1. + eps*sin(8 * pi * eps * x)) - y);
	// solid(cs, fs, 1.e-6 + 1. - y);
	cleansmallcell(cs, fs, csTL);

	if(thetac == 90.) {
		fraction(f, 1 - x);
	} else {
		fraction(f, 0.5*epsInv + (1 - sqrt(1 - sq(y * cos(thetar)))) / cos(thetar) - x);
	}

	foreach() {
		f[] *= cs[];
		contact_angle[] = thetac;
	}
}

event early_end(i++, t <= t_end) {
	double xx = epsInv - 0.1;
	double dyy = 0.1;
	for (double yy = dyy; yy <= 2; yy += dyy) {
		if (interpolate(cs, xx, yy) * interpolate(f, xx, yy) > 1e-2) {
			fprintf(stderr, "Fluid has reached end of channel at time t = %g\n", t);
			return 1;
		}
	}
}

event logfile(i += 100) {
	// fprintf(stderr, "%d %g %d %d %g\n", i, t, mgp.i, mgu.i, du);
	int n_cells = 0;
	foreach(reduction(+:n_cells)){
		n_cells += 1;
	}
	timing s = timer_timing(perf.gt, i, perf.tnc, NULL);
	fprintf(stderr, "%06d %8.4g %9.5g ;; %9.7g %9.7g ;; %5d %5d %7d\n", i, t, normf(u.x).rms, s.real, s.cpu, mgu.i, mgp.i, n_cells);
}

event images(t += 0.01) {
	scalar u_r[], l[];
	foreach() {
		u_r[] = sqrt(sq(u.x[]) + sq(u.y[]));
		l[] = level;
	}

	// output_ppm(u_r, file = "movie-u.mp4",   box={{0.,0.}, {length,depth*2}}, linear=true);
	// output_ppm(u.x, file = "movie-u_x.mp4", box={{0.,0.}, {length,depth*2}}, linear=true);
	// output_ppm(u.y, file = "movie-u_y.mp4", box={{0.,0.}, {length,depth*2}}, linear=true);
	// output_ppm(f,   file = "movie-f.mp4",   box={{0.,0.}, {length,depth*2}}, linear=true);
	// output_ppm(p,   file = "movie-p.mp4",   box={{0.,0.}, {length,depth*2}}, linear=true);
	// output_ppm(l,   file = "movie-l.mp4",   box={{0.,0.}, {length,depth*2}}, linear=true);

	output_ppm(u_r, file = "movie-u.mp4",   box={{0.,0.}, {epsInv,2}});
	output_ppm(u.x, file = "movie-u_x.mp4", box={{0.,0.}, {epsInv,2}});
	output_ppm(u.y, file = "movie-u_y.mp4", box={{0.,0.}, {epsInv,2}});
	output_ppm(f,   file = "movie-f.mp4",   box={{0.,0.}, {epsInv,2}});
	output_ppm(p,   file = "movie-p.mp4",   box={{0.,0.}, {epsInv,2}});
	output_ppm(l,   file = "movie-l.mp4",   box={{0.,0.}, {epsInv,2}});
	
	// Nice video of the fluid fraction
	clear();
	view(tx=-0.5, ty = 0.,
		 fov = 6,
		 width = 1080, height = 216);

	double grey = 0.8;

	draw_vof("cs", "fs", filled=-1, fc={grey, grey, grey});
	// draw_vof("f", filled=1, fc={0.5, 0.5, 1.}, lw=2, min=0., max=1.);
	squares("p", linear=true);
	// draw_vof("f", lw=3);
	isoline("f", 0.5);
	// cells();
	mirror (n = {0,1}) {
		draw_vof("cs", "fs", filled=-1, fc={grey, grey, grey});
		squares("u.x", linear=true);
		// draw_vof("f", lw=3);
		isoline("f", 0.5);
		// cells();
	}
	save("flow.mp4");

}

#if TREE
event adapt (i++) {
  scalar sf1[];
  foreach() {
    sf1[] = (8. * tmp_c[] +
	     4. * (tmp_c[-1] + tmp_c[1] +
		   tmp_c[0, 1] + tmp_c[0, -1] +
		   tmp_c[0, 0, 1] + tmp_c[0, 0, -1]) +
	     2. * (tmp_c[-1, 1] + tmp_c[-1, 0, 1] + tmp_c[-1, 0, -1] + tmp_c[-1, -1] +
		   tmp_c[0, 1, 1] + tmp_c[0, 1, -1] + tmp_c[0, -1, 1] + tmp_c[0, -1, -1] +
		   tmp_c[1, 1] + tmp_c[1, 0, 1] + tmp_c[1, -1] + tmp_c[1, 0, -1]) +
	     tmp_c[1, -1, 1] + tmp_c[-1, 1, 1] + tmp_c[-1, 1, -1] + tmp_c[1, 1, 1] +
	     tmp_c[1, 1, -1] + tmp_c[-1, -1, -1] + tmp_c[1, -1, -1] + tmp_c[-1, -1, 1]) / 64.;
    sf1[] += cs[];
  }
  adapt_wavelet ({sf1}, (double[]){1e-5}, minlevel = MIN_LEVEL, maxlevel = MAX_LEVEL);
}
#endif

event fields(t = 0.; t += 0.01) {

	// Name of directory for field outputs and other variables
	char* out_dir = "fields";
	char* spec_dir;
	char* file_name;
	int file_name_len;
	FILE* fp;

	// Interface positions
	spec_dir = "interface";

	file_name_len = snprintf(NULL, 0, "%s/%s/%s-%09.6f.dat", out_dir, spec_dir, spec_dir, t)+1;
	file_name = malloc(file_name_len);
	snprintf(file_name, file_name_len, "%s/%s/%s-%09.6f.dat", out_dir, spec_dir, spec_dir, t);

	// fp = fopen(file_name, "w");
	// output_facets(f, fp);
	// fclose(fp);
	output_facets_mpi(f, file_name);

	// Pressure field
	spec_dir = "pressure";

	file_name_len = snprintf(NULL, 0, "%s/%s/%s-%09.6f.dat", out_dir, spec_dir, spec_dir, t)+1;
	file_name = malloc(file_name_len);
	snprintf(file_name, file_name_len, "%s/%s/%s-%09.6f.dat", out_dir, spec_dir, spec_dir, t);

	fp = fopen(file_name, "w");

	int data_points = 2000;
	double dxx = epsInv / (double)(data_points+2);
	for(int nx = 1; nx <= data_points; nx++) {

		// Basically extracting the pressure from the centreline
		double xx = (double)nx * dxx;
		double yy = 0.01;
		fprintf(fp, "%f %f\n", xx, interpolate(p, xx, yy));

	}
	fclose(fp);

}

event status(t += 0.1) {

	// Name of directory for field outputs and other variables
	char* out_dir = "fields";
	char* spec_dir;
	char* file_name;
	int file_name_len;
	FILE* fp;

	// Phase field
	spec_dir = "phase";

	file_name_len = snprintf(NULL, 0, "%s/%s/%s-%09.6f.png", out_dir, spec_dir, spec_dir, t)+1;
	file_name = malloc(file_name_len);
	snprintf(file_name, file_name_len, "%s/%s/%s-%09.6f.png", out_dir, spec_dir, spec_dir, t);

	fp = fopen(file_name, "w");
	output_ppm(f, fp = fp, spread=-1, box={{0.,0.}, {epsInv,2}});
	fclose(fp);

	scalar l[];
	foreach() {
		l[] = level;
	}

	// Mesh level
	spec_dir = "level";

	file_name_len = snprintf(NULL, 0, "%s/%s/%s-%09.6f.png", out_dir, spec_dir, spec_dir, t)+1;
	file_name = malloc(file_name_len);
	snprintf(file_name, file_name_len, "%s/%s/%s-%09.6f.png", out_dir, spec_dir, spec_dir, t);

	fp = fopen(file_name, "w");
	output_ppm(l, fp = fp, spread=-1, box={{0.,0.}, {epsInv,2}});
	fclose(fp);

}