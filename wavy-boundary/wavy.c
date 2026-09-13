#include "grid/quadtree.h"
#include "EBM_VOF/myembed.h"
#include "navier-stokes/centered.h"
#include "EBM_VOF/embed_contact.h"
#include "EBM_VOF/embed_two-phase.h"
#include "EBM_VOF/embed_tension.h"
#include "view.h"
#include "output-mpi.h"

////// Mesh size variable
#define MESH_LEVEL	8

////// Dimensional quanities
#define rhoLeft		1000.		// Density of left  phase (kg/m^3)
#define muLeft		5.e-2		// Dynamic viscosity of left  phase (kg/m/s)

#define rhoRight	1.2			// Density of right phase (kg/m^3)
#define muRight		1.8e-5		// Dynamic viscosity of right phase (kg/m/s)

#define gamma		0.035		// Surface tension between phases (kg/s^2)

#define p_in		1.e5		// Atmospheric pressure at inlet (Pa)
#define p_out		1.e5		// Atmospheric pressure at outlet (Pa)

#define t_end		100.		// End time of simulation (s)

////// Dimensioned scales
#define L_ref		0.1			// Reference length of channel
#define R_ref		0.001		// Reference half-width of channel

#define U_ref		(gamma*R_ref/(muLeft*L_ref)) 	// Reference velocity scale (From Ca = R/L)
#define T_ref		(R_ref/U_ref) 					// Reference time scale
#define P_ref		(gamma/R_ref)					// Reference pressure scale
#define nd_P(p)		((p - p_out)/P_ref)				// Nondimensionalised pressure


////// Corrensponding dimensionless quantities
#define eps			(R_ref/L_ref)					// Ratio of length scales
#define epsInv		(L_ref/R_ref)					// Inverse ratio of length scales
#define Re			(rhoLeft*U_ref*R_ref/muLeft)	// Reynolds number
#define Ca			(muLeft*U_ref/gamma)			// Capillary number
#define rhoRatio	(rhoRight/rhoLeft)				// Ratio of densities
#define muRatio		(muRight/muLeft)				// Ratio of viscosities

double thetac		= 60.;							// Contact angle in degrees
#define thetar		(thetac*pi/180.)				// Contact angle in degrees

#define t1			(t/T_ref)						// Dimensionless time units
#define t1_end		(t_end/T_ref)					// Dimensionless end time
#define t_spacing	0.2

#define growth		(2*L_ref*(cos(thetar)/R_ref + (p_in-p_out)/gamma)/3)	// Dimensionless leading-order growth rate
#define p_jump		(cos(thetar))											// Dimensionless leading-order pressure jump over interface

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


p[left] = dirichlet(nd_P(p_in));
pf[left] = dirichlet_face(nd_P(p_in)); // Unsure if dirichlet_face is required here, as opposed to dirichlet

p[right] = dirichlet(nd_P(p_out));
pf[right] = dirichlet_face(nd_P(p_out));

f[embed] = neumann(0.);
f[bottom] = neumann(0.);
f[left] = dirichlet(1.);
f[right] = dirichlet(0.);

cs[left] = neumann(0.);
cs[right] = neumann(0.);



int main() {

	size(epsInv);
	origin(0., 0.);
	
	// Fluid properties (1: left, 2: right)
	rho1 = Re*Ca, rho2 = Re*Ca*rhoRatio;
	mu1 = Ca, mu2 = Ca*muRatio;
	f.sigma = 1.;

	// Attributes required for EBM
	tmp_c.height = tmp_h;
	tmp_c.hnew1 = hnew1;
	tmp_c.oxyi = o_interface;
	tmp_c.nc = ncc;

	fprintf(stdout, "\n\nMesh level: %2d\n", MESH_LEVEL);
	
	fprintf(stdout, "\nReynolds: %0.9g\n", Re);
	fprintf(stdout, "epsReynolds: %0.9g\n", eps*Re);
	fprintf(stdout, "Capillary: %0.9g\n", Ca);
	fprintf(stdout, "Density ratio: %0.9g\n", rhoRatio);
	fprintf(stdout, "Viscosity ratio: %0.9g\n", muRatio);

	fprintf(stdout, "\nVelocity scale: %0.9g\n", U_ref);
	fprintf(stdout, "Time scale: %0.9g\n", T_ref);
	fprintf(stdout, "Pressure scale: %0.9g\n", P_ref);
	
	fprintf(stdout, "\nDimensionless final time: %0.9g\n", t1_end);
	fprintf(stdout, "Leading-order growth rate: %0.9g\n", growth);
	fprintf(stdout, "Leading-order pressure jump: %0.9g\n", p_jump);

	fprintf(stdout, "\n%6s %8s %9s\n", "i", "t", "CPU time");
	init_grid(1 << MESH_LEVEL);
	run();
	
}





event init(t = 0) {
	solid(cs, fs, (1. + eps*(sin(8 * pi * eps * x) + eps)) - y);
	// solid(cs, fs, 1.e-6 + 1. - y);
	cleansmallcell(cs, fs, csTL);

	if(thetac == 90.) {
		fraction(f, 1 - x);
	} else {
		fraction(f, 1. + (1. - sqrt(1 - sq(y * cos(thetar)))) / cos(thetar) - x);
	}

	foreach() {
		f[] *= cs[];
		contact_angle[] = thetac;
	}
}

event early_end(i++, t1 <= t1_end) {
	double xx = epsInv - 0.1;
	double dyy = 0.1;
	for (double yy = dyy; yy <= 2; yy += dyy) {
		if (interpolate(cs, xx, yy) * interpolate(f, xx, yy) > 1e-2) {
			fprintf(stdout, "Fluid has reached end of channel at time t = %g\n", t1);
			return 1;
		}
	}
}

event log(i += 100) {
	
	timing s = timer_timing(perf.gt, i, perf.tnc, NULL);
	fprintf(stdout, "%06d %8.4g %9.7g\n", i, t1, s.real);
}

event images(t += 10.*t_spacing*T_ref) {
	scalar u_r[], l[];
	foreach() {
		u_r[] = sqrt(sq(u.x[]) + sq(u.y[]));
		l[] = level;
	}

	output_ppm(u_r, file = "movies/movie-u.mp4",   box={{0.,0.}, {epsInv,2}});
	output_ppm(u.x, file = "movies/movie-u_x.mp4", box={{0.,0.}, {epsInv,2}});
	output_ppm(u.y, file = "movies/movie-u_y.mp4", box={{0.,0.}, {epsInv,2}});
	output_ppm(f,   file = "movies/movie-f.mp4",   box={{0.,0.}, {epsInv,2}});
	output_ppm(p,   file = "movies/movie-p.mp4",   box={{0.,0.}, {epsInv,2}});
	output_ppm(l,   file = "movies/movie-l.mp4",   box={{0.,0.}, {epsInv,2}});
	
	// Nice video of the fluids, with pressure shown in the top half and horizontal velocity in the bottom half
	clear();
	view(tx=-0.5, ty = 0.,
		 fov = 1.,
		 width = 1080, height = 50);

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
	save("movies/flow-p-u_x.mp4");

}

event fields(t += t_spacing*T_ref) {

	// Name of directory for field outputs and other variables
	char* spec_dir = "interface";
	char* file_name;
	int file_name_len;

	file_name_len = snprintf(NULL, 0, "%s/%s-%09.3f.dat", spec_dir, spec_dir, t1)+1;
	file_name = malloc(file_name_len);
	snprintf(file_name, file_name_len, "%s/%s-%09.3f.dat", spec_dir, spec_dir, t1);

	output_facets_mpi(f, file_name);

}