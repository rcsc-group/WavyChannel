#include "grid/quadtree.h"
#include "EBM_VOF/myembed.h"
#include "navier-stokes/centered.h"
#include "EBM_VOF/embed_contact.h"
#include "EBM_VOF/embed_two-phase.h"
#include "EBM_VOF/embed_tension.h"
#include "view.h"


// Mesh size variables
int level0 = 7, level_max = 10;

// Domain size variables
double length = 1., epsilon = 0.01;

// Simulation time
double t_end = 3.0;

vector tmp_h[], o_interface[], ncc[], hnew1[];
double thetac = 60.;
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


p[left] = dirichlet(0.);
pf[left] = dirichlet_face(0.); // Unsure if dirichlet_face is required here, as opposed to dirichlet

p[right] = dirichlet(0.);
pf[right] = dirichlet_face(0.);

f[embed] = neumann(0.);
f[bottom] = neumann(0.);
f[left] = dirichlet(1.);
f[right] = dirichlet(0.);

cs[left] = neumann(0.);
cs[right] = neumann(0.);


int main() {

	size(length);
	origin(0., 0.);
	
	// 1: left
	// 2: right
	rho1 = 1., rho2 = 1.e-3;
	mu1 = 0.1, mu2 = 0.0;

	f.sigma = 0.01;

	tmp_c.height = tmp_h;
	tmp_c.hnew1 = hnew1;
	tmp_c.oxyi = o_interface;
	tmp_c.nc = ncc;

	fprintf(stderr, "\n\nLevel: %d\n%6s %8s ;; %9s %9s ;; %5s %5s %7s\n", level_max, "i", "t", "L^2 norm", "CPU time", "mgu.i", "mgp.i", "N");
	init_grid(1 << level_max);
	run();
	
}





event init(t=0) {
	double depth = length * epsilon;
	solid(cs, fs, depth * (1. + 2*epsilon*sin(4 * pi * x / length)) - y);
	// solid(cs, fs, 1.e-6 + 1. - y);
	cleansmallcell(cs, fs, csTL);

	fraction(f, depth - x);
	foreach() {
		f[] *= cs[];
		contact_angle[] = thetac;
	}
}

event logfile(i += 100) {
	// fprintf(stderr, "%d %g %d %d %g\n", i, t, mgp.i, mgu.i, du);
	int n_cells = 0;
	foreach(reduction(+:n_cells)){
		n_cells += 1;
	}
	timing s = timer_timing(perf.gt, i, perf.tnc, NULL);
	fprintf(stderr, "%06d %8.4g ;; %9.5g %9.7g ;; %5d %5d %7d\n", i, t, normf(u.x).rms, s.cpu, mgu.i, mgp.i, n_cells);
}

event images(t += 0.01, t <= t_end) {
	scalar u_r[], l[];
	foreach() {
		u_r[] = sqrt(sq(u.x[]) + sq(u.y[]));
		l[] = level;
	}

	double depth = length * epsilon;

	// output_ppm(u_r, file = "movie-u.mp4",   box={{0.,0.}, {length,depth*2}}, linear=true);
	// output_ppm(u.x, file = "movie-u_x.mp4", box={{0.,0.}, {length,depth*2}}, linear=true);
	// output_ppm(u.y, file = "movie-u_y.mp4", box={{0.,0.}, {length,depth*2}}, linear=true);
	// output_ppm(f,   file = "movie-f.mp4",   box={{0.,0.}, {length,depth*2}}, linear=true);
	// output_ppm(p,   file = "movie-p.mp4",   box={{0.,0.}, {length,depth*2}}, linear=true);
	// output_ppm(l,   file = "movie-l.mp4",   box={{0.,0.}, {length,depth*2}}, linear=true);

	output_ppm(u_r, file = "movie-u.mp4",   box={{0.,0.}, {length,depth*2}});
	output_ppm(u.x, file = "movie-u_x.mp4", box={{0.,0.}, {length,depth*2}});
	output_ppm(u.y, file = "movie-u_y.mp4", box={{0.,0.}, {length,depth*2}});
	output_ppm(f,   file = "movie-f.mp4",   box={{0.,0.}, {length,depth*2}});
	output_ppm(p,   file = "movie-p.mp4",   box={{0.,0.}, {length,depth*2}});
	output_ppm(l,   file = "movie-l.mp4",   box={{0.,0.}, {length,depth*2}});
	
	// Nice video of the fluid fraction
	view(tx=-0.5, ty = 0.,
		 fov = 2,
		 width = 1080, height = 108);

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
  adapt_wavelet ({sf1}, (double[]){1e-5}, minlevel = max(3, level_max - 7), maxlevel = level_max);
}
#endif

event fields(t += 0.1) {

	// Name of directory for field outputs and other variables
	char* out_dir = "fields";
	char* spec_dir;
	char* file_name;
	int file_name_len;
	FILE* fp;

	// Interface positions
	spec_dir = "interface";

	file_name_len = snprintf(NULL, 0, "%s/%s/%s-%05.2f.dat", out_dir, spec_dir, spec_dir, t)+1;
	file_name = malloc(file_name_len);
	snprintf(file_name, file_name_len, "%s/%s/%s-%05.2f.dat", out_dir, spec_dir, spec_dir, t);

	fp = fopen(file_name, "w");
	output_facets(f, fp);
	fclose(fp);

	// Pressure field
	spec_dir = "pressure";

	file_name_len = snprintf(NULL, 0, "%s/%s/%s-%05.2f.dat", out_dir, spec_dir, spec_dir, t)+1;
	file_name = malloc(file_name_len);
	snprintf(file_name, file_name_len, "%s/%s/%s-%05.2f.dat", out_dir, spec_dir, spec_dir, t);

	fp = fopen(file_name, "w");
	for(double xx = 0.; xx <= length; xx += epsilon * 0.05) {
		for(double yy = 0.; yy <= epsilon; yy += epsilon * 0.05) {
			fprintf(fp, "%f %f %f\n", xx, yy, interpolate(p, xx, yy));
		}
	}
	fclose(fp);

}
