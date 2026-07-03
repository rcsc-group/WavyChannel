#include "embed.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "two-phase.h"
#include "tension.h"
#include "utils.h"

int level = 8;

double length = 10., depth = 1.;
double finalTime = 5.;

// Height field for contact angle
vector h[];

////// Boundary conditions
// Free slip, no flux
u.t[embed] = neumann(0.);
u.n[embed] = dirichlet(0.);

// Inlet and outlet pressures
p[left] = dirichlet(10.);
pf[left] = dirichlet(10.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

// Explicit conditions on tracer, including contact angle
h.t[embed] = contact_angle(60 * pi/180.);
f[embed] = neumann(0.);
f[bottom] = neumann(0.);
f[left] = dirichlet(1.);
f[right] = dirichlet(0.);

event init(t=0) {
	solid(cs, fs, depth - y + 1.e-6);
	fraction(f, 1. - x);
}

// Runtime statistics
event logfile(i += 10, t <= finalTime) {
	timing s = timer_timing(perf.gt, i, perf.tnc, NULL);
	fprintf(stderr, "%03d %8.3g ;; %9.5g %9.7g ;; %5d %5d\n", i, t, normf(u.x).rms, s.cpu, mgu.i, mgp.i);
}

// Movie output
event images(t += 0.1) {
	output_ppm(u.x, file = "u_x.mp4", linear=true, box={{0.,0.}, {length,depth}});
	output_ppm(u.y, file = "u_y.mp4", linear=true, box={{0.,0.}, {length,depth}});
	output_ppm(f,   file = "f.mp4",   linear=true, box={{0.,0.}, {length,depth}});
	output_ppm(p,   file = "p.mp4",   linear=true, box={{0.,0.}, {length,depth}});
}

// Storing a few profiles of quantities of interest from different sections of the domain
event profiles(t = end) {

	FILE* fp = fopen("u_x_prof.dat", "w");
	for(double yy = 0.; yy <= depth; yy += depth * 0.01) {
		fprintf(fp, "%g %g %g %g\n", yy, interpolate(u.x, 0.25 * length, yy), interpolate(u.x, 0.5 * length, yy), interpolate(u.x, 0.75 * length, yy));
	}
	fclose(fp);

	fp = fopen("p_prof.dat", "w");
	for(double xx = 0.; xx <= length; xx += length * 0.01) {
		fprintf(fp, "%g %g %g %g\n", xx, interpolate(p, xx, 0.25 * depth), interpolate(p, xx, 0.5 * depth), interpolate(p, xx, 0.75 * depth));
	}
	fclose(fp);

	fp = fopen("f_prof.dat", "w");
	for(double xx = 0.; xx <= length; xx += length * 0.01) {
		fprintf(fp, "%g %g %g %g\n", xx, interpolate(f, xx, 0.25 * depth), interpolate(f, xx, 0.5 * depth), interpolate(f, xx, 0.75 * depth));
	}
	fclose(fp);

}

int main() {

	size(length);
	origin(0., 0.);
	
	// 1: left
	// 2: right
	rho1 = 1., rho2 = 1.;
	mu1 = 1., mu2 = 1.;

	f.height = h;
	f.sigma = 0.1;

	// Runtime statistics
	fprintf(stderr, "\n\nLevel: %d\n%3s %8s ;; %9s %9s ;; %5s %5s\n", level, "i", "t", "L^2 norm", "CPU time", "mgu.i", "mgp.i");

	init_grid(1 << level);
	run();
}
