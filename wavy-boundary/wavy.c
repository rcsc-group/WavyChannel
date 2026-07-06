#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "contact-embed.h"
#include "two-phase.h"
#include "tension.h"
#include "utils.h"

// Mesh refinement level
int level = 8;

// Domain dimensions
double length = 10., depth = 1.;

// Contact angle
double theta0 = 30.;

// Boundary conditions
u.t[embed] = neumann(0.); // Even though it's set to be free slip, the u_x profile shows no slip
u.n[embed] = dirichlet(0.);

p[left] = dirichlet(0.);
pf[left] = dirichlet_face(0.); // Unsure if dirichlet_face is required here, as opposed to dirichlet

p[right] = dirichlet(0.);
pf[right] = dirichlet_face(0.);

f[embed] = neumann(0.);
f[bottom] = neumann(0.);
f[left] = dirichlet(1.);
f[right] = dirichlet(0.);

event init(t=0) {
	solid(cs, fs, 1.e-6 + 1. - y);
	fraction(f, 1. - x);
}

event logfile(i += 10, t <= 100.) {
	timing s = timer_timing(perf.gt, i, perf.tnc, NULL);
	fprintf(stderr, "%03d %8.3g ;; %9.5g %9.7g ;; %5d %5d\n", i, t, normf(u.x).rms, s.cpu, mgu.i, mgp.i);
}

// Movie frame is a bit bigger than the fluid domain
event images(t += 0.2) {
	output_ppm(u.x, file = "u_x.mp4", linear=true, box={{0.,0.}, {length,depth+0.1}});
	output_ppm(u.y, file = "u_y.mp4", linear=true, box={{0.,0.}, {length,depth+0.1}});
	output_ppm(f,   file = "f.mp4",   linear=true, box={{0.,0.}, {length,depth+0.1}});
	output_ppm(p,   file = "p.mp4",   linear=true, box={{0.,0.}, {length,depth+0.1}});
}

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

	// Surface tension and contact angle
	f.sigma = 0.1;
	const scalar c[] = theta0 * pi / 180.;
	contact_angle = c;

	// Runtime statistics
	fprintf(stderr, "\n\nLevel: %d\n%3s %8s ;; %9s %9s ;; %5s %5s\n", level, "i", "t", "L^2 norm", "CPU time", "mgu.i", "mgp.i");

	init_grid(1 << level);
	run();
	
}
