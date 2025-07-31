@if _XOPEN_SOURCE < 700
  @undef _XOPEN_SOURCE
  @define _XOPEN_SOURCE 700
@endif
@if _GNU_SOURCE
@include <stdint.h>
@include <string.h>
@include <fenv.h>
@endif
#define _CATCH
#define dimension 3
#define BGHOSTS 2
#include "common.h"
#ifndef BASILISK_HEADER_0
#define BASILISK_HEADER_0
#line 1 "BaseFlow.c"
#include "grid/multigrid3D.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
//#include "reduced.h"            // Re-express the gravity as the interfacial force in 2-phases
#include "tension.h"
#include "view.h"
#include "lambda2.h"
#include "maxruntime.h"
#include "navier-stokes/perfs.h"

// Constant Parameter
#define density_ratio 1000.0  // Density ratio 
#define viscosity_ratio 54.0  // Viscosity ratio
#define eotvos 79.0           // Eotvos number
#define diameter 1.0          // pipe diameter 
#define pipe_length 2.0      // pipe length (Shouldn't be integer - interface intersect the grid)

// Input Parameter
double froude_liquid = 0.02;     // Froude number on liquid phase
double froude_gas = 4.12;        // Froude number on gas phase
double reynold_liquid = 240.0;   // Reynold number on liquid phase
double reynold_gas = 3284.0;     // Reynold number on gas phase

// Guessing value 
double h_L_D_init = 0.2;             // Liquid height to diameter ratio
double forcing = 0.1;                   // forcing term
double U1s_init = 0.01;         // Initialize liquid superficial velocity
double U2s_init = 2.0;          // Inititalize gas superficial velocity 

// Simulation Parameter 
scalar f0[];                    // volume fraction initially
double U_LS;
double U_GS;
double liquid_area;
double froude_liquid_instant;
double FROUDE;

// Function declaration
double calculate_froude_liquid(double U_LS_local);
void calculate_superficial_vel(double *U_LS, double *U_GS, double *liquid_area);

// Boundary Condition 
// No slip, no penetration BC on the embedded
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.r[embed] = dirichlet(0.);

int main(){
  size(pipe_length);                 // setting the physical size if the x-dir (axial dir)
    dimensions (nx = pipe_length, ny = diameter * 2, nz = diameter * 2);    // domain size 
    init_grid (512);               // grid number without refinement
    
    rho1 = density_ratio;         // Scaled density of phase 1 (Liquid)  
    rho2 = 1.0;                   // Scaled density of phase 2 (Gas)
    mu1 = viscosity_ratio;        // Scaled dynamyic viscosity of phase 1 (Liquid)  
    mu2 = 1.0;                    // Scaled dynamic viscosity of phase 2 (Gas)
           
    origin (0, -(diameter)/2.0, -(diameter)/2.0);      // center point
    FROUDE = froude_liquid;                        // Initialize Froude Number 
    f.sigma = (1./eotvos) * (density_ratio - 1.0);                           // Surface tension coefficient sigma 

    periodic(right);               // Periodic BC at the axial direction
    run();
}

// Initial condition
event init (t = 0) {
  solid(cs,fs, -sq(y) - sq(z) + pow(diameter/2,2));         // Define solid pipe geometry
  fractions_cleanup (cs, fs);

  fraction(f0, y < (h_L_D_init - diameter/2) ? 
          -sq(y) - sq(z) + pow(diameter/2,2) :-1);    // initialize liquid holdup

  foreach() {
    f[] = f0[];                         // Initialize volume fraction at initial time 
    u.x[] = U1s_init*f[] + U2s_init*(1-f[]);    // Initialize inlet BC at initial time 
  }

  view(camera="left",fov=0,tx=0,ty=0);
  clear();
  draw_vof("cs",filled=-1,fc={1,1,1});
  draw_vof("f0",filled=-1,fc={1,0,1});
  cells();
  save("grid_t0_left.jpg"); 

  view(camera="top",fov=0,tx=0,ty=0);
  clear();
  draw_vof("cs",filled=-1,fc={1,1,1});
  draw_vof("f0",filled=-1,fc={1,0,1});
  cells();
  save("grid_t0_top.jpg"); 
  
  view(fov=0,tx=0,ty=0, theta = 3.14/2 + 0.5,  phi = 0.4,  psi = 0.);
  clear();
  draw_vof("f",filled=-1,fc={1,0,1});
  cells();
  save("grid_t0_side.jpg"); 
  printf("Hi!");    
}

event solute_movie (i += 2) {
  view(camera="front",fov=0,tx=0,ty=0);
  clear();
  draw_vof("cs",filled=1,fc={1,1,1});
  draw_vof("f", lw=3,lc={1,1,0}, min = -0.1, max = 0.1);
  scalar vel[];
  foreach()
    vel[] = sqrt(sq(u.x[])+sq(u.y[]));
  squares("vel",min=0.,max=1.,map=jet);
  cells();
  save("velocity_xy.mp4");

  view(camera="front",fov=0,tx=0,ty=0);
  clear();
  draw_vof("cs",filled=1,fc={1,1,1});
  draw_vof("f", lw=3,lc={1,1,0}, min = -0.1, max = 0.1);
  foreach()
    vel[] = u.z[];
  squares("vel",min=0.,max=1.,map=jet);
  cells();
  save("uz.mp4");


  view(fov=0,tx=0,ty=0, theta =-(3.14/2 + 0.2) ,  phi = 0.4,  psi = 0.);
  clear();
  draw_vof("f", lw=3,lc={1,1,0}, min = -0.1, max = 0.1);
  save("front_vof_vid.mp4");
}

event end(i=10){
}

event acceleration (i++){
  face vector av = a;                        // acceleration term in vertical direction
  foreach_face(y){
    av.z[] -= 1.0;         // gravitational acceleration
  }
  foreach_face(x){
    av.x[] += forcing;                        // forcing term 
  }
}

event logfile (i += 2) {
  calculate_superficial_vel(&U_LS, &U_GS, &liquid_area);
  fprintf(stderr, "t = %g, U_LS = %g, U_GS = %g, liquid area = %g\n",
          t, U_LS, U_GS, liquid_area);
}
  
void calculate_superficial_vel(double *U_LS, double *U_GS, double *liquid_area) {
  *U_LS = 0.0;
  *U_GS = 0.0;
  *liquid_area = 0.0;

  foreach_boundary(right) {
    *U_LS += f[] * u.x[] * sq(Delta) / (M_PI * sq(0.5 * diameter));
    *U_GS += (1. - f[]) * u.x[] * sq(Delta) / (M_PI * sq(0.5 * diameter));
    *liquid_area += f[] * sq(Delta);
  }
}

#endif
