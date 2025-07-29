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

// Input Parameter
double U1s = 0.023097;          // superficial velocity of phase 1 (heavier phase)                
double U2s = 0.3217195;         // superficial velocity of phase 2 (lighter phase)
double h_L_D = 0.2;             // Liquid height to diameter ratio
double pressure_drop = -10.0;    // Pressure drop correponding to the liquid holdup
double diameter = 1.0;          // pipe diameter 
double pipe_length = 3.0;       // pipe length (Shouldn't be integer - interface intersect the grid)
double Eo = 40.0;               // Eotvos number
double We = 4.0;                // Weber number 
double density_1 = 1000.0;      // Water density 
double density_2 = 1.0;         // Air density 
double dyn_visc_1 = 0.001;      // Water dynamic viscosity 
double dyn_visc_2 = 1.818e-5;   // Air dynamic viscosity 
double U_m
double FROUDE

// Simulation Parameter
int LEVEL = 7;                  // Level of grid refinement 
scalar f0[];                    // volume fraction initially 
face vector av[];               // forcing term via pressure drop 

// Function 
double liquid_area();  

// Boundary Condition 
// No slip, no penetration BC on the embedded
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.r[embed] = dirichlet(0.);

int main(){

    dimensions (nx = diameter*2.0, ny = diameter*2.0, nz = pipe_length);    // domain size 
    init_grid (32);               // grid number without refinement 

    rho1 = density_1/density_2;   // Scaled density of phase 1  
    rho2 = 1.0;                   // Scaled density of phase 2 
    mu1 = dyn_visc_1/dyn_visc_2;  // Scaled dynamic viscosity of phase 1  
    mu2 = 1.0;                    // Scaled dynamic viscosity of phase 2  
           
    origin (-(diameter*2.0)/2, -(diameter*2.0)/2,0);            // center point

    U_m = U1s + U2s;                                          // mixture velocity 
    FROUDE = (We/Eo)*((density_1-density_2)/density_2);       // Froude Number 
    f.sigma = 1./We;                           // Surface tension coefficient sigma 
//    G.y = - sq(U_m/FROUDE)/(diameter);         // gravitational acceleration
    a = av;                        // acceleration term in NS eq

    periodic(front);               // Periodic BC at the axial direction
    run();
}

// Initial condition
event init (t = 0) {
  solid(cs,fs, -sq(x) - sq(y) + pow(diameter/2,2));         // Define solid pipe geometry
  fractions_cleanup (cs, fs);

  fraction(f0, y<(h_L_D - diameter/2) ? 
          -sq(x) - sq(y) + pow(diameter/2,2) :-1);    // initialize liquid holdup

  foreach() {
    f[] = f0[];                         // Initialize volume fraction at initial time 
    u.z[] = U1s*f0[] + U2s*(1-f0[]);    // Initialize inlet BC at initial time 
  }

  view(camera="front",fov=0,tx=0,ty=0);
  clear();
  draw_vof("cs",filled=-1,fc={1,1,1});
  draw_vof("f0",filled=-1,fc={1,0,1});
  cells();
  save("grid_t0.jpg"); 

  view(fov=0,tx=0,ty=0, theta = 3.14/2 + 0.5,  phi = 0.4,  psi = 0.);
  clear();
  draw_vof("f",filled=-1,fc={1,0,1});
  cells();
  save("grid_1.jpg"); 
  printf("Hi!");    
}

event logfile (i++){
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
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
  save("c.mp4");

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
  save("yayayaya.mp4");
}

event end(i=30){
}

double liquid_area() {
    double Area = 0.;
    foreach_boundary (front){
        Area += f[] * sq(Delta);
    }
  return Area;
}

event acceleration (i++){ 
  foreach_face(y){
    av.y[] = - sq(U_m/FROUDE)/(diameter);         // gravitational acceleration
  }
  foreach_face(z){
    av.z[] = pressure_drop;
  }
}

#endif
