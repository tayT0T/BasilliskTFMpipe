#include "grid/multigrid3D.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "reduced.h"
#include "tension.h"
#include "view.h"
#include "lambda2.h"
#include "maxruntime.h"
#include "navier-stokes/perfs.h"

double U1s = 0.023097;          // superficial velocity of phase 1 (heavier phase)                
double U2s = 0.3217195;         // superficial velocity of phase 2 (lighter phase)
double h_L_D = 0.2;             // Liquid holdup  
double diameter = 1.0;          // pipe diameter 
double pipe_length = 3.0;       // pipe length (Shouldn't be integer - interface intersect the grid)
double Eo = 40.0;               // Eotvos number
double We = 4.0;                // Weber number 
double density_1 = 1000.0;      // Water density 
double density_2 = 1.0;         // Air density 
double dyn_visc_1 = 0.001;      // Water dynamic viscosity 
double dyn_visc_2 = 1.818e-5;   // Air dynamic viscosity 

int LEVEL = 7;                  // Level of grid refinement 
scalar f0[];                    // volume fraction initially 

double liquid_area();           

int main(){

    dimensions (nx = diameter+0.2, ny = diameter+0.2, nz = pipe_length);    // domain size 
    init_grid (32);               // grid number without refinement 

    rho1 = density_1/density_2;   // Scaled density of phase 1  
    rho2 = 1.0;                   // Scaled density of phase 2 
    mu1 = dyn_visc_1/dyn_visc_2;  // Scaled dynamic viscosity of phase 1  
    mu2 = 1.0;                    // Scaled dynamic viscosity of phase 2  
           
    origin (-(diameter+0.2)/2, -(diameter+0.2)/2,0);            // center point

    double U_m = U1s + U2s;                                          // mixture velocity 
    double FROUDE = (We/Eo)*((density_1-density_2)/density_2);       // Froude Number 
    f.sigma = 1./We;                           // Surface tension coefficient sigma 
    G.y = - sq(U_m/FROUDE)/(diameter);         // gravitational force 
    run();
}

// Flow rate condition at the inlet 
//u.n[front] = dirichlet(U1s*f[] + U2s*(1-f[]));
u.n[front] = dirichlet( (U1s*liquid_area()/pi)*f[] + (U2s*(pi-liquid_area()/pi)*(1-f[]))); // liquid-vel = superfic * holdup
p[front]   = neumann(0.);
pf[front]  = neumann(0.);

// outlet BC
u.n[back]  = neumann(0.);
p[back]    = dirichlet(0.);
pf[back]   = dirichlet(0.);

// No slip, no penetration BC
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.r[embed] = dirichlet(0.);

event init (t = 0) {
  solid(cs,fs, -sq(x) - sq(y) + pow(diameter/2,2));         // Define solid pipe geometry
  fractions_cleanup (cs, fs);

  fraction(f0, y<(h_L_D - diameter/2) ? 
          -sq(x) - sq(y) + pow(diameter/2,2) :-1);    // initialize liquid holdup

  foreach() {
    f[] = f0[];                         // Initialize volume fraction at initial time 
    u.z[] = (U1s*liquid_area()/pi)*f0[] + (U2s*(pi-liquid_area()/pi)*(1-f0[]));    // Initialize inlet BC at initial time 
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

event end(i=10){
}

double liquid_area() {
    double Area = 0.;
    foreach_boundary (front){
        Area += f[] * sq(Delta);
    }
  return Area;
}
