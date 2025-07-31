#include "grid/multigrid3D.h"
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
#define width 1.0          // square channel width  
#define pipe_length 3.0      // pipe length 

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
void calculate_superficial_vel(double *U_LS, double *U_GS, double *liquid_area);

// Boundary Condition 
// No slip, no penetration BC on the square pipe
u.n[front] = dirichlet(0.);
u.t[front] = dirichlet(0.);
u.r[front] = dirichlet(0.);

u.n[back] = dirichlet(0.);
u.t[back] = dirichlet(0.);
u.r[back] = dirichlet(0.);

u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);
u.r[top] = dirichlet(0.);

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.r[bottom] = dirichlet(0.);

int main(){
  size(pipe_length);                 // setting the physical size if the x-dir (axial dir)
    dimensions (nx = pipe_length, ny = width, nz = width);    // domain size 
    init_grid (512);               // grid number without refinement
    
    rho1 = density_ratio;         // Scaled density of phase 1 (Liquid)  
    rho2 = 1.0;                   // Scaled density of phase 2 (Gas)
    mu1 = viscosity_ratio;        // Scaled dynamyic viscosity of phase 1 (Liquid)  
    mu2 = 1.0;                    // Scaled dynamic viscosity of phase 2 (Gas)
           
    origin (0, -width/2.0, -width/2.0);      // center point
    FROUDE = froude_liquid;                        // Initialize Froude Number 
    f.sigma = (1./eotvos) * (density_ratio - 1.0);                           // Surface tension coefficient sigma 

    periodic(right);               // Periodic BC at the axial direction
    run();
}

// Initial condition
event init (t = 0) {

  fraction(f0, -y + (h_L_D_init - width/2));    // initialize liquid holdup

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
    *U_LS += f[] * u.x[] * sq(Delta) / (M_PI * sq(width));
    *U_GS += (1. - f[]) * u.x[] * sq(Delta) / (M_PI * sq(width));
    *liquid_area += f[] * sq(Delta);
  }
}
