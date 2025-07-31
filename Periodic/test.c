#include "grid/multigrid.h"

int main() {
  init_grid(32);
  run();
}

event init (t = 0) {
  foreach()
    printf("%g %g %g\n", x, y, z);
}
