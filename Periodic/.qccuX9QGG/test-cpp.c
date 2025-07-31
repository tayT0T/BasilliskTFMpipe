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
#define dimension 2
#include "common.h"
#ifndef BASILISK_HEADER_0
#define BASILISK_HEADER_0
#line 1 "test.c"
#include "grid/multigrid.h"

int main() {
  init_grid(32);
  run();
}

event init (t = 0) {
  foreach()
    printf("%g %g %g\n", x, y, z);
}

#endif
