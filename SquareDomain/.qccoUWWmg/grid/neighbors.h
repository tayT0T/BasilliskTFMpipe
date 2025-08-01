#ifndef BASILISK_HEADER_56
#define BASILISK_HEADER_56
#line 1 "/home/pwachara/basilisk/src/grid/neighbors.h"
#if dimension == 1
macro1 foreach_neighbor (int _s = GHOSTS,
			 Point point = point, break = (_k = _nn + 1)) {
  {
    const int _nn = _s;
    const int _i = point.i;
    for (int _k = - _nn; _k <= _nn; _k++) {
      point.i = _i + _k;
      POINT_VARIABLES();
      {...}
    }
    point.i = _i;
  }
}
#elif dimension == 2
macro1 foreach_neighbor (int _s = GHOSTS,
			 Point point = point, break = (_k = _l = _nn + 1)) {
  {
    const int _nn = _s;
    const int _i = point.i, _j = point.j;
    for (int _k = - _nn; _k <= _nn; _k++) {
      point.i = _i + _k;
      for (int _l = - _nn; _l <= _nn; _l++) {
	point.j = _j + _l;
	POINT_VARIABLES();
	{...}
      }
    }
    point.i = _i; point.j = _j;
  }
}
#else // dimension == 3
macro1 foreach_neighbor (int _s = GHOSTS,
			 Point point = point, break = (_l = _m = _n = _nn + 1)) {
  {
    const int _nn = _s;
    const int _i = point.i, _j = point.j, _k = point.k;
    for (int _l = - _nn; _l <= _nn; _l++) {
      point.i = _i + _l;
      for (int _m = - _nn; _m <= _nn; _m++) {
	point.j = _j + _m;
	for (int _n = - _nn; _n <= _nn; _n++) {
	  point.k = _k + _n;
	  POINT_VARIABLES();
	  {...}
	}
      }
    }
    point.i = _i; point.j = _j; point.k = _k;
  }
}
#endif // dimension == 3

#endif
