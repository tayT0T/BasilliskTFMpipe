#line 0 "BaseFlow-cpp.c"
#line 0 "<built-in>"
#line 0 "<command-line>"
#line 1 "/usr/include/stdc-predef.h"
#line 0 "<command-line>"
#line 1 "BaseFlow-cpp.c"
#if _XOPEN_SOURCE < 700
  #undef _XOPEN_SOURCE
  #define _XOPEN_SOURCE 700
#endif
#if _GNU_SOURCE
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#endif



#line 1 "/home/pwachara/basilisk/src/common.h"
typedef double number; 






number macro_min (number a, number b) { return a < b ? a : b; } 
number macro_sq (number x) { return x*x; } 
number macro_cube (number x) { return x*x*x; } 
int macro_sign (number x) { return (int)(x > 0 ? 1 : -1); } 
int macro_sign2 (number x) { return (int)(x > 0 ? 1 : x < 0 ? -1 : 0); } 
number macro_clamp (number x, number a, number b) {
  return x < a ? a : x > b ? b : x;
} 



#line 1 "/home/pwachara/basilisk/src/grid/config.h"
#line 1 "/home/pwachara/basilisk/src/ast/std/stdlib.h"
#include <stdlib.h>
#line 2 "/home/pwachara/basilisk/src/grid/config.h"
#line 1 "/home/pwachara/basilisk/src/ast/std/stdio.h"
#include <stdio.h>
#line 3 "/home/pwachara/basilisk/src/grid/config.h"
#line 1 "/home/pwachara/basilisk/src/ast/std/stddef.h"
#include <stddef.h>
#line 4 "/home/pwachara/basilisk/src/grid/config.h"
#line 1 "/home/pwachara/basilisk/src/ast/std/stdint.h"
#include <stdint.h>
#line 5 "/home/pwachara/basilisk/src/grid/config.h"
#line 1 "/home/pwachara/basilisk/src/ast/std/stdbool.h"
#include <stdbool.h>
#line 6 "/home/pwachara/basilisk/src/grid/config.h"
#line 1 "/home/pwachara/basilisk/src/ast/std/stdarg.h"
#include <stdarg.h>
#line 7 "/home/pwachara/basilisk/src/grid/config.h"
#line 1 "/home/pwachara/basilisk/src/ast/std/string.h"
#include <string.h>
#line 8 "/home/pwachara/basilisk/src/grid/config.h"
#line 1 "/home/pwachara/basilisk/src/ast/std/float.h"
#include <float.h>
#line 9 "/home/pwachara/basilisk/src/grid/config.h"
#line 1 "/home/pwachara/basilisk/src/ast/std/limits.h"
#include <limits.h>
#line 10 "/home/pwachara/basilisk/src/grid/config.h"
#line 1 "/home/pwachara/basilisk/src/ast/std/math.h"
#include <math.h>
#line 11 "/home/pwachara/basilisk/src/grid/config.h"
#line 1 "/home/pwachara/basilisk/src/ast/std/time.h"
#include <time.h>
#line 12 "/home/pwachara/basilisk/src/grid/config.h"
#line 1 "/home/pwachara/basilisk/src/ast/std/sys/time.h"
#include <sys/time.h>
#line 13 "/home/pwachara/basilisk/src/grid/config.h"
#line 1 "/home/pwachara/basilisk/src/ast/std/sys/resource.h"
#include <sys/resource.h>
#line 14 "/home/pwachara/basilisk/src/grid/config.h"




#define unmap(x,y)
#define trash(x)


void macro_BEGIN_FOREACH() {;}
#line 41 "/home/pwachara/basilisk/src/grid/config.h"
# define OMP(x)
void macro_OMP_SERIAL() {;}

# include <mpi.h>
static int mpi_rank, mpi_npe;
# define tid() mpi_rank
# define pid() mpi_rank
# define npe() mpi_npe
#line 69 "/home/pwachara/basilisk/src/grid/config.h"
#define _NVARMAX 65536
#define is_constant(v) ((v).i >= _NVARMAX)
#define constant(v) (is_constant(v) ? _constant[(v).i - _NVARMAX] : 1e30f)

#define systderr stderr
#define systdout stdout

FILE * qstderr (void);
FILE * qstdout (void);
FILE * ferr = NULL, * fout = NULL;
#define not_mpi_compatible()\
do {\
  if (npe() > 1) {\
    fprintf (ferr, "%s() is not compatible with MPI (yet)\n", __func__);\
    exit (1);\
  }\
} while(0)\

#line 86

# define system(command) (pid() == 0 ? system(command) : 0)
#line 97 "/home/pwachara/basilisk/src/grid/config.h"
static inline void qassert (const char * file, int line, const char * cond) {
  fprintf (ferr, "%s:%d: Assertion `%s' failed.\n", file, line, cond);
  abort();
}







#define sysmalloc malloc
#define syscalloc calloc
#define sysrealloc realloc
#define sysfree free
#define systrdup strdup




# define pmalloc(s,func,file,line) malloc(s)
# define pcalloc(n,s,func,file,line) calloc(n,s)
# define prealloc(p,s,func,file,line) realloc(p,s)
# define pfree(p,func,file,line) free(p)
# define pstrdup(s,func,file,line) strdup(s)






#line 1 "/home/pwachara/basilisk/src/grid/array.h"


typedef struct {
  void * p;
  long max, len;
} Array;

Array * array_new()
{
  Array * a = ((Array *) pmalloc ((1)*sizeof(Array),__func__,__FILE__,__LINE__));
  a->p = NULL;
  a->max = a->len = 0;
  return a;
}

void array_free (Array * a)
{
  pfree (a->p,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
}

void * array_append (Array * a, void * elem, size_t size)
{
  if (a->len + size >= a->max) {
    a->max += ( size > 4096 ? size : 4096);
    a->p = prealloc (a->p, a->max,__func__,__FILE__,__LINE__);
  }
  memcpy (((char *)a->p) + a->len, elem, size);
  a->len += size;
  return (void *)(((char *)a->p) + a->len - size);
}

void * array_shrink (Array * a)
{
  void * p = prealloc (a->p, a->len,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
  return p;
}
#line 129 "/home/pwachara/basilisk/src/grid/config.h"
#line 372 "/home/pwachara/basilisk/src/grid/config.h"
# define tracing(...)
# define end_tracing(...)
#line 388 "/home/pwachara/basilisk/src/grid/config.h"
static bool in_prof = false;
static double prof_start, _prof;
#define prof_start(name)\
  if (!(!in_prof)) qassert ("/home/pwachara/basilisk/src/grid/config.h", 391, "!in_prof"); in_prof = true;\
  prof_start = MPI_Wtime();\

#line 393

#define prof_stop()\
  if (!(in_prof)) qassert ("/home/pwachara/basilisk/src/grid/config.h", 395, "in_prof"); in_prof = false;\
  _prof = MPI_Wtime();\
  mpi_time += _prof - prof_start;\

#line 398






     
int mpi_all_reduce0 (void *sendbuf, void *recvbuf, int count,
       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{tracing("mpi_all_reduce0","/home/pwachara/basilisk/src/grid/config.h",405);
  { int _ret= MPI_Allreduce (sendbuf, recvbuf, count, datatype, op, comm);end_tracing("mpi_all_reduce0","/home/pwachara/basilisk/src/grid/config.h",408);return _ret;}
end_tracing("mpi_all_reduce0","/home/pwachara/basilisk/src/grid/config.h",409);}

#define mpi_all_reduce(v,type,op) {\
  prof_start ("mpi_all_reduce");\
  union { int a; float b; double c;} global;\
  mpi_all_reduce0 (&(v), &global, 1, type, op, MPI_COMM_WORLD);\
  memcpy (&(v), &global, sizeof (v));\
  prof_stop();\
}\

#line 418


     
void mpi_all_reduce_array (void * v, MPI_Datatype datatype, MPI_Op op, int elem)
{tracing("mpi_all_reduce_array","/home/pwachara/basilisk/src/grid/config.h",421);
  size_t size;
  if (datatype == MPI_DOUBLE) size = sizeof (double);
  else if (datatype == MPI_INT) size = sizeof (int);
  else if (datatype == MPI_LONG) size = sizeof (long);
  else if (datatype == MPI_C_BOOL) size = sizeof (bool);
  else if (datatype == MPI_UNSIGNED_CHAR) size = sizeof (unsigned char);
  else {
    fprintf (ferr, "unknown reduction type\n");
    fflush (ferr);
    abort();
  }
  void * global = pmalloc (elem*size,__func__,__FILE__,__LINE__), * tmp = pmalloc (elem*size,__func__,__FILE__,__LINE__);
  memcpy (tmp, v, elem*size);
  mpi_all_reduce0 (tmp, global, elem, datatype, op, MPI_COMM_WORLD);
  memcpy (v, global, elem*size);
  pfree (global,__func__,__FILE__,__LINE__), pfree (tmp,__func__,__FILE__,__LINE__);
end_tracing("mpi_all_reduce_array","/home/pwachara/basilisk/src/grid/config.h",439);}


#define QFILE FILE

FILE * qstderr (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "log-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systderr;
  }
  return fp;
}

FILE * qstdout (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "out-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systdout;
  }
  return fp;
}

static void finalize (void)
{
  MPI_Finalize();
}

void mpi_init()
{
  int initialized;
  MPI_Initialized (&initialized);
  if (!initialized) {
    MPI_Init (NULL, NULL);
    MPI_Comm_set_errhandler (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
    atexit (finalize);
  }
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_npe);
  srand (mpi_rank + 1);
  if (ferr == NULL) {
    if (mpi_rank > 0) {
      ferr = fopen ("/dev/null", "w");
      fout = fopen ("/dev/null", "w");
    }
    else {
      ferr = systderr;
      fout = systdout;
    }
    char * etrace = getenv ("MALLOC_TRACE"), name[80];
    if (etrace && mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      setenv ("MALLOC_TRACE", name, 1);
    }
#line 519 "/home/pwachara/basilisk/src/grid/config.h"
  }
}
#line 532 "/home/pwachara/basilisk/src/grid/config.h"
void macro_OMP_PARALLEL() {;}
#define OMP_PARALLEL(...) OMP(omp parallel __VA_ARGS__)

#define NOT_UNUSED(x) (void)(x)

void macro1_VARIABLES() { ; }
#define _index(a,m) (a.i)
#define val(a,k,l,m) data(k,l,m)[_index(a,m)]

double _val_higher_dimension = 0.;
#line 550 "/home/pwachara/basilisk/src/grid/config.h"
#if (_GNU_SOURCE || __APPLE__) && !_OPENMP && !_CADNA
double undefined;
# if __APPLE__
# include <stdint.h>
# include "fp_osx.h"
# endif
#if _GPU
# define enable_fpe(flags)
#else
# define enable_fpe(flags) feenableexcept (flags)
#endif
# define disable_fpe(flags) fedisableexcept (flags)
static void set_fpe (void) {
  int64_t lnan = 0x7ff0000000000001;
  if (!(sizeof (int64_t) == sizeof (double))) qassert ("/home/pwachara/basilisk/src/grid/config.h", 564, "sizeof (int64_t) == sizeof (double)");
  memcpy (&undefined, &lnan, sizeof (double));
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
#else
# define undefined ((double) DBL_MAX)
# define enable_fpe(flags)
# define disable_fpe(flags)
static void set_fpe (void) {}
#endif



static FILE ** qpopen_pipes = NULL;

FILE * qpopen (const char * command, const char * type)
{
  if (pid() > 0)
    return fopen ("/dev/null", type);
  FILE * fp = popen (command, type);
  if (fp) {
    FILE ** i = qpopen_pipes;
    int n = 0;
    while (i && *i) { n++; i++; }
    qpopen_pipes = (FILE * *) prealloc (qpopen_pipes, (n + 2)*sizeof(FILE *),__func__,__FILE__,__LINE__);
    qpopen_pipes[n] = fp;
    qpopen_pipes[n+1] = NULL;
  }
  return fp;
}

int qpclose (FILE * fp)
{
  if (pid() > 0)
    return fclose (fp);
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i == fp)
      *i = (FILE *) 1;
    i++;
  }
  return pclose (fp);
}

static void qpclose_all()
{
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i != (FILE *) 1)
      pclose (*i);
    i++;
  }
  pfree (qpopen_pipes,__func__,__FILE__,__LINE__);
  qpopen_pipes = NULL;
}






FILE * lfopen (const char * name, const char * mode)
{
  char fname[80];
  sprintf (fname, "%s-%d", name, pid());
  return fopen (fname, mode);
}

#line 1 "/home/pwachara/basilisk/src/grid/../ast/symbols.h"

const char * symbol_name (int sym);
int token_symbol (int token);
enum yysymbol_kind_t
{
  sym_YYEMPTY = -2,
  sym_YYEOF = 0,
  sym_YYerror = 1,
  sym_YYUNDEF = 2,
  sym_IDENTIFIER = 3,
  sym_I_CONSTANT = 4,
  sym_F_CONSTANT = 5,
  sym_STRING_LITERAL = 6,
  sym_FUNC_NAME = 7,
  sym_SIZEOF = 8,
  sym_PTR_OP = 9,
  sym_INC_OP = 10,
  sym_DEC_OP = 11,
  sym_LEFT_OP = 12,
  sym_RIGHT_OP = 13,
  sym_LE_OP = 14,
  sym_GE_OP = 15,
  sym_EQ_OP = 16,
  sym_NE_OP = 17,
  sym_AND_OP = 18,
  sym_OR_OP = 19,
  sym_MUL_ASSIGN = 20,
  sym_DIV_ASSIGN = 21,
  sym_MOD_ASSIGN = 22,
  sym_ADD_ASSIGN = 23,
  sym_SUB_ASSIGN = 24,
  sym_LEFT_ASSIGN = 25,
  sym_RIGHT_ASSIGN = 26,
  sym_AND_ASSIGN = 27,
  sym_XOR_ASSIGN = 28,
  sym_OR_ASSIGN = 29,
  sym_TYPEDEF_NAME = 30,
  sym_ENUMERATION_CONSTANT = 31,
  sym_TYPEDEF = 32,
  sym_EXTERN = 33,
  sym_STATIC = 34,
  sym_AUTO = 35,
  sym_REGISTER = 36,
  sym_INLINE = 37,
  sym_CONST = 38,
  sym_RESTRICT = 39,
  sym_VOLATILE = 40,
  sym_BOOL = 41,
  sym_CHAR = 42,
  sym_SHORT = 43,
  sym_INT = 44,
  sym_LONG = 45,
  sym_SIGNED = 46,
  sym_UNSIGNED = 47,
  sym_FLOAT = 48,
  sym_DOUBLE = 49,
  sym_VOID = 50,
  sym_COMPLEX = 51,
  sym_IMAGINARY = 52,
  sym_STRUCT = 53,
  sym_UNION = 54,
  sym_ENUM = 55,
  sym_ELLIPSIS = 56,
  sym_CASE = 57,
  sym_DEFAULT = 58,
  sym_IF = 59,
  sym_ELSE = 60,
  sym_SWITCH = 61,
  sym_WHILE = 62,
  sym_DO = 63,
  sym_FOR = 64,
  sym_GOTO = 65,
  sym_CONTINUE = 66,
  sym_BREAK = 67,
  sym_RETURN = 68,
  sym_ALIGNAS = 69,
  sym_ALIGNOF = 70,
  sym_ATOMIC = 71,
  sym_GENERIC = 72,
  sym_NORETURN = 73,
  sym_STATIC_ASSERT = 74,
  sym_THREAD_LOCAL = 75,
  sym_MAYBECONST = 76,
  sym_NEW_FIELD = 77,
  sym_TRACE = 78,
  sym_FOREACH_DIMENSION = 79,
  sym_REDUCTION = 80,
  sym_MACRO = 81,
  sym_ELLIPSIS_MACRO = 82,
  sym_MACRODEF = 83,
  sym_foreach_statement = 84,
  sym_85_ = 85,
  sym_86_ = 86,
  sym_87_ = 87,
  sym_88_ = 88,
  sym_89_ = 89,
  sym_90_ = 90,
  sym_91_ = 91,
  sym_92_ = 92,
  sym_93_ = 93,
  sym_94_ = 94,
  sym_95_ = 95,
  sym_96_ = 96,
  sym_97_ = 97,
  sym_98_ = 98,
  sym_99_ = 99,
  sym_100_ = 100,
  sym_101_ = 101,
  sym_102_ = 102,
  sym_103_ = 103,
  sym_104_ = 104,
  sym_105_ = 105,
  sym_106_ = 106,
  sym_107_ = 107,
  sym_108_ = 108,
  sym_YYACCEPT = 109,
  sym_translation_unit = 110,
  sym_primary_expression = 111,
  sym_expression_error = 112,
  sym_constant = 113,
  sym_enumeration_constant = 114,
  sym_string = 115,
  sym_generic_selection = 116,
  sym_generic_assoc_list = 117,
  sym_generic_association = 118,
  sym_postfix_expression = 119,
  sym_postfix_initializer = 120,
  sym_array_access = 121,
  sym_function_call = 122,
  sym_member_identifier = 123,
  sym_argument_expression_list = 124,
  sym_argument_expression_list_item = 125,
  sym_unary_expression = 126,
  sym_unary_operator = 127,
  sym_cast_expression = 128,
  sym_multiplicative_expression = 129,
  sym_additive_expression = 130,
  sym_shift_expression = 131,
  sym_relational_expression = 132,
  sym_equality_expression = 133,
  sym_and_expression = 134,
  sym_exclusive_or_expression = 135,
  sym_inclusive_or_expression = 136,
  sym_logical_and_expression = 137,
  sym_logical_or_expression = 138,
  sym_conditional_expression = 139,
  sym_assignment_expression = 140,
  sym_assignment_operator = 141,
  sym_expression = 142,
  sym_constant_expression = 143,
  sym_declaration = 144,
  sym_declaration_specifiers = 145,
  sym_init_declarator_list = 146,
  sym_init_declarator = 147,
  sym_storage_class_specifier = 148,
  sym_type_specifier = 149,
  sym_types = 150,
  sym_struct_or_union_specifier = 151,
  sym_struct_or_union = 152,
  sym_struct_declaration_list = 153,
  sym_struct_declaration = 154,
  sym_specifier_qualifier_list = 155,
  sym_struct_declarator_list = 156,
  sym_struct_declarator = 157,
  sym_enum_specifier = 158,
  sym_enumerator_list = 159,
  sym_enumerator = 160,
  sym_atomic_type_specifier = 161,
  sym_type_qualifier = 162,
  sym_function_specifier = 163,
  sym_alignment_specifier = 164,
  sym_declarator = 165,
  sym_direct_declarator = 166,
  sym_generic_identifier = 167,
  sym_pointer = 168,
  sym_type_qualifier_list = 169,
  sym_parameter_type_list = 170,
  sym_parameter_list = 171,
  sym_parameter_declaration = 172,
  sym_identifier_list = 173,
  sym_type_name = 174,
  sym_abstract_declarator = 175,
  sym_direct_abstract_declarator = 176,
  sym_type_not_specified = 177,
  sym_initializer = 178,
  sym_initializer_list = 179,
  sym_designation = 180,
  sym_designator_list = 181,
  sym_designator = 182,
  sym_static_assert_declaration = 183,
  sym_statement = 184,
  sym_labeled_statement = 185,
  sym_compound_statement = 186,
  sym_187_1 = 187,
  sym_block_item_list = 188,
  sym_block_item = 189,
  sym_expression_statement = 190,
  sym_selection_statement = 191,
  sym_for_scope = 192,
  sym_iteration_statement = 193,
  sym_for_declaration_statement = 194,
  sym_jump_statement = 195,
  sym_external_declaration = 196,
  sym_function_declaration = 197,
  sym_function_definition = 198,
  sym_declaration_list = 199,
  sym_basilisk_statements = 200,
  sym_macro_statement = 201,
  sym_reduction_list = 202,
  sym_reduction = 203,
  sym_reduction_operator = 204,
  sym_reduction_array = 205,
  sym_foreach_dimension_statement = 206,
  sym_forin_declaration_statement = 207,
  sym_forin_statement = 208,
  sym_forin_arguments = 209,
  sym_event_definition = 210,
  sym_event_parameters = 211,
  sym_event_parameter = 212,
  sym_boundary_definition = 213,
  sym_external_foreach_dimension = 214,
  sym_attribute = 215,
  sym_new_field = 216,
  sym_root = 217
};
#line 633 "/home/pwachara/basilisk/src/grid/config.h"

enum typedef_kind_t {
  sym_SCALAR = sym_root + 1,
  sym_VECTOR,
  sym_TENSOR,
  sym_COORD,
  sym__COORD,
  sym_VEC4,
  sym_IVEC
};

#define attroffset(x) (offsetof(_Attributes,x))




typedef int Reduce;

void macro_foreach_face (char flags, Reduce reductions,
   const char * order)
{;}
void macro_einstein_sum() {;}
void macro_diagonalize (int a) {;}




#define dimensional(...)

#define show_dimension_internal(...)
#define display_value(...)
#define interpreter_verbosity(...)
#line 20 "/home/pwachara/basilisk/src/common.h"

static inline double noise() { return 1. - 2.*rand()/(double)RAND_MAX; }


typedef struct {
  long n;
  long tn;
  int depth;
  int maxdepth;
} Grid;
Grid * grid = NULL;

double X0 = 0., Y0 = 0., Z0 = 0.;

double L0 = 1.;




int N = 16;


typedef struct { int i; } scalar;

typedef struct {
  scalar x;

  scalar y;


  scalar z;

} vector;

typedef struct {
  scalar * x;

  scalar * y;


  scalar * z;

} vectorl;

typedef struct {
  vector x;

  vector y;


  vector z;

} tensor;

struct { int x, y, z; } Period = {false, false, false};

typedef struct {
  double x, y, z;
} coord;

OMP(omp declare reduction (+ : coord :
      omp_out.x += omp_in.x,
      omp_out.y += omp_in.y,
      omp_out.z += omp_in.z))
#line 96 "/home/pwachara/basilisk/src/common.h"
void normalize (coord * n)
{
  double norm = 0.;
  
    norm += ( (n->x)*(n->x));    norm += ( (n->y)*(n->y));    norm += ( (n->z)*(n->z));
  norm = sqrt(norm);
  
    n->x /= norm;    n->y /= norm;    n->z /= norm;
}

void origin (double x, double y, double z) {
  X0 = x; Y0 = y; Z0 = z;
}

void size (double L) {
  L0 = L;
}

double zero (double s0, double s1, double s2) { return 0.; }
#line 123 "/home/pwachara/basilisk/src/common.h"
  enum { right, left, top, bottom, front, back };

int nboundary = 2*3;



double * _constant = NULL;
size_t datasize = 0;
typedef struct _Point Point;

#line 1 "/home/pwachara/basilisk/src/grid/boundaries.h"


typedef struct _Boundary Boundary;

struct _Boundary {
  void (* destroy) (Boundary * b);
  void (* level) (const Boundary * b, scalar * list, int l);

  void (* restriction) (const Boundary * b, scalar * list, int l);
};

static Boundary ** boundaries = NULL;

void add_boundary (Boundary * b) {
  int len = 0;
  if (boundaries) {
    Boundary ** i = boundaries;
    while (*i++) len++;
  }
  boundaries = (Boundary * *) prealloc (boundaries, (len + 2)*sizeof(Boundary *),__func__,__FILE__,__LINE__);
  boundaries[len] = b;
  boundaries[len+1] = NULL;
}

void free_boundaries() {
  if (!boundaries)
    return;
  Boundary ** i = boundaries, * b;
  while ((b = *i++))
    if (b->destroy)
      b->destroy (b);
    else
      pfree (b,__func__,__FILE__,__LINE__);
  pfree (boundaries,__func__,__FILE__,__LINE__);
  boundaries = NULL;
}
#line 47 "/home/pwachara/basilisk/src/grid/boundaries.h"
typedef struct {
  Boundary parent;
  int d;
} BoxBoundary;
#line 134 "/home/pwachara/basilisk/src/common.h"



typedef struct {
  int x;

  int y;


  int z;

} ivec;
typedef double (* BoundaryFunc) (Point, Point, scalar, bool *);
typedef struct {
  BoundaryFunc * boundary;
  BoundaryFunc * boundary_homogeneous;
  double (* gradient) (double, double, double);
  void (* delete) (scalar);
  char * name;
  ivec d;
  vector v;
  int face;
  bool nodump, freed;
  int block;
  scalar * depends;  
#line 19 "/home/pwachara/basilisk/src/grid/stencils.h"
bool input, output, nowarning;
  int width;
  int dirty;  
#line 21 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
void (* prolongation) (Point, scalar);
  void (* restriction) (Point, scalar);
  
#line 178 "/home/pwachara/basilisk/src/embed.h"
bool third;  
#line 28 "/home/pwachara/basilisk/src/vof.h"
scalar * tracers, c;
  bool inverse;  
#line 21 "/home/pwachara/basilisk/src/iforce.h"
scalar phi;
  
#line 460 "/home/pwachara/basilisk/src/heights.h"
vector height;  
#line 22 "/home/pwachara/basilisk/src/tension.h"
double sigma;

#line 159 "/home/pwachara/basilisk/src/common.h"
} _Attributes;

static _Attributes * _attribute = NULL;
#line 171 "/home/pwachara/basilisk/src/common.h"
ivec Dimensions = {1,1,1};




int list_len (scalar * list)
{
  if (!list) return 0;
  int ns = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ ns++;}}
  return ns;
}

scalar * list_append (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,__LINE__);
  list[len] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_prepend (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,__LINE__);
  for (int i = len; i >= 1; i--)
    list[i] = list[i-1];
  list[0] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_add (scalar * list, scalar s)
{
  {scalar*_i=(scalar*)( list);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){
    if (t.i == s.i)
      return list;}}
  return list_append (list, s);
}

int list_lookup (scalar * l, scalar s)
{
  if (l != NULL)
    {scalar*_i=(scalar*)( l);if(_i)for(scalar s1=*_i;(&s1)->i>=0;s1=*++_i){
      if (s1.i == s.i)
 return true;}}
  return false;
}

scalar * list_copy (scalar * l)
{
  scalar * list = NULL;
  if (l != NULL)
    {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      list = list_append (list, s);}}
  return list;
}

scalar * list_concat (scalar * l1, scalar * l2)
{
  scalar * l3 = list_copy (l1);
  {scalar*_i=(scalar*)( l2);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    l3 = list_append (l3, s);}}
  return l3;
}

void list_print (scalar * l, FILE * fp)
{
  int i = 0;
  {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    fprintf (fp, "%s%s", i++ == 0 ? "{" : ",", _attribute[s.i].name);}}
  fputs (i > 0 ? "}\n" : "{}\n", fp);
}

int vectors_len (vector * list)
{
  if (!list) return 0;
  int nv = 0;
  {vector*_i=(vector*)( list);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ nv++;}}
  return nv;
}

vector * vectors_append (vector * list, vector v)
{
  int len = vectors_len (list);
  list = (vector *) prealloc (list, (len + 2)*sizeof(vector),__func__,__FILE__,__LINE__);
  list[len] = v;
  list[len + 1] = (vector){{-1}};
  return list;
}

vector * vectors_add (vector * list, vector v)
{
  {vector*_i=(vector*)( list);if(_i)for(vector w=*_i;(&w)->x.i>=0;w=*++_i){ {
    bool id = true;
    
      if (w.x.i != v.x.i)
 id = false;      
#line 268
if (w.y.i != v.y.i)
 id = false;      
#line 268
if (w.z.i != v.z.i)
 id = false;
    if (id)
      return list;
  }}}
  return vectors_append (list, v);
}

vector * vectors_copy (vector * l)
{
  vector * list = NULL;
  if (l != NULL)
    {vector*_i=(vector*)( l);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
      list = vectors_append (list, v);}}
  return list;
}

vector * vectors_from_scalars (scalar * s)
{
  vector * list = NULL;
  while (s->i >= 0) {
    vector v;
     {
      if (!(s->i >= 0)) qassert ("/home/pwachara/basilisk/src/common.h", 291, "s->i >= 0");
      v.x = *s++;
    } 
#line 290
{
      if (!(s->i >= 0)) qassert ("/home/pwachara/basilisk/src/common.h", 291, "s->i >= 0");
      v.y = *s++;
    } 
#line 290
{
      if (!(s->i >= 0)) qassert ("/home/pwachara/basilisk/src/common.h", 291, "s->i >= 0");
      v.z = *s++;
    }
    list = vectors_append (list, v);
  }
  return list;
}

int tensors_len (tensor * list)
{
  if (!list) return 0;
  int nt = 0;
  {tensor*_i=(tensor*)( list);if(_i)for(tensor t=*_i;(&t)->x.x.i>=0;t=*++_i){ nt++;}}
  return nt;
}

tensor * tensors_append (tensor * list, tensor t)
{
  int len = tensors_len (list);
  list = (tensor *) prealloc (list, (len + 2)*sizeof(tensor),__func__,__FILE__,__LINE__);
  list[len] = t;
  list[len + 1] = (tensor){{{-1}}};
  return list;
}

tensor * tensors_from_vectors (vector * v)
{
  tensor * list = NULL;
  while (v->x.i >= 0) {
    tensor t;
     {
      if (!(v->x.i >= 0)) qassert ("/home/pwachara/basilisk/src/common.h", 322, "v->x.i >= 0");
      t.x = *v++;
    } 
#line 321
{
      if (!(v->y.i >= 0)) qassert ("/home/pwachara/basilisk/src/common.h", 322, "v->x.i >= 0");
      t.y = *v++;
    } 
#line 321
{
      if (!(v->z.i >= 0)) qassert ("/home/pwachara/basilisk/src/common.h", 322, "v->x.i >= 0");
      t.z = *v++;
    }
    list = tensors_append (list, t);
  }
  return list;
}

static inline bool is_vertex_scalar (scalar s)
{
  
    if (_attribute[s.i].d.x != -1)
      return false;    
#line 333
if (_attribute[s.i].d.y != -1)
      return false;    
#line 333
if (_attribute[s.i].d.z != -1)
      return false;
  return true;
}

scalar * all = NULL;
scalar * baseblock = NULL;



scalar (* init_scalar) (scalar, const char *);
scalar (* init_vertex_scalar) (scalar, const char *);
vector (* init_vector) (vector, const char *);
vector (* init_face_vector) (vector, const char *);
tensor (* init_tensor) (tensor, const char *);
void (* scalar_clone) (scalar, scalar);






static double mpi_time = 0.;


typedef struct {
  clock_t c;
  struct timeval tv;
  double tm;
} timer;

timer timer_start (void)
{
  timer t;
  t.c = clock();
  gettimeofday (&t.tv, NULL);

  t.tm = mpi_time;

  return t;
}

double timer_elapsed (timer t)
{
  struct timeval tvend;
  gettimeofday (&tvend, NULL);
  return ((tvend.tv_sec - t.tv.tv_sec) +
   (tvend.tv_usec - t.tv.tv_usec)/1e6);
}



const vector zerof = {{_NVARMAX+4},{_NVARMAX+5},{_NVARMAX+6}};
const vector unityf = {{_NVARMAX+7},{_NVARMAX+8},{_NVARMAX+9}};
const scalar unity = {_NVARMAX+10};
const scalar zeroc = {_NVARMAX+11};



        vector fm = {{_NVARMAX+12},{_NVARMAX+13},{_NVARMAX+14}};
        scalar cm = {_NVARMAX+15};
#line 407 "/home/pwachara/basilisk/src/common.h"
void * matrix_new (int n, int p, size_t size)
{
  void ** m = ((void * *) pmalloc ((n)*sizeof(void *),__func__,__FILE__,__LINE__));
  char * a = ((char *) pmalloc ((n*p*size)*sizeof(char),__func__,__FILE__,__LINE__));
  for (int i = 0; i < n; i++)
    m[i] = a + i*p*size;
  return m;
}

double matrix_inverse (double ** m, int n, double pivmin)
{
  int indxc[n], indxr[n], ipiv[n];
  int i, icol = 0, irow = 0, j, k, l, ll;
  double big, dum, pivinv, minpiv = 1e30f;

  for (j = 0; j < n; j++)
    ipiv[j] = -1;

  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 0)
 for (k = 0; k < n; k++) {
   if (ipiv[k] == -1) {
     if (fabs (m[j][k]) >= big) {
       big = fabs (m[j][k]);
       irow = j;
       icol = k;
     }
   }
 }
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < n; l++)
 do { double _tmp_ = m[irow][l]; m[irow][l] = m[icol][l]; m[icol][l] = _tmp_; } while(false);
    indxr[i] = irow;
    indxc[i] = icol;
    if (fabs (m[icol][icol]) <= pivmin)
      return 0.;
    if (fabs (m[icol][icol]) < minpiv)
      minpiv = fabs (m[icol][icol]);
    pivinv = 1.0/m[icol][icol];
    m[icol][icol] = 1.0;
    for (l = 0; l < n; l++) m[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
 dum = m[ll][icol];
 m[ll][icol] = 0.0;
 for (l = 0; l < n; l++)
   m[ll][l] -= m[icol][l]*dum;
      }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
 do { double _tmp_ = m[k][indxr[l]]; m[k][indxr[l]] = m[k][indxc[l]]; m[k][indxc[l]] = _tmp_; } while(false);
  }
  return minpiv;
}

void matrix_free (void * m)
{
  pfree (((void **) m)[0],__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
}



typedef void (* free_solver_func) (void);

static Array * free_solver_funcs = NULL;

void free_solver_func_add (free_solver_func func)
{
  if (!free_solver_funcs)
    free_solver_funcs = array_new();
  array_append (free_solver_funcs, &func, sizeof(free_solver_func));
}



static char * display_defaults = NULL;

static void free_display_defaults() {
  pfree (display_defaults,__func__,__FILE__,__LINE__);
}

void display (const char * commands, bool overwrite)
{
  if (display_defaults == NULL)
    free_solver_func_add (free_display_defaults);
  if (overwrite) {
    pfree (display_defaults,__func__,__FILE__,__LINE__);
    display_defaults = pmalloc (strlen(commands) + 2,__func__,__FILE__,__LINE__);
    strcpy (display_defaults, "@");
    strcat (display_defaults, commands);
  }
  else {
    if (!display_defaults)
      display_defaults = pstrdup ("@",__func__,__FILE__,__LINE__);
    display_defaults =
      prealloc (display_defaults,
        strlen(display_defaults) + strlen(commands) + 1,__func__,__FILE__,__LINE__);
    strcat (display_defaults, commands);
  }
}



typedef struct {
  double x;

  double y;


  double z;

} _coord;



typedef struct {
  float r, g, b, a;
} vec4;
#line 540 "/home/pwachara/basilisk/src/common.h"
typedef struct {
  coord x, y, z;
} mat3;

OMP(omp declare reduction (+ : mat3 :
      omp_out.x.x += omp_in.x.x,
      omp_out.x.y += omp_in.x.y,
      omp_out.x.z += omp_in.x.z,
      omp_out.y.x += omp_in.y.x,
      omp_out.y.y += omp_in.y.y,
      omp_out.y.z += omp_in.y.z,
      omp_out.z.x += omp_in.z.x,
      omp_out.z.y += omp_in.z.y,
      omp_out.z.z += omp_in.z.z
      ))

typedef struct {
  uint32_t s;
} Adler32Hash;

static
inline void a32_hash_init (Adler32Hash * hash)
{
  hash->s = 0;
}

static
inline void a32_hash_add (Adler32Hash * hash, const void * data, size_t size)
{
  const uint8_t * buffer = (const uint8_t*) data;
  for (size_t n = 0; n < size; n++, buffer++)
    hash->s = *buffer + (hash->s << 6) + (hash->s << 16) - hash->s;
}

static
inline uint32_t a32_hash (const Adler32Hash * hash)
{
  return hash->s;
}
#line 14 "BaseFlow-cpp.c"
#line 1 "BaseFlow.c"
#line 1 "grid/multigrid3D.h"
#line 1 "/home/pwachara/basilisk/src/grid/multigrid3D.h"

#line 1 "grid/multigrid.h"
#line 1 "/home/pwachara/basilisk/src/grid/multigrid.h"



typedef double real;
#line 29 "/home/pwachara/basilisk/src/grid/multigrid.h"
int Dimensions_scale = 1;


typedef struct {
  Grid g;
  char * d;
  size_t * shift;
} Multigrid;

struct _Point {
  int i;

  int j;


  int k;

  int level;





  struct { int x, y, z; } n;





  #define _BLOCK_INDEX

};
static Point last_point;
#line 85 "/home/pwachara/basilisk/src/grid/multigrid.h"
#undef val
#define val(a,l,m,o) (((real *)((Multigrid *)grid)->d)[point.k + (o) +\
       (((size_t)(1 << point.level)) + 2*2)*\
       (point.j + (m) +\
        (point.i + (l))*(((size_t)(1 << point.level)) + 2*2)) +\
       ((Multigrid *)grid)->shift[point.level] +\
       _index(a,0)*((Multigrid *)grid)->shift[depth() + 1]])\

#line 92

#line 121 "/home/pwachara/basilisk/src/grid/multigrid.h"
#define allocated(a,l,m) (point.i+(a) >= 0 &&\
         point.i+(a) < ((size_t)(1 << point.level)) + 2*2 &&\
         point.j+(l) >= 0 &&\
         point.j+(l) < ((size_t)(1 << point.level)) + 2*2 &&\
         point.k+(m) >= 0 &&\
         point.k+(m) < ((size_t)(1 << point.level)) + 2*2)\

#line 127


#define allocated_child(a,l,m) (level < depth() &&\
         point.i > 0 && point.i <= ((size_t)(1 << point.level)) + 2 &&\
         point.j > 0 && point.j <= ((size_t)(1 << point.level)) + 2 &&\
         point.k > 0 && point.k <= ((size_t)(1 << point.level)) + 2)\

#line 133




#define depth() (grid->depth)
#line 186 "/home/pwachara/basilisk/src/grid/multigrid.h"
#define fine(a,l,m,o)\
(((real *)((Multigrid *)grid)->d)[2*point.k - 2 + (o) +\
   (((size_t)(1 << point.level))*2 + 2*2)*\
   (2*point.j - 2 + (m) +\
    (2*point.i - 2 + (l))*(((size_t)(1 << point.level))*2 + 2*2)) +\
   ((Multigrid *)grid)->shift[point.level + 1] +\
   _index(a,0)*((Multigrid *)grid)->shift[depth() + 1]])\

#line 193

#define coarse(a,l,m,o)\
(((real *)((Multigrid *)grid)->d)[(point.k + 2)/2 + (o) +\
   (((size_t)(1 << point.level))/2 + 2*2)*\
   ((point.j + 2)/2 + (m) +\
    ((point.i + 2)/2 + (l))*(((size_t)(1 << point.level))/2 + 2*2)) +\
   ((Multigrid *)grid)->shift[point.level - 1] +\
   _index(a,0)*((Multigrid *)grid)->shift[depth() + 1]])\

#line 201


void macro_POINT_VARIABLES (Point point)
{ 
  
#line 537 "/home/pwachara/basilisk/src/grid/config.h"
;  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
}
#line 244 "/home/pwachara/basilisk/src/grid/multigrid.h"
void macro_foreach_level (int l, char flags, Reduce reductions) {
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = l;
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
    for (_k = 2; _k < point.n.x + 2; _k++) {
      point.i = _k;

      for (point.j = 2; point.j < point.n.y + 2; point.j++)

 for (point.k = 2; point.k < point.n.z + 2; point.k++)


   ;
    }
  }
}

void macro_foreach (char flags, Reduce reductions) {
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
      for (_k = 2; _k < point.n.x + 2; _k++) {
 point.i = _k;

 for (point.j = 2; point.j < point.n.y + 2; point.j++)

   for (point.k = 2; point.k < point.n.z + 2; point.k++)


     ;
      }
  }
}

#define is_active(cell) (true)
#define is_leaf(cell) (point.level == depth())
#define is_local(cell) (true)
#define leaf 2
#define refine_cell(...) do {\
  fprintf (stderr, "grid depths do not match. Aborting.\n");\
  if (!(0)) qassert ("/home/pwachara/basilisk/src/grid/multigrid.h", 292, "0");\
} while (0)\

#line 294

#define tree ((Multigrid *)grid)
#line 1 "grid/foreach_cell.h"
#line 1 "/home/pwachara/basilisk/src/grid/foreach_cell.h"
#line 66 "/home/pwachara/basilisk/src/grid/foreach_cell.h"
void macro_foreach_cell_root (Point root)
{
  {
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point = {0};





    int kg = 0; NOT_UNUSED(kg);
    struct { int l, i, j, k, stage; } stack[20];

    int _s = -1;
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].k = root.k; stack[_s].stage = 0; };
    while (_s >= 0) {
      int stage;
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; stage = stack[_s].stage; _s--; };
      if (!allocated (0,0,0))
 continue;
      switch (stage) {
      case 0: {

 ;

 if (point.level < grid->depth) {
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 1; };
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };
 }
 break;
      }
#line 106 "/home/pwachara/basilisk/src/grid/foreach_cell.h"
      case 1: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 2; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;
      case 2: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 3; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;
      case 3: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 4; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;
      case 4: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 5; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;
      case 5: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 6; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;
      case 6: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 7; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;
      case 7: { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;

      }
    }
  }
}

void macro_foreach_cell()
{
  {





    Point root = {2,2,2,0};
#line 67
{
  {
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point = {0};





    int kg = 0; NOT_UNUSED(kg);
    struct { int l, i, j, k, stage; } stack[20];

    int _s = -1;
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].k = root.k; stack[_s].stage = 0; };
    while (_s >= 0) {
      int stage;
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; stage = stack[_s].stage; _s--; };
      if (!allocated (0,0,0))
 continue;
      switch (stage) {
      case 0: {
      
#line 136
{ 
  
#line 537 "/home/pwachara/basilisk/src/grid/config.h"
;  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 136 "/home/pwachara/basilisk/src/grid/foreach_cell.h"
;} 
#line 91
if (point.level < grid->depth) {
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 1; };
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };
 }
 break;
      }
#line 106 "/home/pwachara/basilisk/src/grid/foreach_cell.h"
      case 1: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 2; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;
      case 2: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 3; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;
      case 3: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 4; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;
      case 4: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 5; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;
      case 5: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 6; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;
      case 6: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 7; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;
      case 7: { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;

      }
    }
  }
}
  
#line 137
}
}

void macro_foreach_cell_all() {
  {
    Point root = {0};
    for (root.i = 2*Period.x; root.i <= 2*(2 - Period.x); root.i++)

      for (root.j = 2*Period.y; root.j <= 2*(2 - Period.y); root.j++)


 for (root.k = 2*Period.z; root.k <= 2*(2 - Period.z); root.k++) 
#line 22 "/home/pwachara/basilisk/src/grid/config.h"
{

#line 67 "/home/pwachara/basilisk/src/grid/foreach_cell.h"
{
  {
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point = {0};





    int kg = 0; NOT_UNUSED(kg);
    struct { int l, i, j, k, stage; } stack[20];

    int _s = -1;
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].k = root.k; stack[_s].stage = 0; };
    while (_s >= 0) {
      int stage;
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; stage = stack[_s].stage; _s--; };
      if (!allocated (0,0,0))
 continue;
      switch (stage) {
      case 0: {
     
#line 151
{ 
  
#line 537 "/home/pwachara/basilisk/src/grid/config.h"
;  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 151 "/home/pwachara/basilisk/src/grid/foreach_cell.h"
;} 
#line 91
if (point.level < grid->depth) {
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 1; };
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };
 }
 break;
      }
#line 106 "/home/pwachara/basilisk/src/grid/foreach_cell.h"
      case 1: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 2; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;
      case 2: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 3; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;
      case 3: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 4; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;
      case 4: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 5; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;
      case 5: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 6; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;
      case 6: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 7; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;
      case 7: { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;

      }
    }
  }
}
#line 22 "/home/pwachara/basilisk/src/grid/config.h"
}
  
#line 152 "/home/pwachara/basilisk/src/grid/foreach_cell.h"
}
}

void macro_foreach_cell_post_root (bool condition, Point root)
{
  {
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point = {0};





    int kg = 0; NOT_UNUSED(kg);
    struct { int l, i, j, k, stage; } stack[20];

    int _s = -1;
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].k = root.k; stack[_s].stage = 0; };
    while (_s >= 0) {
      int stage;
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; stage = stack[_s].stage; _s--; };
      if (!allocated (0,0,0))
 continue;
      switch (stage) {
      case 0: {
 if (point.level == grid->depth) {
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 8; };
 }
 else {
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 1; };
   if (condition)
     { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };
 }
 break;
      }
#line 210 "/home/pwachara/basilisk/src/grid/foreach_cell.h"
      case 1:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 2; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };
 break;
      case 2:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 3; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };
 break;
      case 3:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 4; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };
 break;
      case 4:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 5; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };
 break;
      case 5:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 6; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };
 break;
      case 6:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 7; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };
 break;
      case 7:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 8; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };
 break;

      default:
 ;

      }
    }
  }
}

void macro_foreach_cell_post (bool condition)
{
  {





    Point root = {2,2,2,0};
#line 156
{
  {
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point = {0};





    int kg = 0; NOT_UNUSED(kg);
    struct { int l, i, j, k, stage; } stack[20];

    int _s = -1;
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].k = root.k; stack[_s].stage = 0; };
    while (_s >= 0) {
      int stage;
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; stage = stack[_s].stage; _s--; };
      if (!allocated (0,0,0))
 continue;
      switch (stage) {
      case 0: {
 if (point.level == grid->depth) {
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 8; };
 }
 else {
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 1; };
   if (condition)
     { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };
 }
 break;
      }
#line 210 "/home/pwachara/basilisk/src/grid/foreach_cell.h"
      case 1:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 2; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };
 break;
      case 2:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 3; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };
 break;
      case 3:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 4; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };
 break;
      case 4:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 5; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };
 break;
      case 5:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 6; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };
 break;
      case 6:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 7; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };
 break;
      case 7:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 8; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };
 break;

      default:
      
#line 265
{ 
  
#line 537 "/home/pwachara/basilisk/src/grid/config.h"
;  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 265 "/home/pwachara/basilisk/src/grid/foreach_cell.h"
;}      
#line 249
}
    }
  }
}
  
#line 266
}
}

void macro_foreach_cell_post_all (bool condition)
{
  {
    Point root = {0};
    for (root.i = 0; root.i <= 2*2; root.i++)

      for (root.j = 0; root.j <= 2*2; root.j++)


 for (root.k = 0; root.k <= 2*2; root.k++) 
#line 22 "/home/pwachara/basilisk/src/grid/config.h"
{

#line 156 "/home/pwachara/basilisk/src/grid/foreach_cell.h"
{
  {
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point = {0};





    int kg = 0; NOT_UNUSED(kg);
    struct { int l, i, j, k, stage; } stack[20];

    int _s = -1;
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].k = root.k; stack[_s].stage = 0; };
    while (_s >= 0) {
      int stage;
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; stage = stack[_s].stage; _s--; };
      if (!allocated (0,0,0))
 continue;
      switch (stage) {
      case 0: {
 if (point.level == grid->depth) {
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 8; };
 }
 else {
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 1; };
   if (condition)
     { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };
 }
 break;
      }
#line 210 "/home/pwachara/basilisk/src/grid/foreach_cell.h"
      case 1:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 2; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };
 break;
      case 2:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 3; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };
 break;
      case 3:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 4; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };
 break;
      case 4:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 5; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };
 break;
      case 5:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 6; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };
 break;
      case 6:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 7; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };
 break;
      case 7:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 8; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };
 break;

      default:
     
#line 281
{ 
  
#line 537 "/home/pwachara/basilisk/src/grid/config.h"
;  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 281 "/home/pwachara/basilisk/src/grid/foreach_cell.h"
;}      
#line 249
}
    }
  }
}
#line 22 "/home/pwachara/basilisk/src/grid/config.h"
}
  
#line 282 "/home/pwachara/basilisk/src/grid/foreach_cell.h"
}
}

void macro_foreach_leaf()
{
#line 126
{
  {





    Point root = {2,2,2,0};
#line 67
{
  {
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point = {0};





    int kg = 0; NOT_UNUSED(kg);
    struct { int l, i, j, k, stage; } stack[20];

    int _s = -1;
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].k = root.k; stack[_s].stage = 0; };
    while (_s >= 0) {
      int stage;
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; stage = stack[_s].stage; _s--; };
      if (!allocated (0,0,0))
 continue;
      switch (stage) {
      case 0: {
    
#line 288
{ 
  
#line 537 "/home/pwachara/basilisk/src/grid/config.h"
;  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 288 "/home/pwachara/basilisk/src/grid/foreach_cell.h"
if (is_leaf (cell)) {
      if (is_active(cell) && is_local(cell))
 ;
      continue;
    }} 
#line 91
if (point.level < grid->depth) {
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 1; };
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };
 }
 break;
      }
#line 106 "/home/pwachara/basilisk/src/grid/foreach_cell.h"
      case 1: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 2; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;
      case 2: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 3; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;
      case 3: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 4; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;
      case 4: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 5; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;
      case 5: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 6; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;
      case 6: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 7; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;
      case 7: { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;

      }
    }
  }
}
  
#line 137
}
}

#line 293
}
#line 297 "/home/pwachara/basilisk/src/grid/multigrid.h"

void macro_foreach_face_generic (char flags, Reduce reductions,
    const char * order)
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
      for (_k = 2; _k <= point.n.x + 2; _k++) {
 point.i = _k;

 for (point.j = 2; point.j <= point.n.y + 2; point.j++)

   for (point.k = 2; point.k <= point.n.z + 2; point.k++)


     ;
      }
  }
}

#define is_coarse() (point.level < depth())
#line 385 "/home/pwachara/basilisk/src/grid/multigrid.h"
void macro_is_face_x (Point p) {
  if (p.j < p.n.y + 2 && p.k < p.n.z + 2) {
    int ig = -1; NOT_UNUSED(ig);
    ;
  }
}

void macro_is_face_y (Point p) {
  if (p.i < p.n.x + 2 && p.k < p.n.z + 2) {
    int jg = -1; NOT_UNUSED(jg);
    ;
  }
}

void macro_is_face_z (Point p) {
  if (p.i < p.n.x + 2 && p.j < p.n.y + 2) {
    int kg = -1; NOT_UNUSED(kg);
    ;
  }
}

void macro_foreach_child (Point point)
{
  {
    int _i = 2*point.i - 2;
    int _j = 2*point.j - 2;
    int _k = 2*point.k - 2;
    point.level++;
    point.n.x *= 2, point.n.y *= 2, point.n.z *= 2;
    for (int _l = 0; _l < 2; _l++)
      for (int _m = 0; _m < 2; _m++)
 for (int _n = 0; _n < 2; _n++) {
   point.i = _i + _l; point.j = _j + _m; point.k = _k + _n; 
  
   
#line 537 "/home/pwachara/basilisk/src/grid/config.h"
;  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
   
#line 419
;
 }
    point.i = (_i + 2)/2;
    point.j = (_j + 2)/2;
    point.k = (_k + 2)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2, point.n.z /= 2;
  }
}


#if TRASH
# undef trash
# define trash(list) reset(list, undefined)
#endif

#line 1 "grid/neighbors.h"
#line 1 "/home/pwachara/basilisk/src/grid/neighbors.h"
#line 33 "/home/pwachara/basilisk/src/grid/neighbors.h"
void macro_foreach_neighbor (int _s,
    Point point) {
  {
    const int _nn = _s;
    const int _i = point.i, _j = point.j, _k = point.k;
    for (int _l = - _nn; _l <= _nn; _l++) {
      point.i = _i + _l;
      for (int _m = - _nn; _m <= _nn; _m++) {
 point.j = _j + _m;
 for (int _n = - _nn; _n <= _nn; _n++) {
   point.k = _k + _n; 
  
   
#line 537 "/home/pwachara/basilisk/src/grid/config.h"
;  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;   
#line 45 "/home/pwachara/basilisk/src/grid/neighbors.h"
;
 }
      }
    }
    point.i = _i; point.j = _j; point.k = _k;
  }
}
#line 436 "/home/pwachara/basilisk/src/grid/multigrid.h"

void reset (void * alist, double val)
{
  scalar * list = (scalar *) alist;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s))
      for (int b = 0; b < _attribute[s.i].block; b++) {
 real * data = (real *) ((Multigrid *)grid)->d;
 data += (s.i + b)*((Multigrid *)grid)->shift[depth() + 1];
 for (size_t i = 0; i < ((Multigrid *)grid)->shift[depth() + 1]; i++, data++)
   *data = val;
      }}}
}
#line 519 "/home/pwachara/basilisk/src/grid/multigrid.h"
void macro_foreach_boundary_dir (int l, int d, Reduce reductions) {
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = l < 0 ? depth() : l;
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int * _i = &point.j, * _j = &point.k;
    int _n[2] = { point.n.y, point.n.z };
    if (d == left) {
      point.i = 2;
      ig = -1;
    }
    else if (d == right) {
      point.i = point.n.x + 2 - 1;
      ig = 1;
    }
    else if (d == bottom) {
      point.j = 2;
      _i = &point.i, _n[0] = point.n.x;
      jg = -1;
    }
    else if (d == top) {
      point.j = point.n.y + 2 - 1;
      _i = &point.i, _n[0] = point.n.x;
      jg = 1;
    }
    else if (d == back) {
      point.k = 2;
      _i = &point.i; _j = &point.j;
      _n[0] = point.n.x, _n[1] = point.n.y;
      kg = -1;
    }
    else if (d == front) {
      point.k = point.n.z + 2 - 1;
      _i = &point.i; _j = &point.j;
      _n[0] = point.n.x, _n[1] = point.n.y;
      kg = 1;
    }
    int _l;
    OMP(omp for schedule(static))
      for (_l = 0; _l < _n[0] + 2*2; _l++) {
 *_i = _l;
 for (int _m = 0; _m < _n[1] + 2*2; _m++) {
   *_j = _m;
   ;
 }
      }
  }
}

#define neighbor(o,p,q)\
  ((Point){point.i+o, point.j+p, point.k+q, point.level, point.n _BLOCK_INDEX})\

#line 571

#define is_boundary(point) (point.i < 2 || point.i >= point.n.x + 2 ||\
    point.j < 2 || point.j >= point.n.y + 2 ||\
    point.k < 2 || point.k >= point.n.z + 2)\

#line 575




extern double (* default_scalar_bc[]) (Point, Point, scalar, bool *);
static double periodic_bc (Point point, Point neighbor, scalar s, bool * data);

void macro_foreach_boundary (int b, Reduce reductions)
{
  if (default_scalar_bc[b] != periodic_bc) 
#line 22 "/home/pwachara/basilisk/src/grid/config.h"
{ 
#line 519 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = (
#line 585
depth()
#line 523
) < 0 ? depth() : (
#line 585
depth()
#line 523
);
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int * _i = &point.j, * _j = &point.k;
    int _n[2] = { point.n.y, point.n.z };
    if (b == left) {
      point.i = 2;
      ig = -1;
    }
    else if (b == right) {
      point.i = point.n.x + 2 - 1;
      ig = 1;
    }
    else if (b == bottom) {
      point.j = 2;
      _i = &point.i, _n[0] = point.n.x;
      jg = -1;
    }
    else if (b == top) {
      point.j = point.n.y + 2 - 1;
      _i = &point.i, _n[0] = point.n.x;
      jg = 1;
    }
    else if (b == back) {
      point.k = 2;
      _i = &point.i; _j = &point.j;
      _n[0] = point.n.x, _n[1] = point.n.y;
      kg = -1;
    }
    else if (b == front) {
      point.k = point.n.z + 2 - 1;
      _i = &point.i; _j = &point.j;
      _n[0] = point.n.x, _n[1] = point.n.y;
      kg = 1;
    }
    int _l;
    OMP(omp for schedule(static))
      for (_l = 0; _l < _n[0] + 2*2; _l++) {
 *_i = _l;
 for (int _m = 0; _m < _n[1] + 2*2; _m++) {
   *_j = _m;
      
#line 586
{  
#line 537 "/home/pwachara/basilisk/src/grid/config.h"
;  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 586
if (!is_boundary(point))
 ;} 
#line 564
}
      }
  }
}
#line 22 "/home/pwachara/basilisk/src/grid/config.h"
}

#line 588 "/home/pwachara/basilisk/src/grid/multigrid.h"
}

#define neighborp(k,l,o) neighbor(k,l,o)

static void box_boundary_level (const Boundary * b, scalar * scalars, int l)
{
  disable_fpe (FE_DIVBYZERO|FE_INVALID);
  for (int d = 0; d < 2*3; d++)
    if (default_scalar_bc[d] == periodic_bc)
      {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 if (!is_constant(s) && _attribute[s.i].block > 0) {
   if (is_vertex_scalar (s))
     _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = NULL;
   else if (_attribute[s.i].face) {
     vector v = _attribute[s.i].v;
     _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
   }
 }}}
  for (int bghost = 1; bghost <= 2; bghost++)
    for (int d = 0; d < 2*3; d++) {

      scalar * list = NULL, * listb = NULL;
      {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 if (!is_constant(s) && _attribute[s.i].block > 0) {
   scalar sb = s;

   if (_attribute[s.i].v.x.i >= 0) {

     int j = 0;
     while ((&_attribute[s.i].v.x)[j].i != s.i) j++;
     sb = (&_attribute[s.i].v.x)[(j - d/2 + 3) % 3];
   }

   if (_attribute[sb.i].boundary[d] && _attribute[sb.i].boundary[d] != periodic_bc) {
     list = list_append (list, s);
     listb = list_append (listb, sb);
   }
 }}}

      if (list) { 
#line 519
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = l < 0 ? depth() : l;
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int * _i = &point.j, * _j = &point.k;
    int _n[2] = { point.n.y, point.n.z };
    if (d == left) {
      point.i = 2;
      ig = -1;
    }
    else if (d == right) {
      point.i = point.n.x + 2 - 1;
      ig = 1;
    }
    else if (d == bottom) {
      point.j = 2;
      _i = &point.i, _n[0] = point.n.x;
      jg = -1;
    }
    else if (d == top) {
      point.j = point.n.y + 2 - 1;
      _i = &point.i, _n[0] = point.n.x;
      jg = 1;
    }
    else if (d == back) {
      point.k = 2;
      _i = &point.i; _j = &point.j;
      _n[0] = point.n.x, _n[1] = point.n.y;
      kg = -1;
    }
    else if (d == front) {
      point.k = point.n.z + 2 - 1;
      _i = &point.i; _j = &point.j;
      _n[0] = point.n.x, _n[1] = point.n.y;
      kg = 1;
    }
    int _l;
    OMP(omp for schedule(static))
      for (_l = 0; _l < _n[0] + 2*2; _l++) {
 *_i = _l;
 for (int _m = 0; _m < _n[1] + 2*2; _m++) {
   *_j = _m; 
#line 628
{  
#line 537 "/home/pwachara/basilisk/src/grid/config.h"
;  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 628
{
   scalar s, sb;
   {scalar*_i0=listb;scalar*_i1= list;if(_i0)for(sb=*_i0,s=*_i1;_i0->i>= 0;sb=*++_i0,s=*++_i1){ {
     if ((_attribute[s.i].face && sb.i == _attribute[s.i].v.x.i) || is_vertex_scalar (s)) {

       if (bghost == 1)
 
    val(s,(ig + 1)/2,(jg + 1)/2,(kg + 1)/2) =
    _attribute[sb.i].boundary[d] (point, neighborp(ig,jg,kg), s, NULL);
     }
     else

      
  val(s,bghost*ig,bghost*jg,bghost*kg) =
  _attribute[sb.i].boundary[d] (neighborp((1 - bghost)*ig,
       (1 - bghost)*jg,
       (1 - bghost)*kg),
    neighborp(bghost*ig,bghost*jg,bghost*kg),
    s, NULL);
   }}}
 }} 
#line 564
}
      }
  }
}
 
#line 649
pfree (list,__func__,__FILE__,__LINE__);
 pfree (listb,__func__,__FILE__,__LINE__);
      }
    }
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
#line 786 "/home/pwachara/basilisk/src/grid/multigrid.h"
void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  Multigrid * m = ((Multigrid *)grid);
  pfree (m->d,__func__,__FILE__,__LINE__);
  pfree (m->shift,__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
  grid = NULL;
}

int log_base2 (int n) {
  int m = n, r = 0;
  while (m > 1)
    m /= 2, r++;
  return (1 << r) < n ? r + 1 : r;
}

void init_grid (int n)
{
  free_grid();
  Multigrid * m = ((Multigrid *) pmalloc ((1)*sizeof(Multigrid),__func__,__FILE__,__LINE__));
  grid = (Grid *) m;
  grid->depth = grid->maxdepth = log_base2(n/Dimensions.x);
  N = (1 << grid->depth)*Dimensions.x;
#line 820 "/home/pwachara/basilisk/src/grid/multigrid.h"
  Boundary * b = ((Boundary *) pcalloc (1, sizeof(Boundary),__func__,__FILE__,__LINE__));
  b->level = box_boundary_level;
  add_boundary (b);

  Boundary * mpi_boundary_new();
  mpi_boundary_new();
#line 835 "/home/pwachara/basilisk/src/grid/multigrid.h"
  m->shift = ((size_t *) pmalloc ((depth() + 2)*sizeof(size_t),__func__,__FILE__,__LINE__));
  size_t totalsize = 0;
  for (int l = 0; l <= depth() + 1; l++) {
    m->shift[l] = totalsize;
    struct _Point point = { .level = l };
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    size_t size = 1;
    
      size *= point.n.x + 2*2;      size *= point.n.y + 2*2;      size *= point.n.z + 2*2;
    totalsize += size;
  }
  m->d = (char *) pmalloc(m->shift[depth() + 1]*datasize,__func__,__FILE__,__LINE__);
  reset (all, 0.);
}

void realloc_scalar (int size)
{
  Multigrid * p = ((Multigrid *)grid);
  datasize += size;
  p->d = (char *) prealloc (p->d, (p->shift[depth() + 1]*datasize)*sizeof(char),__func__,__FILE__,__LINE__);
}


int mpi_coords[3];
#line 867 "/home/pwachara/basilisk/src/grid/multigrid.h"
Point locate (double xp, double yp, double zp)
{
  Point point = {0};
  point.level = depth();
  point.n.x = point.n.y = point.n.z = 1 << point.level;
  point.level = -1;

  point.i = (xp - X0)/L0*point.n.x*Dimensions.x + 2 - mpi_coords[0]*point.n.x;
  if (point.i < 2 || point.i >= point.n.x + 2)
    return point;

  point.j = (yp - Y0)/L0*point.n.x*Dimensions.x + 2 - mpi_coords[1]*point.n.x;
  if (point.j < 2 || point.j >= point.n.y + 2)
    return point;


  point.k = (zp - Z0)/L0*point.n.x*Dimensions.x + 2 - mpi_coords[2]*point.n.x;
  if (point.k < 2 || point.k >= point.n.z + 2)
    return point;
#line 902 "/home/pwachara/basilisk/src/grid/multigrid.h"
  point.level = depth();
  return point;
}




#line 1 "grid/multigrid-common.h"
#line 1 "/home/pwachara/basilisk/src/grid/multigrid-common.h"


#line 1 "grid/variables.h"
#line 1 "/home/pwachara/basilisk/src/grid/variables.h"
void macro2_VARIABLES (Point point, int _ig, int _jg, int _kg)
{
  double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((_ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((_jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((_kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
}
#line 4 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
#line 1 "grid/cartesian-common.h"
#line 1 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
#line 1 "grid/events.h"
#line 1 "/home/pwachara/basilisk/src/grid/events.h"
typedef struct _Event Event;
typedef int (* Expr) (int *, double *, Event *);

struct _Event {
  int last, nexpr;
  int (* action) (const int, const double, Event *);
  Expr expr[3];
  int * arrayi;
  double * arrayt;
  char * file;
  int line;
  char * name;
  double t;
  int i, a;
  void * data;
  Event * next;
};

static Event * Events = NULL;

int iter = 0, inext = 0;
double t = 0, tnext = 0;
void init_events (void);
void event_register (Event event);
static void _init_solver (void);





static int END_EVENT = 1234567890;
static double TEND_EVENT = 1234567890;
static double TEPS = 1e-9;

static void event_error (Event * ev, const char * s)
{
  fprintf (ferr, "%s:%d: error: %s\n", ev->file, ev->line, s);
  exit (1);
}

static void init_event (Event * ev)
{
  if (ev->arrayi || ev->arrayt) {
    ev->i = -1; ev->t = - TEND_EVENT;
    if (ev->arrayi)
      ev->i = ev->arrayi[0];
    else
      ev->t = ev->arrayt[0];
    ev->a = 1;
    ev->expr[1] = NULL;
  }
  else {
    if (ev->nexpr > 0) {
      Expr init = NULL, cond = NULL, inc = NULL;
      for (int j = 0; j < ev->nexpr; j++) {
 int i = -123456; double t = - TEND_EVENT;
 (* ev->expr[j]) (&i, &t, ev);
 if (i == -123456 && t == - TEND_EVENT) {

   if (cond)
     event_error (ev, "events can only use a single condition");
   cond = ev->expr[j];
 }
 else {

   int i1 = i; double t1 = t;
   (* ev->expr[j]) (&i1, &t1, ev);
   if (i1 == i && t1 == t) {


     if (init)
       event_error (ev, "events can only use a single initialisation");
     init = ev->expr[j];
   }
   else {

     if (inc)
       event_error (ev, "events can only use a single increment");
     inc = ev->expr[j];
   }
 }
      }
      ev->expr[0] = init;
      ev->expr[1] = cond;
      ev->expr[2] = inc;
      ev->nexpr = 0;
    }
    ev->i = -1; ev->t = - TEND_EVENT;
    if (ev->expr[0]) {
      (* ev->expr[0]) (&ev->i, &ev->t, ev);
      if (ev->i == END_EVENT || ev->t == TEND_EVENT) {
 ev->i = END_EVENT; ev->t = - TEND_EVENT;
      }
    }
    else if (ev->expr[2]) {
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (ev->i != -1)
 ev->i = 0;
      if (ev->t != - TEND_EVENT)
 ev->t = 0;
    }
  }
}

enum { event_done, event_alive, event_stop };

static int event_finished (Event * ev)
{
  ev->i = -1; ev->t = - TEND_EVENT;
  return event_done;
}

void event_register (Event event) {
  if (!(Events)) qassert ("/home/pwachara/basilisk/src/grid/events.h", 114, "Events");
  if (!(!event.last)) qassert ("/home/pwachara/basilisk/src/grid/events.h", 115, "!event.last");
  int n = 0, parent = -1;
  for (Event * ev = Events; !ev->last; ev++) {
    if (!strcmp (event.name, ev->name)) {
      if (!(parent < 0)) qassert ("/home/pwachara/basilisk/src/grid/events.h", 119, "parent < 0");
      parent = n;
    }
    n++;
  }
  if (parent < 0) {
    Events = (Event *) prealloc (Events, (n + 2)*sizeof(Event),__func__,__FILE__,__LINE__);
    Events[n] = event;
    Events[n].next = NULL;
    Events[n + 1].last = true;
    init_event (&Events[n]);
  }
  else {
    Event * ev = ((Event *) pcalloc (1, sizeof(Event),__func__,__FILE__,__LINE__));
    *ev = Events[parent];
    Events[parent] = event;
    Events[parent].next = ev;
    init_event (&Events[parent]);
  }
}

static int event_cond (Event * ev, int i, double t)
{
  if (!ev->expr[1])
    return true;
  return (* ev->expr[1]) (&i, &t, ev);
}
#line 162 "/home/pwachara/basilisk/src/grid/events.h"
static bool overload_event() { return true; }

static int event_do (Event * ev, bool action)
{
  if ((iter > ev->i && t > ev->t) || !event_cond (ev, iter, t))
    return event_finished (ev);
  if (!overload_event() || iter == ev->i || fabs (t - ev->t) <= TEPS*t) {
    if (action) {
      bool finished = false;
      for (Event * e = ev; e; e = e->next) {



 if ((* e->action) (iter, t, e))
   finished = true;
      }
      if (finished) {
 event_finished (ev);
 return event_stop;
      }
    }
    if (ev->arrayi) {
      ev->i = ev->arrayi[ev->a++];
      if (ev->i < 0)
 return event_finished (ev);
    }
    if (ev->arrayt) {
      ev->t = ev->arrayt[ev->a++];
      if (ev->t < 0)
 return event_finished (ev);
    }
    else if (ev->expr[2]) {
      int i0 = ev->i;
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (i0 == -1 && ev->i != i0)
 ev->i += iter + 1;
      if (!event_cond (ev, iter + 1, ev->t))
 return event_finished (ev);
    }
    else if (ev->expr[0] && !ev->expr[1])
      return event_finished (ev);
  }
  return event_alive;
}

static void end_event_do (bool action)
{




  for (Event * ev = Events; !ev->last; ev++)
    if (ev->i == END_EVENT && action)
      for (Event * e = ev; e; e = e->next) {



 e->action (iter, t, e);
      }
}

int events (bool action)
{





  if (iter == 0)
    for (Event * ev = Events; !ev->last; ev++)
      init_event (ev);

  int cond = 0, cond1 = 0;
  inext = END_EVENT; tnext = 1e30f;
  for (Event * ev = Events; !ev->last && !cond; ev++)
    if (ev->i != END_EVENT &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond = 1;
  for (Event * ev = Events; !ev->last; ev++) {
    int status = event_do (ev, action);
    if (status == event_stop) {
      end_event_do (action);
      return 0;
    }
    if (status == event_alive && ev->i != END_EVENT &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond1 = 1;
    if (ev->t > t && ev->t < tnext)
      tnext = ev->t;
    if (ev->i > iter && ev->i < inext)
      inext = ev->i;
  }
  if (overload_event() && (!cond || cond1) && (tnext != 1e30f || inext != END_EVENT)) {
    inext = iter + 1;
    return 1;
  }
  end_event_do (action);
  return 0;
}

void event (const char * name)
{
  for (Event * ev = Events; !ev->last; ev++)
    if (!strcmp (ev->name, name))
      for (Event * e = ev; e; e = e->next) {



 (* e->action) (0, 0, e);
      }
}

double dtnext (double dt)
{
  if (tnext != 1e30f && tnext > t) {
    if (!(dt > 0.)) qassert ("/home/pwachara/basilisk/src/grid/events.h", 277, "dt > 0.");
    unsigned int n = (tnext - t)/dt;
    if (!(n < INT_MAX)) qassert ("/home/pwachara/basilisk/src/grid/events.h", 279, "n < INT_MAX");
    if (n == 0)
      dt = tnext - t;
    else {
      double dt1 = (tnext - t)/n;
      if (dt1 > dt*(1. + TEPS))
 dt = (tnext - t)/(n + 1);
      else if (dt1 < dt)
 dt = dt1;
      tnext = t + dt;
    }
  }
  else
    tnext = t + dt;
  return dt;
}

void init_solver()
{
  Events = pmalloc (sizeof (Event),__func__,__FILE__,__LINE__);
  Events[0].last = 1;
  _attribute = pcalloc (datasize/sizeof(real), sizeof (_Attributes),__func__,__FILE__,__LINE__);
  int n = datasize/sizeof(real);
  all = (scalar *) pmalloc (sizeof (scalar)*(n + 1),__func__,__FILE__,__LINE__);
  baseblock = (scalar *) pmalloc (sizeof (scalar)*(n + 1),__func__,__FILE__,__LINE__);
  for (int i = 0; i < n; i++)
    baseblock[i].i = all[i].i = i;
  baseblock[n].i = all[n].i = -1;




  mpi_init();





}
#line 2 "/home/pwachara/basilisk/src/grid/cartesian-common.h"

void (* debug) (Point);

#define _val_constant(a,k,l,m) ((const double) _constant[a.i -_NVARMAX])
#define val_diagonal(a,k,l,m) ((k) == 0 && (l) == 0 && (m) == 0)

#line 1 "grid/fpe.h"
#line 1 "/home/pwachara/basilisk/src/grid/fpe.h"


#include <signal.h>
#include <unistd.h>

static int gdb()
{
  if (last_point.level >= 0) {
    debug (last_point);
    fputc ('\n', ferr);
    fflush (ferr);
  }
  char command[80];
  sprintf (command, "exec xterm -e 'gdb -p %d' & xterm -e 'gnuplot plot -'",
    getpid());
  return system (command);
}

static void caught_abort (int sig)
{
  fprintf (ferr, "Caught signal %d (Aborted)\n", sig);
  gdb();
}

static void caught_fpe (int sig)
{
  fprintf (ferr, "Caught signal %d (Floating Point Exception)\n", sig);
  gdb();
  exit (1);
}

static void caught_segfault (int sig)
{
  fprintf (ferr, "Caught signal %d (Segmentation Fault)\n", sig);
  gdb();
  exit (2);
}

void catch_fpe (void)
{
  struct sigaction act;
  act.sa_handler = caught_fpe;
  sigemptyset (&act.sa_mask);
  act.sa_flags = 0;
  last_point.level = -1;
  sigaction (8, &act, NULL);
  act.sa_handler = caught_segfault;
  sigaction (11, &act, NULL);
  act.sa_handler = caught_abort;
  act.sa_flags = SA_RESETHAND;
  sigaction (6, &act, NULL);
}
#line 9 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
#line 1 "grid/stencils.h"
#line 1 "/home/pwachara/basilisk/src/grid/stencils.h"
#line 17 "/home/pwachara/basilisk/src/grid/stencils.h"










typedef struct _External External;

struct _External {
  char * name;
  void * pointer;
  int type;
  int nd;
  char reduct;
  char global;
  void * data;
  scalar s;
  External * externals, * next;
  int used;
};

typedef struct {
  const char * fname;
  int line;
  int first;
  int face;
  bool vertex;
  int parallel;
  scalar * listc;
  vectorl listf;
  scalar * dirty;
  void * data;
} ForeachData;







static inline bool scalar_is_dirty (scalar s)
{
  if (_attribute[s.i].dirty)
    return true;
  scalar * depends = _attribute[s.i].depends;
  {scalar*_i=(scalar*)( depends);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
    if (_attribute[d.i].dirty)
      return true;}}
  return false;
}




static inline bool scalar_depends_from (scalar a, scalar b)
{
  scalar * depends = _attribute[a.i].depends;
  {scalar*_i=(scalar*)( depends);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (s.i == b.i)
      return true;}}
  return false;
}







void boundary_internal (scalar * list, const char * fname, int line);
void (* boundary_face) (vectorl);







void check_stencil (ForeachData * loop)
{
  loop->listf = (vectorl){NULL};




  {scalar*_i=(scalar*)( baseblock);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    bool write = _attribute[s.i].output, read = _attribute[s.i].input;




    {





      if (read && scalar_is_dirty (s)) {





 if (_attribute[s.i].face) {
   if (_attribute[s.i].width > 0)
     loop->listc = list_append (loop->listc, s);
   else if (!write) {
     scalar sn = _attribute[s.i].v.x.i >= 0 ? _attribute[s.i].v.x : s;
     
       if (_attribute[s.i].v.x.i == s.i) {




  if (_attribute[sn.i].boundary[left] || _attribute[sn.i].boundary[right])
    loop->listc = list_append (loop->listc, s);
  else if (_attribute[s.i].dirty != 2)
    loop->listf.x = list_append (loop->listf.x, s);
       }       
#line 130
if (_attribute[s.i].v.y.i == s.i) {




  if (_attribute[sn.i].boundary[bottom] || _attribute[sn.i].boundary[top])
    loop->listc = list_append (loop->listc, s);
  else if (_attribute[s.i].dirty != 2)
    loop->listf.y = list_append (loop->listf.y, s);
       }       
#line 130
if (_attribute[s.i].v.z.i == s.i) {




  if (_attribute[sn.i].boundary[back] || _attribute[sn.i].boundary[front])
    loop->listc = list_append (loop->listc, s);
  else if (_attribute[s.i].dirty != 2)
    loop->listf.z = list_append (loop->listf.z, s);
       }
   }
 }





 else if (_attribute[s.i].width > 0)
   loop->listc = list_append (loop->listc, s);
      }





      if (write) {
 if (3 > 1 && !loop->vertex && loop->first && !_attribute[s.i].nowarning) {
   bool vertex = true;
   
     if (_attribute[s.i].d.x != -1)
       vertex = false;     
#line 159
if (_attribute[s.i].d.y != -1)
       vertex = false;     
#line 159
if (_attribute[s.i].d.z != -1)
       vertex = false;
   if (vertex)
     fprintf (ferr,
       "%s:%d: warning: vertex scalar '%s' should be assigned with"
       " a foreach_vertex() loop\n",
       loop->fname, loop->line, _attribute[s.i].name);
 }
 if (_attribute[s.i].face) {
   if (loop->face == 0 && loop->first && !_attribute[s.i].nowarning)
     fprintf (ferr,
       "%s:%d: warning: face vector '%s' should be assigned with"
       " a foreach_face() loop\n",
       loop->fname, loop->line, _attribute[s.i].name);
 }
 else if (loop->face) {
   if (_attribute[s.i].v.x.i < 0) {
     int d = 1, i = 0;
      {
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.x.i = s.i;
  _attribute[s.i].boundary[left] = _attribute[s.i].boundary[right] = NULL;





       }
       d *= 2, i++;
     } 
#line 177
{
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.y.i = s.i;
  _attribute[s.i].boundary[bottom] = _attribute[s.i].boundary[top] = NULL;





       }
       d *= 2, i++;
     } 
#line 177
{
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.z.i = s.i;
  _attribute[s.i].boundary[back] = _attribute[s.i].boundary[front] = NULL;





       }
       d *= 2, i++;
     }
     if (!_attribute[s.i].face && loop->first && !_attribute[s.i].nowarning)
       fprintf (ferr,
         "%s:%d: warning: scalar '%s' should be assigned with "
         "a foreach_face(x|y|z) loop\n",
         loop->fname, loop->line, _attribute[s.i].name);
   }
   else {
     char * name = NULL;
     if (_attribute[s.i].name) {
       name = pstrdup (_attribute[s.i].name,__func__,__FILE__,__LINE__);
       char * s = name + strlen(name) - 1;
       while (s != name && *s != '.') s--;
       if (s != name) *s = '\0';
     }
     struct { int x, y, z; } input, output;
     vector v = _attribute[s.i].v;

     
       input.x = _attribute[v.x.i].input, output.x = _attribute[v.x.i].output;       input.y = _attribute[v.y.i].input, output.y = _attribute[v.y.i].output;       input.z = _attribute[v.z.i].input, output.z = _attribute[v.z.i].output;

     init_face_vector (v, name);


     
       _attribute[v.x.i].input = input.x, _attribute[v.x.i].output = output.x;       _attribute[v.y.i].input = input.y, _attribute[v.y.i].output = output.y;       _attribute[v.z.i].input = input.z, _attribute[v.z.i].output = output.z;





     pfree (name,__func__,__FILE__,__LINE__);
   }
 }
 else if (loop->vertex) {
   bool vertex = true;
   
     if (_attribute[s.i].d.x != -1)
       vertex = false;     
#line 225
if (_attribute[s.i].d.y != -1)
       vertex = false;     
#line 225
if (_attribute[s.i].d.z != -1)
       vertex = false;
   if (!vertex) {
     char * name = NULL;
     if (_attribute[s.i].name) name = pstrdup (_attribute[s.i].name,__func__,__FILE__,__LINE__);
     init_vertex_scalar (s, name);
     
       _attribute[s.i].v.x.i = -1;       _attribute[s.i].v.y.i = -1;       _attribute[s.i].v.z.i = -1;




     pfree (name,__func__,__FILE__,__LINE__);
   }
 }





 loop->dirty = list_append (loop->dirty, s);
 {scalar*_i=(scalar*)( baseblock);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
   if (scalar_depends_from (d, s))
     loop->dirty = list_append (loop->dirty, d);}}
      }
    }
  }}}
}




void boundary_stencil (ForeachData * loop)
{
  bool flux = false;
  
    if (loop->listf.x)
      flux = true;    
#line 261
if (loop->listf.y)
      flux = true;    
#line 261
if (loop->listf.z)
      flux = true;
  if (flux) {
#line 276 "/home/pwachara/basilisk/src/grid/stencils.h"
    boundary_face (loop->listf);
    
      pfree (loop->listf.x,__func__,__FILE__,__LINE__), loop->listf.x = NULL;      pfree (loop->listf.y,__func__,__FILE__,__LINE__), loop->listf.y = NULL;      pfree (loop->listf.z,__func__,__FILE__,__LINE__), loop->listf.z = NULL;
  }




  if (loop->listc) {






    boundary_internal (loop->listc, loop->fname, loop->line);
    pfree (loop->listc,__func__,__FILE__,__LINE__), loop->listc = NULL;
  }





  if (loop->dirty) {






    {scalar*_i=(scalar*)( loop->dirty);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = true;}}
    pfree (loop->dirty,__func__,__FILE__,__LINE__), loop->dirty = NULL;
  }
}

void macro_foreach_stencil (char flags, Reduce reductions)
{
  {
    static int _first = 1.;
    ForeachData _loop = {
      .fname = __FILE__, .line = __LINE__, .first = _first
    };
    if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
 _attribute[s.i].input = _attribute[s.i].output = _attribute[s.i].nowarning = false;
 _attribute[s.i].width = 0;
      }
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point);

    ;

    check_stencil (&_loop);
    boundary_stencil (&_loop);
    _first = 0;
  }
}

void macro_foreach_vertex_stencil (char flags, Reduce reductions) {  
#line 314
{
    static int _first = 1.;
    ForeachData _loop = {
      .fname = "/home/pwachara/basilisk/src/grid/stencils.h", .line = 335, .first = _first
    };
    if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
 _attribute[s.i].input = _attribute[s.i].output = _attribute[s.i].nowarning = false;
 _attribute[s.i].width = 0;
      }
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point); 
#line 335
{
    _loop.vertex = true;
    ;
  }    
#line 328
check_stencil (&_loop);
    boundary_stencil (&_loop);
    _first = 0;
  }







}

void macro_foreach_face_stencil (char flags, Reduce reductions, const char * order) {  
#line 314
{
    static int _first = 1.;
    ForeachData _loop = {
      .fname = "/home/pwachara/basilisk/src/grid/stencils.h", .line = 342, .first = _first
    };
    if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
 _attribute[s.i].input = _attribute[s.i].output = _attribute[s.i].nowarning = false;
 _attribute[s.i].width = 0;
      }
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point);
    
#line 343
;    
#line 328
check_stencil (&_loop);
    boundary_stencil (&_loop);
    _first = 0;
  }

#line 344
}

void macro_foreach_level_stencil (int l, char flags, Reduce reductions) {
  if (0) {

    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point);
    ;
  }
}

void macro_foreach_coarse_level_stencil (int l, char flags, Reduce reductions) {  
#line 347
if (0) {

    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point);
    





;  
#line 352
}





}

void macro_foreach_level_or_leaf_stencil (int l, char flags, Reduce reductions) {  
#line 347
if (0) {

    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point);
    
#line 362
;  
#line 352
}










}

void macro_foreach_point_stencil (double xp, double yp, double zp, char flags, Reduce reductions)
{  
#line 314
{
    static int _first = 1.;
    ForeachData _loop = {
      .fname = "/home/pwachara/basilisk/src/grid/stencils.h", .line = 367, .first = _first
    };
    if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
 _attribute[s.i].input = _attribute[s.i].output = _attribute[s.i].nowarning = false;
 _attribute[s.i].width = 0;
      }
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point);
    
#line 368
;    
#line 328
check_stencil (&_loop);
    boundary_stencil (&_loop);
    _first = 0;
  }

#line 369
}

void macro_foreach_region_stencil (coord p, coord box[2], coord n, char flags, Reduce reductions)
{  
#line 314
{
    static int _first = 1.;
    ForeachData _loop = {
      .fname = "/home/pwachara/basilisk/src/grid/stencils.h", .line = 373, .first = _first
    };
    if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
 _attribute[s.i].input = _attribute[s.i].output = _attribute[s.i].nowarning = false;
 _attribute[s.i].width = 0;
      }
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point);
    
#line 374
;    
#line 328
check_stencil (&_loop);
    boundary_stencil (&_loop);
    _first = 0;
  }

#line 375
}

void macro__stencil_is_face_x (ForeachData l) { l.face |= (1 << 0); ; }
void macro__stencil_is_face_y (ForeachData l) { l.face |= (1 << 1); ; }
void macro__stencil_is_face_z (ForeachData l) { l.face |= (1 << 2); ; }

void stencil_val (Point p, scalar s, int i, int j, int k,
    const char * file, int line, bool overflow);
void stencil_val_a (Point p, scalar s, int i, int j, int k, bool input,
      const char * file, int line);

#define _stencil_val(a,_i,_j,_k)\
  stencil_val (point, a, _i, _j, _k, __FILE__, __LINE__, false)\

#line 508

#define _stencil_val_o(a,_i,_j,_k)\
  stencil_val (point, a, _i, _j, _k, __FILE__, __LINE__, true)\

#line 511

#define _stencil_val_a(a,_i,_j,_k)\
  stencil_val_a (point, a, _i, _j, _k, false, __FILE__, __LINE__)\

#line 514

#define _stencil_val_r(a,_i,_j,_k)\
  stencil_val_a (point, a, _i, _j, _k, true, __FILE__, __LINE__)\

#line 517


#define _stencil_fine(a,_i,_j,_k) _stencil_val(a,_i,_j,_k)
#define _stencil_fine(a,_i,_j,_k) _stencil_val(a,_i,_j,_k)
#define _stencil_fine_a(a,_i,_j,_k) _stencil_val_a(a,_i,_j,_k)
#define _stencil_fine_r(a,_i,_j,_k) _stencil_val_r(a,_i,_j,_k)

#define _stencil_coarse(a,_i,_j,_k) _stencil_val(a,_i,_j,_k)
#define _stencil_coarse_a(a,_i,_j,_k) _stencil_val_a(a,_i,_j,_k)
#define _stencil_coarse_r(a,_i,_j,_k) _stencil_val_r(a,_i,_j,_k)

#define r_assign(x)
#define _assign(x)

#define _stencil_neighbor(i,j,k)
#define _stencil_child(i,j,k)
#define _stencil_aparent(i,j,k)
#define _stencil_aparent_a(i,j,k)
#define _stencil_aparent_r(i,j,k)

#define _stencil_allocated(i,j,k) true

#define _stencil_neighborp(i,j,k) neighborp(i,j,k)

int _stencil_nop;
#define _stencil_val_higher_dimension (_stencil_nop = 1)
#define _stencil__val_constant(a,_i,_j,_k) (_stencil_nop = 1)
#define _stencil_val_diagonal(a,_i,_j,_k) (_stencil_nop = 1)

typedef void _stencil_undefined;

#define o_stencil -3
#line 10 "/home/pwachara/basilisk/src/grid/cartesian-common.h"

void macro_foreach_point (double _x, double _y, double _z,
        char flags, Reduce reductions)
{
  {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    coord _p = { _x, _y, _z };
    Point point = locate (_p.x, _p.y, _p.z);
    if (point.level >= 0)
      ;
  }
}

void macro_foreach_region (coord p, coord box[2], coord n,
         char flags, Reduce reductions)
{
  {
    if (n.x < 1) n.x = 1;
    if (n.y < 1) n.y = 1;
    if (n.z < 1) n.z = 1;

    for (int _i = 0; _i < (int) n.x; _i++) {
      p.x = box[0].x + (box[1].x - box[0].x)/n.x*(_i + 0.5);
      for (int _j = 0; _j < (int) n.y; _j++) {
 p.y = box[0].y + (box[1].y - box[0].y)/n.y*(_j + 0.5);
 for (int _k = 0; _k < (int) n.z; _k++) {
   p.z = box[0].z + (box[1].z - box[0].z)/n.z*(_k + 0.5);
   Point point = locate (p.x, p.y, p.z);
   if (point.level >= 0) {
     int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
     ;
   }
 }
      }
    }
  }
}




static inline
double dirichlet (double expr, Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 54 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
return 2.*expr - val(s,0,0,0);
}

static inline
double dirichlet_homogeneous (double expr, Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 60 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
return - val(s,0,0,0);
}

static inline
double dirichlet_face (double expr)
{
  return expr;
}

static inline
double dirichlet_face_homogeneous (double expr)
{
  return 0.;
}

static inline
double neumann (double expr, Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 78 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
return Delta*expr + val(s,0,0,0);
}

static inline
double neumann_homogeneous (double expr, Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 84 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
return val(s,0,0,0);
}
#line 145 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
static void init_block_scalar (scalar sb, const char * name, const char * ext,
          int n, int block)
{
  char bname[strlen(name) + strlen(ext) + 10];
  if (n == 0) {
    strcat (strcpy (bname, name), ext);
    _attribute[sb.i].block = block;
    baseblock = list_append (baseblock, sb);
  }
  else {
    sprintf (bname, "%s%d%s", name, n, ext);
    _attribute[sb.i].block = - n;
  }
  _attribute[sb.i].name = pstrdup (bname,__func__,__FILE__,__LINE__);
  all = list_append (all, sb);
}

#define interpreter_set_int(...)
#define interpreter_reset_scalar(...)

scalar alloc_block_scalar (const char * name, const char * ext, int block)
{
  interpreter_set_int (&block);
  int nvar = datasize/sizeof(real);

  scalar s = {0};
  while (s.i < nvar) {
    int n = 0;
    scalar sb = s;
    while (sb.i < nvar && n < block && _attribute[sb.i].freed)
      n++, sb.i++;
    if (n >= block) {
      memset (&_attribute[s.i], 0, block*sizeof (_Attributes));
      for (sb.i = s.i, n = 0; n < block; n++, sb.i++) {
 init_block_scalar (sb, name, ext, n, block);
 interpreter_reset_scalar (sb);
      }
      trash (((scalar []){s, {-1}}));
      return s;
    }
    s.i = sb.i + 1;
  }


  s = (scalar){nvar};
  if (!(nvar + block <= _NVARMAX)) qassert ("/home/pwachara/basilisk/src/grid/cartesian-common.h", 190, "nvar + block <= _NVARMAX");

  if (_attribute == NULL)
    _attribute = (_Attributes *) pcalloc (nvar + block + 1, sizeof (_Attributes),__func__,__FILE__,__LINE__);
  else
    _attribute = (_Attributes *)
      prealloc (_attribute, (nvar + block + 1)*sizeof (_Attributes),__func__,__FILE__,__LINE__);
  memset (&_attribute[nvar], 0, block*sizeof (_Attributes));
  for (int n = 0; n < block; n++, nvar++) {
    scalar sb = (scalar){nvar};
    init_block_scalar (sb, name, ext, n, block);
  }

  realloc_scalar (block*sizeof(real));
  trash (((scalar []){s, {-1}}));
  return s;
}

scalar new_block_scalar (const char * name, const char * ext, int block)
{
  scalar s = alloc_block_scalar (name, ext, block), sb;
  int n = 0;
  for (sb.i = s.i, n = 0; n < block; n++, sb.i++)
    init_scalar (sb, NULL);
  return s;
}

scalar new_scalar (const char * name)
{
  return init_scalar (alloc_block_scalar (name, "", 1), NULL);
}

scalar new_vertex_scalar (const char * name)
{
  return init_vertex_scalar (alloc_block_scalar (name, "", 1), NULL);
}

static vector alloc_block_vector (const char * name, int block)
{
  vector v;
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  
    v.x = alloc_block_scalar (name, ext.x, block);    v.y = alloc_block_scalar (name, ext.y, block);    v.z = alloc_block_scalar (name, ext.z, block);
  return v;
}

vector new_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_vector (v, NULL);
  return v;
}

vector new_face_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_face_vector (v, NULL);
  return v;
}

vector new_block_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    
      vb.x.i = v.x.i + i;      vb.y.i = v.y.i + i;      vb.z.i = v.z.i + i;
    init_vector (vb, NULL);
    
      _attribute[vb.x.i].block = - i;      _attribute[vb.y.i].block = - i;      _attribute[vb.z.i].block = - i;
  }
  
    _attribute[v.x.i].block = block;    _attribute[v.y.i].block = block;    _attribute[v.z.i].block = block;
  return v;
}

vector new_block_face_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    
      vb.x.i = v.x.i + i;      vb.y.i = v.y.i + i;      vb.z.i = v.z.i + i;
    init_face_vector (vb, NULL);
    
      _attribute[vb.x.i].block = - i;      _attribute[vb.y.i].block = - i;      _attribute[vb.z.i].block = - i;
  }
  
    _attribute[v.x.i].block = block;    _attribute[v.y.i].block = block;    _attribute[v.z.i].block = block;
  return v;
}

tensor new_tensor (const char * name)
{
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  tensor t;
   {
    strcat (strcpy (cname, name), ext.x);
    t.x = alloc_block_vector (cname, 1);
  } 
#line 287
{
    strcat (strcpy (cname, name), ext.y);
    t.y = alloc_block_vector (cname, 1);
  } 
#line 287
{
    strcat (strcpy (cname, name), ext.z);
    t.z = alloc_block_vector (cname, 1);
  }
  init_tensor (t, NULL);
  return t;
}

tensor new_symmetric_tensor (const char * name)
{
  struct { char * x, * y, * z; } ext = {".x.x", ".y.y", ".z.z"};
  tensor t;
  
    t.x.x = alloc_block_scalar (name, ext.x, 1);    t.y.y = alloc_block_scalar (name, ext.y, 1);    t.z.z = alloc_block_scalar (name, ext.z, 1);

    t.x.y = alloc_block_scalar (name, ".x.y", 1);
    t.y.x = t.x.y;


    t.x.z = alloc_block_scalar (name, ".x.z", 1);
    t.z.x = t.x.z;
    t.y.z = alloc_block_scalar (name, ".y.z", 1);
    t.z.y = t.y.z;




  init_tensor (t, NULL);
  return t;
}

static int nconst = 0;

void init_const_scalar (scalar s, const char * name, double val)
{
  if (s.i - _NVARMAX >= nconst) {
    _constant = (double *) prealloc (_constant, (s.i - _NVARMAX + 1)*sizeof(double),__func__,__FILE__,__LINE__);
    for (int i = nconst; i < s.i - _NVARMAX; i++)
      _constant[i] = 0.;
    nconst = s.i - _NVARMAX + 1;
  }
  _constant[s.i - _NVARMAX] = val;
}

scalar new_const_scalar (const char * name, int i, double val)
{
  scalar s = (scalar){i + _NVARMAX};
  init_const_scalar (s, name, val);
  return s;
}

void init_const_vector (vector v, const char * name, double * val)
{
  
    init_const_scalar (v.x, name, *val++);    init_const_scalar (v.y, name, *val++);    init_const_scalar (v.z, name, *val++);
}

vector new_const_vector (const char * name, int i, double * val)
{
  vector v;
  
    v.x.i = _NVARMAX + i++;    v.y.i = _NVARMAX + i++;    v.z.i = _NVARMAX + i++;
  init_const_vector (v, name, val);
  return v;
}

static void cartesian_scalar_clone (scalar clone, scalar src)
{
  char * cname = _attribute[clone.i].name;
  BoundaryFunc * boundary = _attribute[clone.i].boundary;
  BoundaryFunc * boundary_homogeneous = _attribute[clone.i].boundary_homogeneous;
  if (!(_attribute[src.i].block > 0 && _attribute[clone.i].block == _attribute[src.i].block)) qassert ("/home/pwachara/basilisk/src/grid/cartesian-common.h", 358, "src.block > 0 && clone.block == src.block");
  pfree (_attribute[clone.i].depends,__func__,__FILE__,__LINE__);
  _attribute[clone.i] = _attribute[src.i];
  _attribute[clone.i].name = cname;
  _attribute[clone.i].boundary = boundary;
  _attribute[clone.i].boundary_homogeneous = boundary_homogeneous;
  for (int i = 0; i < nboundary; i++) {
    _attribute[clone.i].boundary[i] = _attribute[src.i].boundary[i];
    _attribute[clone.i].boundary_homogeneous[i] = _attribute[src.i].boundary_homogeneous[i];
  }
  _attribute[clone.i].depends = list_copy (_attribute[src.i].depends);
}

scalar * list_clone (scalar * l)
{
  scalar * list = NULL;
  int nvar = datasize/sizeof(real), map[nvar];
  for (int i = 0; i < nvar; i++)
    map[i] = -1;
  {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    scalar c = _attribute[s.i].block > 1 ? new_block_scalar("c", "", _attribute[s.i].block) : new_scalar("c");
    scalar_clone (c, s);
    map[s.i] = c.i;
    list = list_append (list, c);
  }}}
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    {
      if (_attribute[s.i].v.x.i >= 0 && map[_attribute[s.i].v.x.i] >= 0)
 _attribute[s.i].v.x.i = map[_attribute[s.i].v.x.i];      
#line 385
if (_attribute[s.i].v.y.i >= 0 && map[_attribute[s.i].v.y.i] >= 0)
 _attribute[s.i].v.y.i = map[_attribute[s.i].v.y.i];      
#line 385
if (_attribute[s.i].v.z.i >= 0 && map[_attribute[s.i].v.z.i] >= 0)
 _attribute[s.i].v.z.i = map[_attribute[s.i].v.z.i];}}}
  return list;
}

void delete (scalar * list)
{
  if (all == NULL)
    return;

  {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){ {
    for (int i = 0; i < _attribute[f.i].block; i++) {
      scalar fb = {f.i + i};
      if (_attribute[f.i].delete)
 _attribute[f.i].delete (fb);
      pfree (_attribute[fb.i].name,__func__,__FILE__,__LINE__); _attribute[fb.i].name = NULL;
      pfree (_attribute[fb.i].boundary,__func__,__FILE__,__LINE__); _attribute[fb.i].boundary = NULL;
      pfree (_attribute[fb.i].boundary_homogeneous,__func__,__FILE__,__LINE__); _attribute[fb.i].boundary_homogeneous = NULL;
      pfree (_attribute[fb.i].depends,__func__,__FILE__,__LINE__); _attribute[fb.i].depends = NULL;
      _attribute[fb.i].freed = true;
    }
  }}}

  if (list == all) {
    all[0].i = -1;
    baseblock[0].i = -1;
    return;
  }

  trash (list);
  {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){ {
    if (_attribute[f.i].block > 0) {
      scalar * s;
      for (s = all; s->i >= 0 && s->i != f.i; s++);
      if (s->i == f.i) {
 for (; s[_attribute[f.i].block].i >= 0; s++)
   s[0] = s[_attribute[f.i].block];
 s->i = -1;
      }
      for (s = baseblock; s->i >= 0 && s->i != f.i; s++);
      if (s->i == f.i) {
 for (; s[1].i >= 0; s++)
   s[0] = s[1];
 s->i = -1;
      }
    }
  }}}
}

void free_solver()
{
  if (!(_val_higher_dimension == 0.)) qassert ("/home/pwachara/basilisk/src/grid/cartesian-common.h", 436, "_val_higher_dimension == 0.");

  if (free_solver_funcs) {
    free_solver_func * a = (free_solver_func *) free_solver_funcs->p;
    for (int i = 0; i < free_solver_funcs->len/sizeof(free_solver_func); i++)
      a[i] ();
    array_free (free_solver_funcs);
  }

  delete (all);
  pfree (all,__func__,__FILE__,__LINE__); all = NULL;
  pfree (baseblock,__func__,__FILE__,__LINE__); baseblock = NULL;
  for (Event * ev = Events; !ev->last; ev++) {
    Event * e = ev->next;
    while (e) {
      Event * next = e->next;
      pfree (e,__func__,__FILE__,__LINE__);
      e = next;
    }
  }

  pfree (Events,__func__,__FILE__,__LINE__); Events = NULL;
  pfree (_attribute,__func__,__FILE__,__LINE__); _attribute = NULL;
  pfree (_constant,__func__,__FILE__,__LINE__); _constant = NULL;




  free_grid();
  qpclose_all();
#if TRACE
  trace_off();
#endif
#if MTRACE
  pmuntrace();
#endif
#if _CADNA
  cadna_end();
#endif
}



void (* boundary_level) (scalar *, int l);
void (* boundary_face) (vectorl);




void boundary_flux (vector * list) __attribute__ ((deprecated));

void boundary_flux (vector * list)
{
  vectorl list1 = {NULL};
  {vector*_i=(vector*)( list);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
    {
      list1.x = list_append (list1.x, v.x);      list1.y = list_append (list1.y, v.y);      list1.z = list_append (list1.z, v.z);}}}
  boundary_face (list1);
  
    pfree (list1.x,__func__,__FILE__,__LINE__);    pfree (list1.y,__func__,__FILE__,__LINE__);    pfree (list1.z,__func__,__FILE__,__LINE__);
}

static scalar * list_add_depends (scalar * list, scalar s)
{
  {scalar*_i=(scalar*)( list);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){
    if (t.i == s.i)
      return list;}}
  scalar * list1 = list;
  {scalar*_i=(scalar*)( _attribute[s.i].depends);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
    if (_attribute[d.i].dirty)
      list1 = list_add_depends (list1, d);}}
  return list_append (list1, s);
}

     
void boundary_internal (scalar * list, const char * fname, int line)
{tracing("boundary_internal","/home/pwachara/basilisk/src/grid/cartesian-common.h",511);
  if (list == NULL)
    {end_tracing("boundary_internal","/home/pwachara/basilisk/src/grid/cartesian-common.h",514);return;}
  scalar * listc = NULL;
  vectorl listf = {NULL};
  bool flux = false;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].block > 0) {
      if (scalar_is_dirty (s)) {
 if (_attribute[s.i].face && _attribute[s.i].dirty != 2)
   {
     if (_attribute[s.i].v.x.i == s.i)
       listf.x = list_add (listf.x, s), flux = true;     
#line 523
if (_attribute[s.i].v.y.i == s.i)
       listf.y = list_add (listf.y, s), flux = true;     
#line 523
if (_attribute[s.i].v.z.i == s.i)
       listf.z = list_add (listf.z, s), flux = true;}
 if (!is_constant(cm) && _attribute[cm.i].dirty)
   listc = list_add_depends (listc, cm);
 if (_attribute[s.i].face != 2)
   listc = list_add_depends (listc, s);
      }




    }}}
  if (flux) {
    boundary_face (listf);
    
      pfree (listf.x,__func__,__FILE__,__LINE__);      pfree (listf.y,__func__,__FILE__,__LINE__);      pfree (listf.z,__func__,__FILE__,__LINE__);
  }
  if (listc) {






    boundary_level (listc, -1);
    {scalar*_i=(scalar*)( listc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = false;}}
    pfree (listc,__func__,__FILE__,__LINE__);
  }
end_tracing("boundary_internal","/home/pwachara/basilisk/src/grid/cartesian-common.h",552);}

void cartesian_boundary_level (scalar * list, int l)
{
  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, l); };
}

void cartesian_boundary_face (vectorl list)
{
  
    {scalar*_i=(scalar*)( list.x);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}    
#line 562
{scalar*_i=(scalar*)( list.y);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}    
#line 562
{scalar*_i=(scalar*)( list.z);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
}

static double symmetry (Point point, Point neighbor, scalar s, bool * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
  
#line 568 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
return val(s,0,0,0);
}

static double antisymmetry (Point point, Point neighbor, scalar s, bool * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
  
#line 573 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
return -val(s,0,0,0);
}

BoundaryFunc default_scalar_bc[] = {
  symmetry, symmetry, symmetry, symmetry, symmetry, symmetry
};

scalar cartesian_init_scalar (scalar s, const char * name)
{

  char * pname;
  if (name) {
    pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
    pname = pstrdup (name,__func__,__FILE__,__LINE__);
  }
  else
    pname = _attribute[s.i].name;
  int block = _attribute[s.i].block;
  BoundaryFunc * boundary = _attribute[s.i].boundary;
  BoundaryFunc * boundary_homogeneous = _attribute[s.i].boundary_homogeneous;
  _attribute[s.i].name = pname;
  if (block < 0)
    _attribute[s.i].block = block;
  else
    _attribute[s.i].block = block > 0 ? block : 1;

  _attribute[s.i].boundary = boundary ? boundary : (BoundaryFunc *) pmalloc (nboundary*sizeof (BoundaryFunc),__func__,__FILE__,__LINE__);
  _attribute[s.i].boundary_homogeneous = boundary_homogeneous ? boundary_homogeneous :
    (BoundaryFunc *) pmalloc (nboundary*sizeof (BoundaryFunc),__func__,__FILE__,__LINE__);
  for (int b = 0; b < nboundary; b++)
    _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] =
      b < 2*3 ? default_scalar_bc[b] : symmetry;
  _attribute[s.i].gradient = NULL;
   {
    _attribute[s.i].d.x = 0;
    _attribute[s.i].v.x.i = -1;
  } 
#line 606
{
    _attribute[s.i].d.y = 0;
    _attribute[s.i].v.y.i = -1;
  } 
#line 606
{
    _attribute[s.i].d.z = 0;
    _attribute[s.i].v.z.i = -1;
  }
  _attribute[s.i].face = false;
  return s;
}

scalar cartesian_init_vertex_scalar (scalar s, const char * name)
{
  cartesian_init_scalar (s, name);
  
    _attribute[s.i].d.x = -1;    _attribute[s.i].d.y = -1;    _attribute[s.i].d.z = -1;
  for (int d = 0; d < nboundary; d++)
    _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = NULL;
  return s;
}

BoundaryFunc default_vector_bc[] = {
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry
};

vector cartesian_init_vector (vector v, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
   {
    if (name) {
      char cname[strlen(name) + 3];
      strcat (strcpy (cname, name), ext.x);
      cartesian_init_scalar (v.x, cname);
    }
    else
      cartesian_init_scalar (v.x, NULL);
    _attribute[v.x.i].v = v;
  } 
#line 633
{
    if (name) {
      char cname[strlen(name) + 3];
      strcat (strcpy (cname, name), ext.y);
      cartesian_init_scalar (v.y, cname);
    }
    else
      cartesian_init_scalar (v.y, NULL);
    _attribute[v.y.i].v = v;
  } 
#line 633
{
    if (name) {
      char cname[strlen(name) + 3];
      strcat (strcpy (cname, name), ext.z);
      cartesian_init_scalar (v.z, cname);
    }
    else
      cartesian_init_scalar (v.z, NULL);
    _attribute[v.z.i].v = v;
  }

  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] =
      d < 2*3 ? default_vector_bc[d] : antisymmetry;
  return v;
}

vector cartesian_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
   {
    _attribute[v.x.i].d.x = -1;
    _attribute[v.x.i].face = true;
  } 
#line 653
{
    _attribute[v.y.i].d.y = -1;
    _attribute[v.y.i].face = true;
  } 
#line 653
{
    _attribute[v.z.i].d.z = -1;
    _attribute[v.z.i].face = true;
  }
  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
  return v;
}

tensor cartesian_init_tensor (tensor t, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
   {
    if (name) {
      char cname[strlen(name) + 3];
      strcat (strcpy (cname, name), ext.x);
      cartesian_init_vector (t.x, cname);
    }
    else
      cartesian_init_vector (t.x, NULL);
  } 
#line 665
{
    if (name) {
      char cname[strlen(name) + 3];
      strcat (strcpy (cname, name), ext.y);
      cartesian_init_vector (t.y, cname);
    }
    else
      cartesian_init_vector (t.y, NULL);
  } 
#line 665
{
    if (name) {
      char cname[strlen(name) + 3];
      strcat (strcpy (cname, name), ext.z);
      cartesian_init_vector (t.z, cname);
    }
    else
      cartesian_init_vector (t.z, NULL);
  }
#line 689 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
    if (!(false)) qassert ("/home/pwachara/basilisk/src/grid/cartesian-common.h", 689, "false");

  return t;
}

void output_cells (FILE * fp, coord c, double size)
{ 
#line 265 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
      for (_k = 2; _k < point.n.x + 2; _k++) {
 point.i = _k;

 for (point.j = 2; point.j < point.n.y + 2; point.j++)

   for (point.k = 2; point.k < point.n.z + 2; point.k++) 
#line 696 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 696 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
{
    bool inside = true;
    coord o = {x,y,z};
    
      if (inside && size > 0. &&
   (o.x > c.x + size || o.x < c.x - size))
 inside = false;      
#line 700
if (inside && size > 0. &&
   (o.y > c.y + size || o.y < c.y - size))
 inside = false;      
#line 700
if (inside && size > 0. &&
   (o.z > c.z + size || o.z < c.z - size))
 inside = false;
    if (inside) {
      Delta /= 2.;
#line 715 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
      for (int i = -1; i <= 1; i += 2) {
 fprintf (fp, "%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n\n",
   x - Delta, y - Delta, z + i*Delta,
   x - Delta, y + Delta, z + i*Delta,
   x + Delta, y + Delta, z + i*Delta,
   x + Delta, y - Delta, z + i*Delta,
   x - Delta, y - Delta, z + i*Delta);
 for (int j = -1; j <= 1; j += 2)
   fprintf (fp, "%g %g %g\n%g %g %g\n\n",
     x + i*Delta, y + j*Delta, z - Delta,
     x + i*Delta, y + j*Delta, z + Delta);
      }

    }
  }}      
#line 282 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
  }
}
  
#line 730 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
fflush (fp);
}
#line 740 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
static char * replace_ (const char * vname)
{
  char * name = pstrdup (vname,__func__,__FILE__,__LINE__), * c = name;
  while (*c != '\0') {
    if (*c == '.')
      *c = '_';
    c++;
  }
  return name;
}

static void debug_plot (FILE * fp, const char * name, const char * cells,
   const char * stencil)
{
  char * vname = replace_ (name);
  fprintf (fp,
    "  load 'debug.plot'\n"
    "  v=%s\n"







    "  splot '%s' w l lc 0, "
    "'%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 1"
           " title columnhead(4+4*v)",

    vname, cells, stencil);
  pfree (vname,__func__,__FILE__,__LINE__);
}

void cartesian_debug (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
  
#line 775 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
char name[80] = "cells";
  if (pid() > 0)
    sprintf (name, "cells-%d", pid());
  FILE * fp = fopen (name, "w");
  output_cells (fp, (coord){x,y,z}, 4.*Delta);
  fclose (fp);

  char stencil[80] = "stencil";
  if (pid() > 0)
    sprintf (stencil, "stencil-%d", pid());
  fp = fopen (stencil, "w");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){





    fprintf (fp, "x y z %s ", _attribute[v.i].name);}}

  fputc ('\n', fp);
#line 821 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
    for (int k = -2; k <= 2; k++)
      for (int l = -2; l <= 2; l++)
 for (int m = -2; m <= 2; m++) {
   {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){ {
     fprintf (fp, "%g %g %g ",
       x + k*Delta + _attribute[v.i].d.x*Delta/2.,
       y + l*Delta + _attribute[v.i].d.y*Delta/2.,
       z + m*Delta + _attribute[v.i].d.z*Delta/2.);
     if (allocated(k,l,m))
       fprintf (fp, "%g ", val(v,k,l,m));
     else
       fputs ("n/a ", fp);
   }}}
   fputc ('\n', fp);
 }

  fclose (fp);

  fp = fopen ("debug.plot", "w");
  fprintf (fp,
    "set term x11\n"
    "set size ratio -1\n"
    "set key outside\n");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    char * name = replace_ (_attribute[s.i].name);
    fprintf (fp, "%s = %d\n", name, s.i);
    pfree (name,__func__,__FILE__,__LINE__);
  }}}
  fclose (fp);

  fprintf (ferr, "Last point stencils can be displayed using (in gnuplot)\n");
  debug_plot (ferr, _attribute[0].name, name, stencil);
  fflush (ferr);

  fp = fopen ("plot", "w");
  debug_plot (fp, _attribute[0].name, name, stencil);
  fclose (fp);
}

void cartesian_methods()
{
  init_scalar = cartesian_init_scalar;
  init_vertex_scalar = cartesian_init_vertex_scalar;
  init_vector = cartesian_init_vector;
  init_face_vector = cartesian_init_face_vector;
  init_tensor = cartesian_init_tensor;
  boundary_level = cartesian_boundary_level;
  boundary_face = cartesian_boundary_face;
  scalar_clone = cartesian_scalar_clone;
  debug = cartesian_debug;
}

tensor init_symmetric_tensor (tensor t, const char * name)
{
  return init_tensor (t, name);
}

static double interpolate_linear (Point point, scalar v,
      double xp, double yp, double zp)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 896 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
  x = (xp - x)/Delta - _attribute[v.i].d.x/2.;
  y = (yp - y)/Delta - _attribute[v.i].d.y/2.;
  z = (zp - z)/Delta - _attribute[v.i].d.z/2.;
  int i = ( (int)(x > 0 ? 1 : -1)), j = ( (int)(y > 0 ? 1 : -1)), k = ( (int)(z > 0 ? 1 : -1));
  x = fabs(x); y = fabs(y); z = fabs(z);

  return (((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +
    (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y)*(1. - z) +
   ((val(v,0,0,k)*(1. - x) + val(v,i,0,k)*x)*(1. - y) +
    (val(v,0,j,k)*(1. - x) + val(v,i,j,k)*x)*y)*z);

}
#line 878
static void _stencil_interpolate_linear (Point point, scalar v,
_stencil_undefined * xp,_stencil_undefined * yp,_stencil_undefined * zp)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;                              
#line 896 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
        
        
        
  
          

_stencil_val(v,0,0,0);_stencil_val(v, o_stencil,0,0);
_stencil_val(v,0,o_stencil,0); _stencil_val(v,o_stencil,o_stencil,0);
_stencil_val(v,0,0,o_stencil); _stencil_val(v,o_stencil,0,o_stencil);
_stencil_val(v,0,o_stencil,o_stencil); _stencil_val(v,o_stencil,o_stencil,o_stencil);  
#line 902
return           
              
    
    ;

}

     
double interpolate (scalar v, double xp, double yp, double zp,
      bool linear)
{tracing("interpolate","/home/pwachara/basilisk/src/grid/cartesian-common.h",910);
  double val = 1e30f;  
#line 314 "/home/pwachara/basilisk/src/grid/stencils.h"
{
    static int _first = 1.;
    ForeachData _loop = {
      .fname = "/home/pwachara/basilisk/src/grid/cartesian-common.h", .line = 914, .first = _first
    };
    if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
 _attribute[s.i].input = _attribute[s.i].output = _attribute[s.i].nowarning = false;
 _attribute[s.i].width = 0;
      }
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point);
    
#line 915 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
{ _stencil_interpolate_linear (point, v, NULL, NULL, NULL); _stencil_val(v,0,0,0);    }    
#line 328 "/home/pwachara/basilisk/src/grid/stencils.h"
check_stencil (&_loop);
    boundary_stencil (&_loop);
    _first = 0;
  }
#line 13 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
{
  {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    coord _p = { xp, yp, zp };
    Point point = locate (_p.x, _p.y, _p.z);
    if (point.level >= 0)
    
#line 915
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 915 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
val = linear ? interpolate_linear (point, v, xp, yp, zp) : val(v,0,0,0);}  
#line 20
}
}
#line 915
{mpi_all_reduce_array(&val,MPI_DOUBLE,MPI_MIN,1);}
  {end_tracing("interpolate","/home/pwachara/basilisk/src/grid/cartesian-common.h",916);return val;}
end_tracing("interpolate","/home/pwachara/basilisk/src/grid/cartesian-common.h",917);}

     
void interpolate_array (scalar * list, coord * a, int n, double * v,
   bool linear)
{tracing("interpolate_array","/home/pwachara/basilisk/src/grid/cartesian-common.h",920);
  int len = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    len++;}}
  for (int i = 0; i < n; i++) {
    double * w = v;
#line 937 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      *(w++) = 1e30f;}}    
#line 314 "/home/pwachara/basilisk/src/grid/stencils.h"
{
    static int _first = 1.;
    ForeachData _loop = {
      .fname = "/home/pwachara/basilisk/src/grid/cartesian-common.h", .line = 939, .first = _first
    };
    if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
 _attribute[s.i].input = _attribute[s.i].output = _attribute[s.i].nowarning = false;
 _attribute[s.i].width = 0;
      }
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point); 
#line 939 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
{   
      
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 { _stencil_val(s,0,0,0); _stencil_interpolate_linear (point, s, NULL, NULL, NULL);    }}}
    }    
#line 328 "/home/pwachara/basilisk/src/grid/stencils.h"
check_stencil (&_loop);
    boundary_stencil (&_loop);
    _first = 0;
  }
#line 13 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
{
  {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    coord _p = { (
#line 939
a[i].x
#line 16
), (
#line 939
a[i].y
#line 16
), (
#line 939
a[i].z
#line 16
) };
    Point point = locate (_p.x, _p.y, _p.z);
    if (point.level >= 0) 
#line 939
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 939 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
{
      int j = 0;
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 v[j++] = !linear ? val(s,0,0,0) : interpolate_linear (point, s, a[i].x, a[i].y, a[i].z);}}
    }}  
#line 20
}
}
#line 943
{mpi_all_reduce_array(v,MPI_DOUBLE,MPI_MIN,len);}

    v = w;
  }
end_tracing("interpolate_array","/home/pwachara/basilisk/src/grid/cartesian-common.h",947);}



typedef int bid;

bid new_bid()
{
  int b = nboundary++;
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    _attribute[s.i].boundary = (BoundaryFunc *) prealloc (_attribute[s.i].boundary, nboundary*sizeof (BoundaryFunc),__func__,__FILE__,__LINE__);
    _attribute[s.i].boundary_homogeneous = (BoundaryFunc *)
      prealloc (_attribute[s.i].boundary_homogeneous, nboundary*sizeof (BoundaryFunc),__func__,__FILE__,__LINE__);
  }}}
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    if (_attribute[s.i].v.x.i < 0)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] = symmetry;
    else if (_attribute[s.i].v.x.i == s.i) {
      vector v = _attribute[s.i].v;
      
 _attribute[v.y.i].boundary[b] = _attribute[v.y.i].boundary_homogeneous[b] = symmetry; _attribute[v.z.i].boundary[b] = _attribute[v.z.i].boundary_homogeneous[b] = symmetry; _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] = symmetry;
      _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] =
 _attribute[v.x.i].face ? NULL : antisymmetry;
    }
  }}}
  return b;
}



static double periodic_bc (Point point, Point neighbor, scalar s, bool * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
  
#line 979 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
return val(s,0,0,0);
}

static void periodic_boundary (int d)
{

  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (is_vertex_scalar (s))
      _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = NULL;
    else
      _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = periodic_bc;}}

  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (_attribute[s.i].face) {
      vector v = _attribute[s.i].v;
      _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
    }}}

  default_scalar_bc[d] = periodic_bc;
  default_vector_bc[d] = periodic_bc;
}

void periodic (int dir)
{





    if (!(dir <= back)) qassert ("/home/pwachara/basilisk/src/grid/cartesian-common.h", 1008, "dir <= back");


  int c = dir/2;
  periodic_boundary (2*c);
  periodic_boundary (2*c + 1);
  (&Period.x)[c] = true;
}


double getvalue (Point point, scalar s, int i, int j, int k)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
  
#line 1020 "/home/pwachara/basilisk/src/grid/cartesian-common.h"
return val(s,i,j,k);
}

void default_stencil (Point p, scalar * list)
{
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    if (_attribute[s.i].v.x.i != -1) {
      vector v = _attribute[s.i].v;
      {scalar*_i=(scalar*)(((vector[]) {v,{{-1},{-1},{-1}}}));if(_i)for(scalar c=*_i;(&c)->i>=0;c=*++_i){
 _attribute[c.i].input = _attribute[c.i].output = _attribute[c.i].nowarning = true, _attribute[c.i].width = 2;}}
    }
    else
      _attribute[s.i].input = _attribute[s.i].output = _attribute[s.i].nowarning = true, _attribute[s.i].width = 2;
  }}}
}




static void write_stencil_index (int * index)
{
  fprintf (qstderr(), "[%d", index[0]);
  for (int d = 1; d < 3; d++)
    fprintf (qstderr(), ",%d", index[d]);
  fputs ("]", qstderr());
}

void stencil_val (Point p, scalar s, int i, int j, int k,
    const char * file, int line, bool overflow)
{
  if (is_constant(s) || s.i < 0)
    return;
  if (_attribute[s.i].block < 0)
    s.i += _attribute[s.i].block;
  if (!_attribute[s.i].name) {
    fprintf (ferr, "%s:%d: error: trying to access a deleted field\n",
      file, line);
    exit (1);
  }
  int index[] = {i, j, k};
  for (int d = 0; d < 3; d++) {
    if (index[d] == o_stencil)
      index[d] = 2;
    else
      index[d] += (&p.i)[d];
  }
  bool central = true;
  for (int d = 0; d < 3; d++) {
    if (!overflow && (index[d] > 2 || index[d] < - 2)) {
      fprintf (qstderr(), "%s:%d: error: stencil overflow: %s",
        file, line, _attribute[s.i].name);
      write_stencil_index (index);
      fprintf (qstderr(), "\n");
      fflush (qstderr());
      abort();
    }
    if (index[d] != 0)
      central = false;
  }
  if (central) {
    if (!_attribute[s.i].output)
      _attribute[s.i].input = true;
  }
  else {
    _attribute[s.i].input = true;
    int d = 0;
     {
      if ((!_attribute[s.i].face || _attribute[s.i].v.x.i != s.i) && abs(index[d]) > _attribute[s.i].width)
 _attribute[s.i].width = abs(index[d]);
      d++;
    } 
#line 1086
{
      if ((!_attribute[s.i].face || _attribute[s.i].v.y.i != s.i) && abs(index[d]) > _attribute[s.i].width)
 _attribute[s.i].width = abs(index[d]);
      d++;
    } 
#line 1086
{
      if ((!_attribute[s.i].face || _attribute[s.i].v.z.i != s.i) && abs(index[d]) > _attribute[s.i].width)
 _attribute[s.i].width = abs(index[d]);
      d++;
    }
  }
}

void stencil_val_a (Point p, scalar s, int i, int j, int k, bool input,
      const char * file, int line)
{
  if (is_constant(s) || s.i < 0) {
    fprintf (ferr, "%s:%d: error: trying to modify a%s field\n",
      file, line, s.i < 0 ? "n undefined" : " constant");
    exit (1);
  }
  if (_attribute[s.i].block < 0)
    s.i += _attribute[s.i].block;
  if (!_attribute[s.i].name) {
    fprintf (ferr, "%s:%d: error: trying to access a deleted field\n",
      file, line);
    exit (1);
  }
  int index[] = {i, j, k};
  for (int d = 0; d < 3; d++)
    index[d] += (&p.i)[d];
  for (int d = 0; d < 3; d++)
    if (index[d] != 0) {
      fprintf (qstderr(), "%s:%d: error: illegal write: %s",
        file, line, _attribute[s.i].name);
      write_stencil_index (index);
      fprintf (qstderr(), "\n");
      fflush (qstderr());
      abort();
    }
  if (input && !_attribute[s.i].output)
    _attribute[s.i].input = true;
  _attribute[s.i].output = true;
} 
#line 5 "/home/pwachara/basilisk/src/grid/multigrid-common.h"

void macro_foreach_level_or_leaf (int l, char flags, Reduce reductions)
{ 
#line 244 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = l;
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
    for (_k = 2; _k < point.n.x + 2; _k++) {
      point.i = _k;

      for (point.j = 2; point.j < point.n.y + 2; point.j++)

 for (point.k = 2; point.k < point.n.z + 2; point.k++)    
#line 9 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 9 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
;}
    
#line 261 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
  }
}
#line 10 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
} 

void macro_foreach_coarse_level (int l, char flags, Reduce reductions)
{ 
#line 244 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = l;
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
    for (_k = 2; _k < point.n.x + 2; _k++) {
      point.i = _k;

      for (point.j = 2; point.j < point.n.y + 2; point.j++)

 for (point.k = 2; point.k < point.n.z + 2; point.k++)    
#line 15 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 15 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
;}
    
#line 261 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
  }
}
#line 16 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
}










void (* restriction) (scalar *);

static inline void restriction_average (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 31 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
double sum = 0.;
  
  
#line 408 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
    int _i = 2*point.i - 2;
    int _j = 2*point.j - 2;
    int _k = 2*point.k - 2;
    point.level++;
    point.n.x *= 2, point.n.y *= 2, point.n.z *= 2;
    for (int _l = 0; _l < 2; _l++)
      for (int _m = 0; _m < 2; _m++)
 for (int _n = 0; _n < 2; _n++) {
   point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;   
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;    
#line 33 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
sum += val(s,0,0,0);
 
#line 420 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
    point.i = (_i + 2)/2;
    point.j = (_j + 2)/2;
    point.k = (_k + 2)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2, point.n.z /= 2;
  }  
#line 34 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
val(s,0,0,0) = sum/(1 << 3);
}
#line 29
static void _stencil_restriction_average (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
     
  
  
#line 408
{
    int _i = 2*point.i - 2;
    int _j = 2*point.j - 2;
    int _k = 2*point.k - 2;
    point.level++;
    point.n.x *= 2, point.n.y *= 2, point.n.z *= 2;
    for (int _l = 0; _l < 2; _l++)
      for (int _m = 0; _m < 2; _m++)
 for (int _n = 0; _n < 2; _n++) {
   point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;   
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;    
#line 33 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{ _stencil_val(s,0,0,0); }
 
#line 420 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
    point.i = (_i + 2)/2;
    point.j = (_j + 2)/2;
    point.k = (_k + 2)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2, point.n.z /= 2;
  }  
#line 34 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
_stencil_val_a(s,0,0,0);    
}

static inline void restriction_volume_average (Point point, scalar s)
{if(!is_constant(cm)){{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 39 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
double sum = 0.;
  
  
#line 408 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
    int _i = 2*point.i - 2;
    int _j = 2*point.j - 2;
    int _k = 2*point.k - 2;
    point.level++;
    point.n.x *= 2, point.n.y *= 2, point.n.z *= 2;
    for (int _l = 0; _l < 2; _l++)
      for (int _m = 0; _m < 2; _m++)
 for (int _n = 0; _n < 2; _n++) {
   point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;   
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;    
#line 41 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
sum += val(cm,0,0,0)*val(s,0,0,0);
 
#line 420 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
    point.i = (_i + 2)/2;
    point.j = (_j + 2)/2;
    point.k = (_k + 2)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2, point.n.z /= 2;
  }  
#line 42 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
val(s,0,0,0) = sum/(1 << 3)/(val(cm,0,0,0) + 1e-30);
}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
#line 38
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 39 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
double sum = 0.;
  
  
#line 408 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
    int _i = 2*point.i - 2;
    int _j = 2*point.j - 2;
    int _k = 2*point.k - 2;
    point.level++;
    point.n.x *= 2, point.n.y *= 2, point.n.z *= 2;
    for (int _l = 0; _l < 2; _l++)
      for (int _m = 0; _m < 2; _m++)
 for (int _n = 0; _n < 2; _n++) {
   point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;   
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;    
#line 41 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
sum += _const_cm*val(s,0,0,0);
 
#line 420 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
    point.i = (_i + 2)/2;
    point.j = (_j + 2)/2;
    point.k = (_k + 2)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2, point.n.z /= 2;
  }  
#line 42 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
val(s,0,0,0) = sum/(1 << 3)/(_const_cm + 1e-30);
}}}

static inline void face_average (Point point, vector v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;   
#line 47 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{







      val(v.x,0,0,0) = (fine(v.x,0,0,0) + fine(v.x,0,1,0) +
        fine(v.x,0,0,1) + fine(v.x,0,1,1))/4.;
      val(v.x,1,0,0) = (fine(v.x,2,0,0) + fine(v.x,2,1,0) +
  fine(v.x,2,0,1) + fine(v.x,2,1,1))/4.;

  } 
#line 47
{







      val(v.y,0,0,0) = (fine(v.y,0,0,0) + fine(v.y,0,0,1) +
        fine(v.y,1,0,0) + fine(v.y,1,0,1))/4.;
      val(v.y,0,1,0) = (fine(v.y,0,2,0) + fine(v.y,0,2,1) +
  fine(v.y,1,2,0) + fine(v.y,1,2,1))/4.;

  } 
#line 47
{







      val(v.z,0,0,0) = (fine(v.z,0,0,0) + fine(v.z,1,0,0) +
        fine(v.z,0,1,0) + fine(v.z,1,1,0))/4.;
      val(v.z,0,0,1) = (fine(v.z,0,0,2) + fine(v.z,1,0,2) +
  fine(v.z,0,1,2) + fine(v.z,1,1,2))/4.;

  }
}

static inline void restriction_face (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 65 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
face_average (point, _attribute[s.i].v);
}

static inline void restriction_vertex (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 70 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
for (int i = 0; i <= 1; i++) {
    val(s,i,0,0) = fine(s,2*i,0,0);

    val(s,i,1,0) = fine(s,2*i,2,0);


    for (int j = 0; j <= 1; j++)
      val(s,i,j,1) = fine(s,2*i,2*j,2);

  }
}

static inline void no_restriction (Point point, scalar s) {}

static inline void no_data (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
  
  
#line 408
{
    int _i = 2*point.i - 2;
    int _j = 2*point.j - 2;
    int _k = 2*point.k - 2;
    point.level++;
    point.n.x *= 2, point.n.y *= 2, point.n.z *= 2;
    for (int _l = 0; _l < 2; _l++)
      for (int _m = 0; _m < 2; _m++)
 for (int _n = 0; _n < 2; _n++) {
   point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;   
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;    
#line 86 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
val(s,0,0,0) = 1e30f;
 
#line 420 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
    point.i = (_i + 2)/2;
    point.j = (_j + 2)/2;
    point.k = (_k + 2)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2, point.n.z /= 2;
  }
#line 87 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
}

void wavelet (scalar s, scalar w)
{
  restriction (((scalar[]){s,{-1}}));
  for (int l = grid->maxdepth - 1; l >= 0; l--) {
  
  
    
#line 347 "/home/pwachara/basilisk/src/grid/stencils.h"
if (0) {

    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point); 
#line 93 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{
  
      
#line 408 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
    int _i = 2*point.i - 2;
    int _j = 2*point.j - 2;
    int _k = 2*point.k - 2;
    point.level++;
    point.n.x *= 2, point.n.y *= 2, point.n.z *= 2;
    for (int _l = 0; _l < 2; _l++)
      for (int _m = 0; _m < 2; _m++)
 for (int _n = 0; _n < 2; _n++) {
   point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;   
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;        
#line 95 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{ _stencil_val(s,0,0,0);_stencil_val_a(w,0,0,0); }
 
#line 420 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
    point.i = (_i + 2)/2;
    point.j = (_j + 2)/2;
    point.k = (_k + 2)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2, point.n.z /= 2;
  }
#line 96 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
default_stencil (      point,((scalar[]){ s,{-1}}));
  
      
#line 408 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
    int _i = 2*point.i - 2;
    int _j = 2*point.j - 2;
    int _k = 2*point.k - 2;
    point.level++;
    point.n.x *= 2, point.n.y *= 2, point.n.z *= 2;
    for (int _l = 0; _l < 2; _l++)
      for (int _m = 0; _m < 2; _m++)
 for (int _n = 0; _n < 2; _n++) {
   point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;   
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2; 
#line 97 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{  
         _stencil_val(s,0,0,0); 
_stencil_val(w,0,0,0);        _stencil_val_a(s,0,0,0); 

        _stencil_val_r(w,0,0,0);  
      }
 
#line 420 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
    point.i = (_i + 2)/2;
    point.j = (_j + 2)/2;
    point.k = (_k + 2)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2, point.n.z /= 2;
  }    
#line 103 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
}
  
#line 352 "/home/pwachara/basilisk/src/grid/stencils.h"
}
#line 13 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{ 
#line 244 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = l;
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
    for (_k = 2; _k < point.n.x + 2; _k++) {
      point.i = _k;

      for (point.j = 2; point.j < point.n.y + 2; point.j++)

 for (point.k = 2; point.k < point.n.z + 2; point.k++) 
#line 93 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 93 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{
  
      
#line 408 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
    int _i = 2*point.i - 2;
    int _j = 2*point.j - 2;
    int _k = 2*point.k - 2;
    point.level++;
    point.n.x *= 2, point.n.y *= 2, point.n.z *= 2;
    for (int _l = 0; _l < 2; _l++)
      for (int _m = 0; _m < 2; _m++)
 for (int _n = 0; _n < 2; _n++) {
   point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;   
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;        
#line 95 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
val(w,0,0,0) = val(s,0,0,0);
 
#line 420 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
    point.i = (_i + 2)/2;
    point.j = (_j + 2)/2;
    point.k = (_k + 2)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2, point.n.z /= 2;
  }      
#line 96 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
_attribute[s.i].prolongation (point, s);
  
      
#line 408 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
    int _i = 2*point.i - 2;
    int _j = 2*point.j - 2;
    int _k = 2*point.k - 2;
    point.level++;
    point.n.x *= 2, point.n.y *= 2, point.n.z *= 2;
    for (int _l = 0; _l < 2; _l++)
      for (int _m = 0; _m < 2; _m++)
 for (int _n = 0; _n < 2; _n++) {
   point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;   
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2; 
#line 97 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{
        double sp = val(s,0,0,0);
        val(s,0,0,0) = val(w,0,0,0);

        val(w,0,0,0) -= sp;
      }
 
#line 420 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
    point.i = (_i + 2)/2;
    point.j = (_j + 2)/2;
    point.k = (_k + 2)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2, point.n.z /= 2;
  }    
#line 103 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
}}
    
#line 261 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
  }
}
#line 16 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
}
    
#line 104
boundary_level (((scalar[]){w,{-1}}), l + 1);
  }
  

  
#line 347 "/home/pwachara/basilisk/src/grid/stencils.h"
if (0) {

    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point);    
#line 108 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{ _stencil_val(s,0,0,0);_stencil_val_a(w,0,0,0); }
  
#line 352 "/home/pwachara/basilisk/src/grid/stencils.h"
} 
#line 244 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = 0;
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
    for (_k = 2; _k < point.n.x + 2; _k++) {
      point.i = _k;

      for (point.j = 2; point.j < point.n.y + 2; point.j++)

 for (point.k = 2; point.k < point.n.z + 2; point.k++)    
#line 108 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 108 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
val(w,0,0,0) = val(s,0,0,0);}
    
#line 261 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
  }
}  
#line 109 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
boundary_level (((scalar[]){w,{-1}}), 0);
}

void inverse_wavelet (scalar s, scalar w)
{
  
  
#line 347 "/home/pwachara/basilisk/src/grid/stencils.h"
if (0) {

    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point);    
#line 115 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{ _stencil_val(w,0,0,0);_stencil_val_a(s,0,0,0); }
  
#line 352 "/home/pwachara/basilisk/src/grid/stencils.h"
} 
#line 244 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = 0;
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
    for (_k = 2; _k < point.n.x + 2; _k++) {
      point.i = _k;

      for (point.j = 2; point.j < point.n.y + 2; point.j++)

 for (point.k = 2; point.k < point.n.z + 2; point.k++)    
#line 115 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 115 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
val(s,0,0,0) = val(w,0,0,0);}
    
#line 261 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
  }
}  
#line 116 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
boundary_level (((scalar[]){s,{-1}}), 0);
  for (int l = 0; l <= grid->maxdepth - 1; l++) {
  
  
    
#line 347 "/home/pwachara/basilisk/src/grid/stencils.h"
if (0) {

    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point); 
#line 118 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{
default_stencil (      point,((scalar[]){ s,{-1}}));
  
      
#line 408 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
    int _i = 2*point.i - 2;
    int _j = 2*point.j - 2;
    int _k = 2*point.k - 2;
    point.level++;
    point.n.x *= 2, point.n.y *= 2, point.n.z *= 2;
    for (int _l = 0; _l < 2; _l++)
      for (int _m = 0; _m < 2; _m++)
 for (int _n = 0; _n < 2; _n++) {
   point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;   
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;        
#line 121 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{ _stencil_val(w,0,0,0);_stencil_val_r(s,0,0,0); }
 
#line 420 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
    point.i = (_i + 2)/2;
    point.j = (_j + 2)/2;
    point.k = (_k + 2)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2, point.n.z /= 2;
  }    
#line 122 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
}
  
#line 352 "/home/pwachara/basilisk/src/grid/stencils.h"
}
#line 13 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{ 
#line 244 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = l;
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
    for (_k = 2; _k < point.n.x + 2; _k++) {
      point.i = _k;

      for (point.j = 2; point.j < point.n.y + 2; point.j++)

 for (point.k = 2; point.k < point.n.z + 2; point.k++) 
#line 118 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 118 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{
      _attribute[s.i].prolongation (point, s);
  
      
#line 408 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
    int _i = 2*point.i - 2;
    int _j = 2*point.j - 2;
    int _k = 2*point.k - 2;
    point.level++;
    point.n.x *= 2, point.n.y *= 2, point.n.z *= 2;
    for (int _l = 0; _l < 2; _l++)
      for (int _m = 0; _m < 2; _m++)
 for (int _n = 0; _n < 2; _n++) {
   point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;   
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;        
#line 121 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
val(s,0,0,0) += val(w,0,0,0);
 
#line 420 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
    point.i = (_i + 2)/2;
    point.j = (_j + 2)/2;
    point.k = (_k + 2)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2, point.n.z /= 2;
  }    
#line 122 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
}}
    
#line 261 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
  }
}
#line 16 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
}
    
#line 123
boundary_level (((scalar[]){s,{-1}}), l + 1);
  }
}

static inline double bilinear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;    
#line 136 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
return (27.*coarse(s,0,0,0) +
     9.*(coarse(s,child.x,0,0) + coarse(s,0,child.y,0) +
  coarse(s,0,0,child.z)) +
     3.*(coarse(s,child.x,child.y,0) + coarse(s,child.x,0,child.z) +
  coarse(s,0,child.y,child.z)) +
     coarse(s,child.x,child.y,child.z))/64.;

}
#line 127
static void _stencil_bilinear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 136 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
_stencil_coarse(s,0,0,0);
_stencil_coarse(s,o_stencil,0,0); _stencil_coarse(s,0,o_stencil,0);
  _stencil_coarse(s,0,0,o_stencil);
_stencil_coarse(s,o_stencil,o_stencil,0); _stencil_coarse(s,o_stencil,0,o_stencil);
  _stencil_coarse(s,0,o_stencil,o_stencil);
     _stencil_coarse(s,o_stencil,o_stencil,o_stencil);    
#line 136
return 
        
         


;

}

static inline void refine_bilinear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
  
  
#line 408
{
    int _i = 2*point.i - 2;
    int _j = 2*point.j - 2;
    int _k = 2*point.k - 2;
    point.level++;
    point.n.x *= 2, point.n.y *= 2, point.n.z *= 2;
    for (int _l = 0; _l < 2; _l++)
      for (int _m = 0; _m < 2; _m++)
 for (int _n = 0; _n < 2; _n++) {
   point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;   
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;    
#line 148 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
val(s,0,0,0) = bilinear (point, s);
 
#line 420 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
    point.i = (_i + 2)/2;
    point.j = (_j + 2)/2;
    point.k = (_k + 2)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2, point.n.z /= 2;
  }
#line 149 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
}

static inline double quadratic (double a, double b, double c)
{
  return (30.*a + 5.*b - 3.*c)/32.;
}

static inline double biquadratic (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 172 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
  if (!(false)) qassert ("/home/pwachara/basilisk/src/grid/multigrid-common.h", 172, "false");
  return 0.;

}

static inline double biquadratic_vertex (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 185 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
if (!(false)) qassert ("/home/pwachara/basilisk/src/grid/multigrid-common.h", 185, "false");
  return 0.;

}

static inline void refine_biquadratic (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
  
  
#line 408
{
    int _i = 2*point.i - 2;
    int _j = 2*point.j - 2;
    int _k = 2*point.k - 2;
    point.level++;
    point.n.x *= 2, point.n.y *= 2, point.n.z *= 2;
    for (int _l = 0; _l < 2; _l++)
      for (int _m = 0; _m < 2; _m++)
 for (int _n = 0; _n < 2; _n++) {
   point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;   
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;    
#line 193 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
val(s,0,0,0) = biquadratic (point, s);
 
#line 420 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
    point.i = (_i + 2)/2;
    point.j = (_j + 2)/2;
    point.k = (_k + 2)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2, point.n.z /= 2;
  }
#line 194 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
}

static inline void refine_linear (Point point, scalar s)
{if(!is_constant(cm)){{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 198 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
coord g;
  if (_attribute[s.i].gradient)
    {
      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));      g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));      g.z = _attribute[s.i].gradient (val(s,0,0,-1), val(s,0,0,0), val(s,0,0,1));}
  else
    {
      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;      g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;      g.z = (val(s,0,0,1) - val(s,0,0,-1))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val(cm,0,0,0), sum = val(cm,0,0,0)*(1 << 3);
  
  
#line 408 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
    int _i = 2*point.i - 2;
    int _j = 2*point.j - 2;
    int _k = 2*point.k - 2;
    point.level++;
    point.n.x *= 2, point.n.y *= 2, point.n.z *= 2;
    for (int _l = 0; _l < 2; _l++)
      for (int _m = 0; _m < 2; _m++)
 for (int _n = 0; _n < 2; _n++) {
   point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;   
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2; 
#line 207 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{
    val(s,0,0,0) = sc;
    
      val(s,0,0,0) += child.x*g.x*val(cm,-child.x,0,0)/cmc;      val(s,0,0,0) += child.y*g.y*val(cm,0,-child.y,0)/cmc;      val(s,0,0,0) += child.z*g.z*val(cm,0,0,-child.z)/cmc;
    sum -= val(cm,0,0,0);
  }
 
#line 420 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
    point.i = (_i + 2)/2;
    point.j = (_j + 2)/2;
    point.k = (_k + 2)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2, point.n.z /= 2;
  }  
#line 213 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
if (!(fabs(sum) < 1e-10)) qassert ("/home/pwachara/basilisk/src/grid/multigrid-common.h", 213, "fabs(sum) < 1e-10");
}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
#line 197
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 198 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
coord g;
  if (_attribute[s.i].gradient)
    {
      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));      g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));      g.z = _attribute[s.i].gradient (val(s,0,0,-1), val(s,0,0,0), val(s,0,0,1));}
  else
    {
      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;      g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;      g.z = (val(s,0,0,1) - val(s,0,0,-1))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*_const_cm, sum = _const_cm*(1 << 3);
  
  
#line 408 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
    int _i = 2*point.i - 2;
    int _j = 2*point.j - 2;
    int _k = 2*point.k - 2;
    point.level++;
    point.n.x *= 2, point.n.y *= 2, point.n.z *= 2;
    for (int _l = 0; _l < 2; _l++)
      for (int _m = 0; _m < 2; _m++)
 for (int _n = 0; _n < 2; _n++) {
   point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;   
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2; 
#line 207 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{
    val(s,0,0,0) = sc;
    
      val(s,0,0,0) += child.x*g.x*_const_cm/cmc;      val(s,0,0,0) += child.y*g.y*_const_cm/cmc;      val(s,0,0,0) += child.z*g.z*_const_cm/cmc;
    sum -= _const_cm;
  }
 
#line 420 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
    point.i = (_i + 2)/2;
    point.j = (_j + 2)/2;
    point.k = (_k + 2)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2, point.n.z /= 2;
  }  
#line 213 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
if (!(fabs(sum) < 1e-10)) qassert ("/home/pwachara/basilisk/src/grid/multigrid-common.h", 213, "fabs(sum) < 1e-10");
}}}

static inline void refine_reset (Point point, scalar v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
  
  
#line 408
{
    int _i = 2*point.i - 2;
    int _j = 2*point.j - 2;
    int _k = 2*point.k - 2;
    point.level++;
    point.n.x *= 2, point.n.y *= 2, point.n.z *= 2;
    for (int _l = 0; _l < 2; _l++)
      for (int _m = 0; _m < 2; _m++)
 for (int _n = 0; _n < 2; _n++) {
   point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;   
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
    
#line 219 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
val(v,0,0,0) = 0.;
 
#line 420 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
    point.i = (_i + 2)/2;
    point.j = (_j + 2)/2;
    point.k = (_k + 2)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2, point.n.z /= 2;
  }
#line 220 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
}

static inline void refine_injection (Point point, scalar v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
  
#line 224 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
double val = val(v,0,0,0);
  
  
#line 408 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
    int _i = 2*point.i - 2;
    int _j = 2*point.j - 2;
    int _k = 2*point.k - 2;
    point.level++;
    point.n.x *= 2, point.n.y *= 2, point.n.z *= 2;
    for (int _l = 0; _l < 2; _l++)
      for (int _m = 0; _m < 2; _m++)
 for (int _n = 0; _n < 2; _n++) {
   point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;   
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
    
#line 226 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
val(v,0,0,0) = val;
 
#line 420 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
    point.i = (_i + 2)/2;
    point.j = (_j + 2)/2;
    point.k = (_k + 2)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2, point.n.z /= 2;
  }
#line 227 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
}

static scalar multigrid_init_scalar (scalar s, const char * name)
{
  s = cartesian_init_scalar (s, name);
  _attribute[s.i].prolongation = refine_bilinear;
  _attribute[s.i].restriction = restriction_average;
  return s;
}

static scalar multigrid_init_vertex_scalar (scalar s, const char * name)
{
  s = cartesian_init_vertex_scalar (s, name);
  _attribute[s.i].restriction = restriction_vertex;
  return s;
}

static void multigrid_setup_vector (vector v)
{
   {
    _attribute[v.x.i].prolongation = refine_bilinear;
    _attribute[v.x.i].restriction = restriction_average;
  } 
#line 246
{
    _attribute[v.y.i].prolongation = refine_bilinear;
    _attribute[v.y.i].restriction = restriction_average;
  } 
#line 246
{
    _attribute[v.z.i].prolongation = refine_bilinear;
    _attribute[v.z.i].restriction = restriction_average;
  }
}

static vector multigrid_init_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
  multigrid_setup_vector (v);
  return v;
}

static vector multigrid_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  
    _attribute[v.y.i].restriction = no_restriction;    _attribute[v.z.i].restriction = no_restriction;    _attribute[v.x.i].restriction = no_restriction;
  _attribute[v.x.i].restriction = restriction_face;
  return v;
}

static tensor multigrid_init_tensor (tensor t, const char * name)
{
  t = cartesian_init_tensor (t, name);
  
    multigrid_setup_vector (t.x);    multigrid_setup_vector (t.y);    multigrid_setup_vector (t.z);
  return t;
}

void multigrid_debug (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
  
#line 278 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
cartesian_debug (point);

  FILE * plot = fopen ("plot", "a");
  if (point.level > 0) {
    char name[80] = "coarse";
    if (pid() > 0)
      sprintf (name, "coarse-%d", pid());
    FILE * fp = fopen (name, "w");
#line 310 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
      double xc = x - child.x*Delta/2., yc = y - child.y*Delta/2.;
      double zc = z - child.z*Delta/2.;
      for (int k = 0; k <= 1; k++)
 for (int l = 0; l <= 1; l++)
   for (int m = 0; m <= 1; m++) {
     {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){
       fprintf (fp, "%g %g %g %g ",
         xc + k*child.x*Delta*2. + _attribute[v.i].d.x*Delta,
         yc + l*child.y*Delta*2. + _attribute[v.i].d.y*Delta,
         zc + m*child.z*Delta*2. + _attribute[v.i].d.z*Delta,
         coarse(v,k*child.x,l*child.y,m*child.z));}}
     fputc ('\n', fp);
   }
      fprintf (ferr, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 3 t ''",
        name);
      fprintf (plot, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 3 t ''",
        name);

    fclose (fp);
  }

  if (is_coarse()) {
    char name[80] = "fine";
    if (pid() > 0)
      sprintf (name, "fine-%d", pid());
    FILE * fp = fopen (name, "w");
#line 367 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
      double xf = x - Delta/4., yf = y - Delta/4., zf = z - Delta/4.;
      for (int k = -2; k <= 3; k++)
 for (int l = -2; l <= 3; l++)
   for (int m = -2; m <= 3; m++) {
     {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){ {
       fprintf (fp, "%g %g %g ",
         xf + k*Delta/2. + _attribute[v.i].d.x*Delta/4.,
         yf + l*Delta/2. + _attribute[v.i].d.y*Delta/4.,
         zf + m*Delta/2. + _attribute[v.i].d.z*Delta/4.);
       if (allocated_child(k,l,m))
  fprintf (fp, "%g ", fine(v,k,l,m));
       else
  fputs ("n/a ", fp);
     }}}
     fputc ('\n', fp);
   }
      fprintf (ferr, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 2 t ''",
        name);
      fprintf (plot, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 2 t ''",
        name);

    fclose (fp);
  }
  fflush (ferr);
  fclose (plot);
}

static void multigrid_restriction (scalar * list)
{
  scalar * listdef = NULL, * listc = NULL, * list2 = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant (s) && _attribute[s.i].block > 0) {
      if (_attribute[s.i].restriction == restriction_average) {
 listdef = list_add (listdef, s);
 list2 = list_add (list2, s);
      }
      else if (_attribute[s.i].restriction != no_restriction) {
 listc = list_add (listc, s);
 if (_attribute[s.i].face)
   {
     list2 = list_add (list2, _attribute[s.i].v.x);     list2 = list_add (list2, _attribute[s.i].v.y);     list2 = list_add (list2, _attribute[s.i].v.z);}
 else
   list2 = list_add (list2, s);
      }
    }}}

  if (listdef || listc) {
    for (int l = depth() - 1; l >= 0; l--) {      
#line 347 "/home/pwachara/basilisk/src/grid/stencils.h"
if (0) {

    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point); 
#line 415 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{
 {scalar*_i=(scalar*)( listdef);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     _stencil_restriction_average (point, s);}}
 {scalar*_i=(scalar*)( listc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {

default_stencil (     point,((scalar[]){ s,{-1}}));
 }}}
      }  
#line 352 "/home/pwachara/basilisk/src/grid/stencils.h"
}
#line 13 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{ 
#line 244 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = l;
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
    for (_k = 2; _k < point.n.x + 2; _k++) {
      point.i = _k;

      for (point.j = 2; point.j < point.n.y + 2; point.j++)

 for (point.k = 2; point.k < point.n.z + 2; point.k++) 
#line 415 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 415 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{
 {scalar*_i=(scalar*)( listdef);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     restriction_average (point, s);}}
 {scalar*_i=(scalar*)( listc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
  
     _attribute[s.i].restriction (point, s);
 }}}
      }}    
#line 261 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
  }
}
#line 16 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
}
      
#line 424
{ Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list2, l); };
    }
    pfree (listdef,__func__,__FILE__,__LINE__);
    pfree (listc,__func__,__FILE__,__LINE__);
    pfree (list2,__func__,__FILE__,__LINE__);
  }
}

void multigrid_methods()
{
  cartesian_methods();
  init_scalar = multigrid_init_scalar;
  init_vertex_scalar = multigrid_init_vertex_scalar;
  init_vector = multigrid_init_vector;
  init_face_vector = multigrid_init_face_vector;
  init_tensor = multigrid_init_tensor;
  restriction = multigrid_restriction;
  debug = multigrid_debug;
}







void subtree_size (scalar size, bool leaves)
{  
#line 314 "/home/pwachara/basilisk/src/grid/stencils.h"
{
    static int _first = 1.;
    ForeachData _loop = {
      .fname = "/home/pwachara/basilisk/src/grid/multigrid-common.h", .line = 456, .first = _first
    };
    if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
 _attribute[s.i].input = _attribute[s.i].output = _attribute[s.i].nowarning = false;
 _attribute[s.i].width = 0;
      }
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point);
    
#line 457 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{_stencil_val_a(size,0,0,0);  }    
#line 328 "/home/pwachara/basilisk/src/grid/stencils.h"
check_stencil (&_loop);
    boundary_stencil (&_loop);
    _first = 0;
  } 
#line 265 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
      for (_k = 2; _k < point.n.x + 2; _k++) {
 point.i = _k;

 for (point.j = 2; point.j < point.n.y + 2; point.j++)

   for (point.k = 2; point.k < point.n.z + 2; point.k++)
    
#line 457 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 457 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
val(size,0,0,0) = 1;}      
#line 282 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
  }
}





  
#line 463 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{ Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b,((scalar[]) {size,{-1}}), depth()); };
  for (int l = depth() - 1; l >= 0; l--) {    
#line 347 "/home/pwachara/basilisk/src/grid/stencils.h"
if (0) {

    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point); 
#line 465 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{      
#line 408 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
    int _i = 2*point.i - 2;
    int _j = 2*point.j - 2;
    int _k = 2*point.k - 2;
    point.level++;
    point.n.x *= 2, point.n.y *= 2, point.n.z *= 2;
    for (int _l = 0; _l < 2; _l++)
      for (int _m = 0; _m < 2; _m++)
 for (int _n = 0; _n < 2; _n++) {
   point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;   
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
 
#line 468 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{ _stencil_val(size,0,0,0); } 
#line 420 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
    point.i = (_i + 2)/2;
    point.j = (_j + 2)/2;
    point.k = (_k + 2)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2, point.n.z /= 2;
  }
      
#line 469 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
_stencil_val_a(size,0,0,0);  
    }  
#line 352 "/home/pwachara/basilisk/src/grid/stencils.h"
}
#line 13 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{ 
#line 244 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = l;
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
    for (_k = 2; _k < point.n.x + 2; _k++) {
      point.i = _k;

      for (point.j = 2; point.j < point.n.y + 2; point.j++)

 for (point.k = 2; point.k < point.n.z + 2; point.k++) 
#line 465 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 465 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
{
      double sum = !leaves;      
#line 408 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
    int _i = 2*point.i - 2;
    int _j = 2*point.j - 2;
    int _k = 2*point.k - 2;
    point.level++;
    point.n.x *= 2, point.n.y *= 2, point.n.z *= 2;
    for (int _l = 0; _l < 2; _l++)
      for (int _m = 0; _m < 2; _m++)
 for (int _n = 0; _n < 2; _n++) {
   point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;   
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
 
#line 468 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
sum += val(size,0,0,0); 
#line 420 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
    point.i = (_i + 2)/2;
    point.j = (_j + 2)/2;
    point.k = (_k + 2)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2, point.n.z /= 2;
  }
      
#line 469 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
val(size,0,0,0) = sum;
    }}    
#line 261 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
  }
}
#line 16 "/home/pwachara/basilisk/src/grid/multigrid-common.h"
}
    
#line 471
{ Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b,((scalar[]) {size,{-1}}), l); };
  }
}
#line 910 "/home/pwachara/basilisk/src/grid/multigrid.h"


void macro_foreach_vertex (char flags, Reduce reductions) {
#line 300
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
      for (_k = 2; _k <= point.n.x + 2; _k++) {
 point.i = _k;

 for (point.j = 2; point.j <= point.n.y + 2; point.j++)

   for (point.k = 2; point.k <= point.n.z + 2; point.k++) 
#line 913
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 913
{
    int ig = -1; NOT_UNUSED(ig);

    int jg = -1; NOT_UNUSED(jg);


    int kg = -1; NOT_UNUSED(kg);

    ;
  }}      
#line 317
}
  }
}

#line 923
}


void macro_foreach_edge (char flags, Reduce reductions) { 
#line 912
{
#line 300
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
      for (_k = 2; _k <= point.n.x + 2; _k++) {
 point.i = _k;

 for (point.j = 2; point.j <= point.n.y + 2; point.j++)

   for (point.k = 2; point.k <= point.n.z + 2; point.k++) 
#line 913
{
    int ig = -1; NOT_UNUSED(ig);

    int jg = -1; NOT_UNUSED(jg);


    int kg = -1; NOT_UNUSED(kg); 







{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 927
{
    struct { int x, y, z; } _a = {point.i, point.j, point.k};
    
      if (_a.x < point.n.x + 2)
 ;      
#line 930
if (_a.y < point.n.y + 2)
 ;      
#line 930
if (_a.z < point.n.z + 2)
 ;
  }}  
#line 922
}      
#line 317
}
  }
}

#line 923
}









}


ivec dimensions (int nx, int ny, int nz)
{
  if (nx != 0 || ny != 0 || nz != 0) {
    Dimensions.x = Dimensions_scale = ( nx > 1 ? nx : 1);

    Dimensions.y = ( ny > 1 ? ny : 1);


    Dimensions.z = ( nz > 1 ? nz : 1);

  }
  return Dimensions;
}


#line 1 "grid/multigrid-mpi.h"
#line 1 "/home/pwachara/basilisk/src/grid/multigrid-mpi.h"
#line 41 "/home/pwachara/basilisk/src/grid/multigrid-mpi.h"
void macro_foreach_slice_x (int start, int end, int l) {
  {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = l; point.n.x = point.n.y = point.n.z = 1 << point.level;
    for (point.i = start; point.i < end; point.i++)
      for (point.j = 0; point.j < point.n.y + 2*2; point.j++)
 for (point.k = 0; point.k < point.n.z + 2*2; point.k++)
   ;
  }
}

void macro_foreach_slice_y (int start, int end, int l) {
  {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = l; point.n.x = point.n.y = point.n.z = 1 << point.level;
    for (point.i = 0; point.i < point.n.x + 2*2; point.i++)
      for (point.j = start; point.j < end; point.j++)
 for (point.k = 0; point.k < point.n.z + 2*2; point.k++)
   ;
  }
}

void macro_foreach_slice_z (int start, int end, int l) {
  {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = l; point.n.x = point.n.y = point.n.z = 1 << point.level;
    for (point.i = 0; point.i < point.n.x + 2*2; point.i++)
      for (point.j = 0; point.j < point.n.y + 2*2; point.j++)
 for (point.k = start; point.k < end; point.k++)
   ;
  }
}



typedef struct {
  Boundary b;
  MPI_Comm cartcomm;
} MpiBoundary;


static void * snd_x (int i, int dst, int tag, int level, scalar * list,
       MPI_Request * req)
{
  if (dst == MPI_PROC_NULL)
    return NULL;
  size_t size = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    size += _attribute[s.i].block;}}
  size *= pow((1 << level) + 2*2, 3 - 1)*2*sizeof(real);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf; 
#line 41
{
  {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = level; point.n.x = point.n.y = point.n.z = 1 << point.level;
    for (point.i = i; point.i < (
#line 95
i + 2
#line 46
); point.i++)
      for (point.j = 0; point.j < point.n.y + 2*2; point.j++)
 for (point.k = 0; point.k < point.n.z + 2*2; point.k++)
    
#line 96
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 96 "/home/pwachara/basilisk/src/grid/multigrid-mpi.h"
{scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      for (scalar sb = s; sb.i < s.i + _attribute[s.i].block; sb.i++, b++)
 memcpy (b, &val(sb,0,0,0), sizeof(real));}}}  
#line 50
}
}
  
#line 99
MPI_Isend (buf, size, MPI_BYTE, dst, tag, MPI_COMM_WORLD, req);
  return buf;
}
#line 85
static void * snd_y (int i, int dst, int tag, int level, scalar * list,
       MPI_Request * req)
{
  if (dst == MPI_PROC_NULL)
    return NULL;
  size_t size = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    size += _attribute[s.i].block;}}
  size *= pow((1 << level) + 2*2, 3 - 1)*2*sizeof(real);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf; 
#line 53
{
  {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = level; point.n.x = point.n.y = point.n.z = 1 << point.level;
    for (point.i = 0; point.i < point.n.x + 2*2; point.i++)
      for (point.j = i; point.j < (
#line 95
i + 2
#line 59
); point.j++)
 for (point.k = 0; point.k < point.n.z + 2*2; point.k++)
    
#line 96
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_y = Delta;

  double Delta_z = Delta;


  double Delta_x = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_y);

  NOT_UNUSED(Delta_z);


  NOT_UNUSED(Delta_x);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 96 "/home/pwachara/basilisk/src/grid/multigrid-mpi.h"
{scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      for (scalar sb = s; sb.i < s.i + _attribute[s.i].block; sb.i++, b++)
 memcpy (b, &val(sb,0,0,0), sizeof(real));}}}  
#line 62
}
}
  
#line 99
MPI_Isend (buf, size, MPI_BYTE, dst, tag, MPI_COMM_WORLD, req);
  return buf;
}
#line 85
static void * snd_z (int i, int dst, int tag, int level, scalar * list,
       MPI_Request * req)
{
  if (dst == MPI_PROC_NULL)
    return NULL;
  size_t size = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    size += _attribute[s.i].block;}}
  size *= pow((1 << level) + 2*2, 3 - 1)*2*sizeof(real);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf; 
#line 65
{
  {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = level; point.n.x = point.n.y = point.n.z = 1 << point.level;
    for (point.i = 0; point.i < point.n.x + 2*2; point.i++)
      for (point.j = 0; point.j < point.n.y + 2*2; point.j++)
 for (point.k = i; point.k < (
#line 95
i + 2
#line 72
); point.k++)
    
#line 96
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_z = Delta;

  double Delta_x = Delta;


  double Delta_y = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_z);

  NOT_UNUSED(Delta_x);


  NOT_UNUSED(Delta_y);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 96 "/home/pwachara/basilisk/src/grid/multigrid-mpi.h"
{scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      for (scalar sb = s; sb.i < s.i + _attribute[s.i].block; sb.i++, b++)
 memcpy (b, &val(sb,0,0,0), sizeof(real));}}}  
#line 74
}
}
  
#line 99
MPI_Isend (buf, size, MPI_BYTE, dst, tag, MPI_COMM_WORLD, req);
  return buf;
}


static void rcv_x (int i, int src, int tag, int level, scalar * list)
{
  if (src == MPI_PROC_NULL)
    return;
  size_t size = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    size += _attribute[s.i].block;}}
  size *= pow((1 << level) + 2*2, 3 - 1)*2*sizeof(real);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf;
  MPI_Status s;
  MPI_Recv (buf, size, MPI_BYTE, src, tag, MPI_COMM_WORLD, &s); 
#line 41
{
  {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = level; point.n.x = point.n.y = point.n.z = 1 << point.level;
    for (point.i = i; point.i < (
#line 115
i + 2
#line 46
); point.i++)
      for (point.j = 0; point.j < point.n.y + 2*2; point.j++)
 for (point.k = 0; point.k < point.n.z + 2*2; point.k++)
    
#line 116
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 116 "/home/pwachara/basilisk/src/grid/multigrid-mpi.h"
{scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      for (scalar sb = s; sb.i < s.i + _attribute[s.i].block; sb.i++, b++)
 memcpy (&val(sb,0,0,0), b, sizeof(real));}}}  
#line 50
}
}
  
#line 119
pfree (buf,__func__,__FILE__,__LINE__);
}
#line 104
static void rcv_y (int i, int src, int tag, int level, scalar * list)
{
  if (src == MPI_PROC_NULL)
    return;
  size_t size = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    size += _attribute[s.i].block;}}
  size *= pow((1 << level) + 2*2, 3 - 1)*2*sizeof(real);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf;
  MPI_Status s;
  MPI_Recv (buf, size, MPI_BYTE, src, tag, MPI_COMM_WORLD, &s); 
#line 53
{
  {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = level; point.n.x = point.n.y = point.n.z = 1 << point.level;
    for (point.i = 0; point.i < point.n.x + 2*2; point.i++)
      for (point.j = i; point.j < (
#line 115
i + 2
#line 59
); point.j++)
 for (point.k = 0; point.k < point.n.z + 2*2; point.k++)
    
#line 116
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_y = Delta;

  double Delta_z = Delta;


  double Delta_x = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_y);

  NOT_UNUSED(Delta_z);


  NOT_UNUSED(Delta_x);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 116 "/home/pwachara/basilisk/src/grid/multigrid-mpi.h"
{scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      for (scalar sb = s; sb.i < s.i + _attribute[s.i].block; sb.i++, b++)
 memcpy (&val(sb,0,0,0), b, sizeof(real));}}}  
#line 62
}
}
  
#line 119
pfree (buf,__func__,__FILE__,__LINE__);
}
#line 104
static void rcv_z (int i, int src, int tag, int level, scalar * list)
{
  if (src == MPI_PROC_NULL)
    return;
  size_t size = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    size += _attribute[s.i].block;}}
  size *= pow((1 << level) + 2*2, 3 - 1)*2*sizeof(real);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf;
  MPI_Status s;
  MPI_Recv (buf, size, MPI_BYTE, src, tag, MPI_COMM_WORLD, &s); 
#line 65
{
  {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = level; point.n.x = point.n.y = point.n.z = 1 << point.level;
    for (point.i = 0; point.i < point.n.x + 2*2; point.i++)
      for (point.j = 0; point.j < point.n.y + 2*2; point.j++)
 for (point.k = i; point.k < (
#line 115
i + 2
#line 72
); point.k++)
    
#line 116
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_z = Delta;

  double Delta_x = Delta;


  double Delta_y = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_z);

  NOT_UNUSED(Delta_x);


  NOT_UNUSED(Delta_y);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 116 "/home/pwachara/basilisk/src/grid/multigrid-mpi.h"
{scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      for (scalar sb = s; sb.i < s.i + _attribute[s.i].block; sb.i++, b++)
 memcpy (&val(sb,0,0,0), b, sizeof(real));}}}  
#line 74
}
}
  
#line 119
pfree (buf,__func__,__FILE__,__LINE__);
}

     
static void mpi_boundary_level (const Boundary * b, scalar * list, int level)
{tracing("mpi_boundary_level","/home/pwachara/basilisk/src/grid/multigrid-mpi.h",123);
  scalar * list1 = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].block > 0)
      list1 = list_add (list1, s);}}
  if (!list1)
    {end_tracing("mpi_boundary_level","/home/pwachara/basilisk/src/grid/multigrid-mpi.h",130);return;}

  prof_start ("mpi_boundary_level");

  if (level < 0) level = depth();
  MpiBoundary * mpi = (MpiBoundary *) b;
  struct { int x, y, z; } dir = {0,1,2};
   {
    int left, right;
    MPI_Cart_shift (mpi->cartcomm, dir.x, 1, &left, &right);
    MPI_Request reqs[2];
    void * buf[2];
    int npl = (1 << level) + 2*2, nr = 0;
    if ((buf[0] = snd_x (npl - 2*2, right, 0, level, list1, &reqs[nr])))
      nr++;
    if ((buf[1] = snd_x (2, left, 1, level, list1, &reqs[nr])))
      nr++;
    rcv_x (0, left, 0, level, list1);
    rcv_x (npl - 2, right, 1, level, list1);
    MPI_Status stats[nr];
    MPI_Waitall (nr, reqs, stats);
    pfree (buf[0],__func__,__FILE__,__LINE__); pfree (buf[1],__func__,__FILE__,__LINE__);
  } 
#line 137
{
    int bottom, top;
    MPI_Cart_shift (mpi->cartcomm, dir.y, 1, &bottom, &top);
    MPI_Request reqs[2];
    void * buf[2];
    int npl = (1 << level) + 2*2, nr = 0;
    if ((buf[0] = snd_y (npl - 2*2, top, 0, level, list1, &reqs[nr])))
      nr++;
    if ((buf[1] = snd_y (2, bottom, 1, level, list1, &reqs[nr])))
      nr++;
    rcv_y (0, bottom, 0, level, list1);
    rcv_y (npl - 2, top, 1, level, list1);
    MPI_Status stats[nr];
    MPI_Waitall (nr, reqs, stats);
    pfree (buf[0],__func__,__FILE__,__LINE__); pfree (buf[1],__func__,__FILE__,__LINE__);
  } 
#line 137
{
    int back, front;
    MPI_Cart_shift (mpi->cartcomm, dir.z, 1, &back, &front);
    MPI_Request reqs[2];
    void * buf[2];
    int npl = (1 << level) + 2*2, nr = 0;
    if ((buf[0] = snd_z (npl - 2*2, front, 0, level, list1, &reqs[nr])))
      nr++;
    if ((buf[1] = snd_z (2, back, 1, level, list1, &reqs[nr])))
      nr++;
    rcv_z (0, back, 0, level, list1);
    rcv_z (npl - 2, front, 1, level, list1);
    MPI_Status stats[nr];
    MPI_Waitall (nr, reqs, stats);
    pfree (buf[0],__func__,__FILE__,__LINE__); pfree (buf[1],__func__,__FILE__,__LINE__);
  }

  pfree (list1,__func__,__FILE__,__LINE__);

  prof_stop();
end_tracing("mpi_boundary_level","/home/pwachara/basilisk/src/grid/multigrid-mpi.h",157);}

static void mpi_boundary_destroy (Boundary * b)
{
  MpiBoundary * m = (MpiBoundary *) b;
  MPI_Comm_free (&m->cartcomm);
  pfree (m,__func__,__FILE__,__LINE__);
}

static void mpi_dimensions_error (int n)
{
  fprintf (ferr,
    "%s:%d: error: the number of MPI processes must be equal to ",
    "/home/pwachara/basilisk/src/grid/multigrid-mpi.h", 170);
  if (n > 1)
    fprintf (ferr, "%dx", n);
  fprintf (ferr, "%d^i\n", 1 << 3);
  exit (1);
}

Boundary * mpi_boundary_new()
{
  MpiBoundary * m = ((MpiBoundary *) pcalloc (1, sizeof(MpiBoundary),__func__,__FILE__,__LINE__));
  int n = 1;
  
    n *= Dimensions.x;    n *= Dimensions.y;    n *= Dimensions.z;
  if (npe() % n)
    mpi_dimensions_error (n);
  int j = npe()/n, i = 0;
  while (j > 1) {
    if (j % (1 << 3))
      mpi_dimensions_error (n);
    j /= 1 << 3;
    i++;
  }
  
    Dimensions.x *= 1 << i;    Dimensions.y *= 1 << i;    Dimensions.z *= 1 << i;
  MPI_Dims_create (npe(), 3, &Dimensions.x);
  MPI_Cart_create (MPI_COMM_WORLD, 3,
     &Dimensions.x, &Period.x, 0, &m->cartcomm);
  MPI_Cart_coords (m->cartcomm, pid(), 3, mpi_coords);


  struct { int x, y, z; } dir = {0,1,2};
   {
    int l, r;
    MPI_Cart_shift (m->cartcomm, dir.x, 1, &l, &r);
    if (l != MPI_PROC_NULL)
      periodic_boundary (left);
    if (r != MPI_PROC_NULL)
      periodic_boundary (right);
  } 
#line 201
{
    int l, r;
    MPI_Cart_shift (m->cartcomm, dir.y, 1, &l, &r);
    if (l != MPI_PROC_NULL)
      periodic_boundary (bottom);
    if (r != MPI_PROC_NULL)
      periodic_boundary (top);
  } 
#line 201
{
    int l, r;
    MPI_Cart_shift (m->cartcomm, dir.z, 1, &l, &r);
    if (l != MPI_PROC_NULL)
      periodic_boundary (back);
    if (r != MPI_PROC_NULL)
      periodic_boundary (front);
  }


  Dimensions_scale = Dimensions.x;
  N /= Dimensions.x;
  int r = 0;
  while (N > 1)
    N /= 2, r++;
  grid->depth = grid->maxdepth = r;
  N = Dimensions.x*(1 << r);
  grid->n = 1 << 3*depth();
  grid->tn = npe()*grid->n;


  Boundary * b = (Boundary *) m;
  b->level = mpi_boundary_level;
  b->destroy = mpi_boundary_destroy;
  add_boundary (b);

  return b;
}

     
double z_indexing (scalar index, bool leaves)
{tracing("z_indexing","/home/pwachara/basilisk/src/grid/multigrid-mpi.h",231);
  long i;
  if (leaves)
    i = pid()*(1 << 3*depth());
  else
    i = pid()*((1 << 3*(depth() + 1)) - 1)/((1 << 3) - 1);
#line 126 "/home/pwachara/basilisk/src/grid/foreach_cell.h"
{
  {





    Point root = {2,2,2,0};
#line 67
{
  {
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point = {0};





    int kg = 0; NOT_UNUSED(kg);
    struct { int l, i, j, k, stage; } stack[20];

    int _s = -1;
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].k = root.k; stack[_s].stage = 0; };
    while (_s >= 0) {
      int stage;
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; stage = stack[_s].stage; _s--; };
      if (!allocated (0,0,0))
 continue;
      switch (stage) {
      case 0: { 
#line 238 "/home/pwachara/basilisk/src/grid/multigrid-mpi.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 238 "/home/pwachara/basilisk/src/grid/multigrid-mpi.h"
{
    if (!leaves || is_leaf(cell))
      val(index,0,0,0) = i++;
    if (is_leaf(cell))
      continue;
  }} 
#line 91 "/home/pwachara/basilisk/src/grid/foreach_cell.h"
if (point.level < grid->depth) {
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 1; };
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };
 }
 break;
      }
#line 106 "/home/pwachara/basilisk/src/grid/foreach_cell.h"
      case 1: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 2; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;
      case 2: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 3; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;
      case 3: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 4; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;
      case 4: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 5; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;
      case 5: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 6; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;
      case 6: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 7; };
 { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;
      case 7: { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;

      }
    }
  }
}
  
#line 137
}
}
  
#line 244 "/home/pwachara/basilisk/src/grid/multigrid-mpi.h"
boundary_internal ((scalar *)((scalar[]){index,{-1}}), "/home/pwachara/basilisk/src/grid/multigrid-mpi.h", 244);
  { double _ret= pid() == 0 ? i*npe() - 1 : -1;end_tracing("z_indexing","/home/pwachara/basilisk/src/grid/multigrid-mpi.h",245);return _ret;}
end_tracing("z_indexing","/home/pwachara/basilisk/src/grid/multigrid-mpi.h",246);}
#line 952 "/home/pwachara/basilisk/src/grid/multigrid.h"
#line 3 "/home/pwachara/basilisk/src/grid/multigrid3D.h"

void multigrid3D_methods() {
  multigrid_methods();
}
#line 2 "BaseFlow.c"
#line 1 "embed.h"
#line 1 "/home/pwachara/basilisk/src/embed.h"
#line 12 "/home/pwachara/basilisk/src/embed.h"
#line 1 "fractions.h"
#line 1 "/home/pwachara/basilisk/src/fractions.h"
#line 12 "/home/pwachara/basilisk/src/fractions.h"
#line 1 "geometry.h"
#line 1 "/home/pwachara/basilisk/src/geometry.h"
#line 35 "/home/pwachara/basilisk/src/geometry.h"
double line_alpha (double c, coord n)
{
  double alpha, n1, n2;

  n1 = fabs (n.x); n2 = fabs (n.y);
  if (n1 > n2)
    do { double _tmp_ = n1; n1 = n2; n2 = _tmp_; } while(false);

  c = ( c < 0. ? 0. : c > 1. ? 1. : c);
  double v1 = n1/2.;
  if (c <= v1/n2)
    alpha = sqrt (2.*c*n1*n2);
  else if (c <= 1. - v1/n2)
    alpha = c*n2 + v1;
  else
    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;

  return alpha - (n.x + n.y)/2.;
}



double plane_alpha (double c, coord n)
{
  double alpha;
  coord n1;

  n1.x = fabs (n.x); n1.y = fabs (n.y); n1.z = fabs (n.z);

  double m1, m2, m3;
  m1 = ( (n1.x) < (n1.y) ? (n1.x) : (n1.y));
  m3 = ( (n1.x) > (n1.y) ? (n1.x) : (n1.y));
  m2 = n1.z;
  if (m2 < m1) {
    double tmp = m1;
    m1 = m2;
    m2 = tmp;
  }
  else if (m2 > m3) {
    double tmp = m3;
    m3 = m2;
    m2 = tmp;
  }
  double m12 = m1 + m2;
  double pr = ( (6.*m1*m2*m3) > 1e-50 ? (6.*m1*m2*m3) : 1e-50);
  double V1 = m1*m1*m1/pr;
  double V2 = V1 + (m2 - m1)/(2.*m3), V3;
  double mm;
  if (m3 < m12) {
    mm = m3;
    V3 = (m3*m3*(3.*m12 - m3) + m1*m1*(m1 - 3.*m3) + m2*m2*(m2 - 3.*m3))/pr;
  }
  else {
    mm = m12;
    V3 = mm/(2.*m3);
  }

  c = ( c < 0. ? 0. : c > 1. ? 1. : c);
  double ch = ( c < (1. - c) ? c : (1. - c));
  if (ch < V1)
    alpha = pow (pr*ch, 1./3.);
  else if (ch < V2)
    alpha = (m1 + sqrt(m1*m1 + 8.*m2*m3*(ch - V1)))/2.;
  else if (ch < V3) {
    double p12 = sqrt (2.*m1*m2);
    double q = 3.*(m12 - 2.*m3*ch)/(4.*p12);
    double teta = acos(( q < (-1.) ? (-1.) : q > 1. ? 1. : q))/3.;
    double cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + m12;
  }
  else if (m12 <= m3)
    alpha = m3*ch + mm/2.;
  else {
    double p = m1*(m2 + m3) + m2*m3 - 1./4., p12 = sqrt(p);
    double q = 3.*m1*m2*m3*(1./2. - ch)/(2.*p*p12);
    double teta = acos(( q < (-1.) ? (-1.) : q > 1. ? 1. : q))/3.;
    double cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + 1./2.;
  }
  if (c > 1./2.) alpha = 1. - alpha;

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;
  if (n.z < 0.)
    alpha += n.z;

  return alpha - (n.x + n.y + n.z)/2.;;
}
#line 163 "/home/pwachara/basilisk/src/geometry.h"
double line_area (double nx, double ny, double alpha)
{
  double a, v, area;

  alpha += (nx + ny)/2.;
  if (nx < 0.) {
    alpha -= nx;
    nx = - nx;
  }
  if (ny < 0.) {
    alpha -= ny;
    ny = - ny;
  }

  if (alpha <= 0.)
    return 0.;

  if (alpha >= nx + ny)
    return 1.;

  if (nx < 1e-10)
    area = alpha/ny;
  else if (ny < 1e-10)
    area = alpha/nx;
  else {
    v = ( alpha*alpha);

    a = alpha - nx;
    if (a > 0.)
      v -= a*a;

    a = alpha - ny;
    if (a > 0.)
      v -= a*a;

    area = v/(2.*nx*ny);
  }

  return ( area < 0. ? 0. : area > 1. ? 1. : area);
}



double plane_volume (coord n, double alpha)
{
  double al = alpha + (n.x + n.y + n.z)/2. +
    ( 0. > (-n.x) ? 0. : (-n.x)) + ( 0. > (-n.y) ? 0. : (-n.y)) + ( 0. > (-n.z) ? 0. : (-n.z));
  if (al <= 0.)
    return 0.;
  double tmp = fabs(n.x) + fabs(n.y) + fabs(n.z);
  if (al >= tmp)
    return 1.;
  if (tmp < 1e-10)
    return 0.;
  double n1 = fabs(n.x)/tmp;
  double n2 = fabs(n.y)/tmp;
  double n3 = fabs(n.z)/tmp;
  al = ( 0. > (( 1. < (al/tmp) ? 1. : (al/tmp))) ? 0. : (( 1. < (al/tmp) ? 1. : (al/tmp))));
  double al0 = ( al < (1. - al) ? al : (1. - al));
  double b1 = ( n1 < n2 ? n1 : n2);
  double b3 = ( n1 > n2 ? n1 : n2);
  double b2 = n3;
  if (b2 < b1) {
    tmp = b1;
    b1 = b2;
    b2 = tmp;
  }
  else if (b2 > b3) {
    tmp = b3;
    b3 = b2;
    b2 = tmp;
  }
  double b12 = b1 + b2;
  double bm = ( b12 < b3 ? b12 : b3);
  double pr = ( (6.*b1*b2*b3) > 1e-50 ? (6.*b1*b2*b3) : 1e-50);
  if (al0 < b1)
    tmp = al0*al0*al0/pr;
  else if (al0 < b2)
    tmp = 0.5*al0*(al0 - b1)/(b2*b3) + b1*b1*b1/pr;
  else if (al0 < bm)
    tmp = (al0*al0*(3.*b12 - al0) + b1*b1*(b1 - 3.*al0) +
    b2*b2*(b2 - 3.*al0))/pr;
  else if (b12 < b3)
    tmp = (al0 - 0.5*bm)/b3;
  else
    tmp = (al0*al0*(3. - 2.*al0) + b1*b1*(b1 - 3.*al0) +
    b2*b2*(b2 - 3.*al0) + b3*b3*(b3 - 3.*al0))/pr;

  double volume = al <= 0.5 ? tmp : 1. - tmp;
  return ( volume < 0. ? 0. : volume > 1. ? 1. : volume);
}
#line 267 "/home/pwachara/basilisk/src/geometry.h"
double rectangle_fraction (coord n, double alpha, coord a, coord b)
{
  coord n1;
   {
    alpha -= n.x*(b.x + a.x)/2.;
    n1.x = n.x*(b.x - a.x);
  } 
#line 270
{
    alpha -= n.y*(b.y + a.y)/2.;
    n1.y = n.y*(b.y - a.y);
  } 
#line 270
{
    alpha -= n.z*(b.z + a.z)/2.;
    n1.z = n.z*(b.z - a.z);
  }
  return plane_volume (n1, alpha);
}
#line 307 "/home/pwachara/basilisk/src/geometry.h"
static coord cube_edge[12][2] = {
  {{0.,0.,0.},{1.,0.,0.}},{{0.,0.,1.},{1.,0.,1.}},
  {{0.,1.,1.},{1.,1.,1.}},{{0.,1.,0.},{1.,1.,0.}},
  {{0.,0.,0.},{0.,1.,0.}},{{0.,0.,1.},{0.,1.,1.}},
  {{1.,0.,1.},{1.,1.,1.}},{{1.,0.,0.},{1.,1.,0.}},
  {{0.,0.,0.},{0.,0.,1.}},{{1.,0.,0.},{1.,0.,1.}},
  {{1.,1.,0.},{1.,1.,1.}},{{0.,1.,0.},{0.,1.,1.}}
};




static int cube_connect[12][2][4] = {
  {{9, 1, 8}, {4, 3, 7}},
  {{6, 2, 5}, {8, 0, 9}},
  {{10, 3, 11}, {5, 1, 6}},
  {{7, 0, 4}, {11, 2, 10}},
  {{3, 7, 0}, {8, 5, 11}},
  {{11, 4, 8}, {1, 6, 2}},
  {{2, 5, 1}, {9, 7, 10}},
  {{10, 6, 9}, {0, 4, 3}},
  {{5, 11, 4}, {0, 9, 1}},
  {{1, 8, 0}, {7, 10, 6}},
  {{6, 9, 7}, {3, 11, 2}},
  {{2, 10, 3}, {4, 8, 5}}
};

int facets (coord n, double alpha, coord v[12], double h)
{
  coord a[12];
  int orient[12];

  for (int i = 0; i < 12; i++) {
    coord e, d;
    double den = 0., t = alpha;
     {
      d.x = h*(cube_edge[i][0].x - 0.5);
      e.x = h*(cube_edge[i][1].x - 0.5);
      den += n.x*(e.x - d.x);
      t -= n.x*d.x;
    } 
#line 342
{
      d.y = h*(cube_edge[i][0].y - 0.5);
      e.y = h*(cube_edge[i][1].y - 0.5);
      den += n.y*(e.y - d.y);
      t -= n.y*d.y;
    } 
#line 342
{
      d.z = h*(cube_edge[i][0].z - 0.5);
      e.z = h*(cube_edge[i][1].z - 0.5);
      den += n.z*(e.z - d.z);
      t -= n.z*d.z;
    }
    orient[i] = -1;
    if (fabs (den) > 1e-10) {
      t /= den;
      if (t >= 0. && t < 1.) {
 double s = - alpha;
  {
   a[i].x = d.x + t*(e.x - d.x);
   s += n.x*e.x;
 } 
#line 353
{
   a[i].y = d.y + t*(e.y - d.y);
   s += n.y*e.y;
 } 
#line 353
{
   a[i].z = d.z + t*(e.z - d.z);
   s += n.z*e.z;
 }
 orient[i] = (s > 0.);
      }
    }
  }

  for (int i = 0; i < 12; i++) {
    int nv = 0, e = i;
    while (orient[e] >= 0) {
      int m = 0, * ne = cube_connect[e][orient[e]];
      v[nv++] = a[e];
      orient[e] = -1;
      while (m < 3 && orient[e] < 0)
 e = ne[m++];
    }
    if (nv > 2)
      return nv;
  }
  return 0;
}






double line_length_center (coord m, double alpha, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }    
#line 388
if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }

  p->x = p->y = p->z = 0.;

  if (alpha <= 0. || alpha >= n.x + n.y)
    return 0.;

  
    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }    
#line 399
if (n.y < 1e-4) {
      p->y = 0.;
      p->x = (m.x < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }

  if (alpha >= n.x) {
    p->x += 1.;
    p->y += (alpha - n.x)/n.y;
  }
  else
    p->x += alpha/n.x;

  double ax = p->x, ay = p->y;
  if (alpha >= n.y) {
    p->y += 1.;
    ay -= 1.;
    p->x += (alpha - n.y)/n.x;
    ax -= (alpha - n.y)/n.x;
  }
  else {
    p->y += alpha/n.y;
    ay -= alpha/n.y;
  }

   {
    p->x /= 2.;
    p->x = ( (p->x) < 0. ? 0. : (p->x) > 1. ? 1. : (p->x));
    if (m.x < 0.)
      p->x = 1. - p->x;
    p->x -= 0.5;
  } 
#line 424
{
    p->y /= 2.;
    p->y = ( (p->y) < 0. ? 0. : (p->y) > 1. ? 1. : (p->y));
    if (m.y < 0.)
      p->y = 1. - p->y;
    p->y -= 0.5;
  }

  return sqrt (ax*ax + ay*ay);
}




double plane_area_center (coord m, double alpha, coord * p)
{
  
    if (fabs (m.x) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.y;
      ((double *)&n)[1] = m.z;
      double length = line_length_center (n, alpha, &q);
      p->x = 0.;
      p->y = ((double *)&q)[0];
      p->z = ((double *)&q)[1];
      return length;
    }    
#line 441
if (fabs (m.y) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.z;
      ((double *)&n)[1] = m.x;
      double length = line_length_center (n, alpha, &q);
      p->y = 0.;
      p->z = ((double *)&q)[0];
      p->x = ((double *)&q)[1];
      return length;
    }    
#line 441
if (fabs (m.z) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.x;
      ((double *)&n)[1] = m.y;
      double length = line_length_center (n, alpha, &q);
      p->z = 0.;
      p->x = ((double *)&q)[0];
      p->y = ((double *)&q)[1];
      return length;
    }

  alpha += (m.x + m.y + m.z)/2.;
  coord n = m;
  
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }    
#line 455
if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }    
#line 455
if (n.z < 0.) {
      alpha -= n.z;
      n.z = - n.z;
    }

  double amax = n.x + n.y + n.z;
  if (alpha < 0. || alpha > amax) {
    p->x = p->y = p->z = 0.;
    return 0.;
  }

  double area = ( alpha*alpha);
  p->x = p->y = p->z = area*alpha;

   {
    double b = alpha - n.x;
    if (b > 0.) {
      area -= b*b;
      p->x -= b*b*(2.*n.x + alpha);
      p->y -= b*b*b;
      p->z -= b*b*b;
    }
  } 
#line 469
{
    double b = alpha - n.y;
    if (b > 0.) {
      area -= b*b;
      p->y -= b*b*(2.*n.y + alpha);
      p->z -= b*b*b;
      p->x -= b*b*b;
    }
  } 
#line 469
{
    double b = alpha - n.z;
    if (b > 0.) {
      area -= b*b;
      p->z -= b*b*(2.*n.z + alpha);
      p->x -= b*b*b;
      p->y -= b*b*b;
    }
  }

  amax = alpha - amax;
   {
    double b = amax + n.x;
    if (b > 0.) {
      area += b*b;
      p->y += b*b*(2.*n.y + alpha - n.z);
      p->z += b*b*(2.*n.z + alpha - n.y);
      p->x += b*b*b;
    }
  } 
#line 480
{
    double b = amax + n.y;
    if (b > 0.) {
      area += b*b;
      p->z += b*b*(2.*n.z + alpha - n.x);
      p->x += b*b*(2.*n.x + alpha - n.z);
      p->y += b*b*b;
    }
  } 
#line 480
{
    double b = amax + n.z;
    if (b > 0.) {
      area += b*b;
      p->x += b*b*(2.*n.x + alpha - n.y);
      p->y += b*b*(2.*n.y + alpha - n.x);
      p->z += b*b*b;
    }
  }

  area *= 3.;
   {
    if (area) {
      p->x /= area*n.x;
      p->x = ( (p->x) < 0. ? 0. : (p->x) > 1. ? 1. : (p->x));
    }
    else
      p->x = 0.;
    if (m.x < 0.) p->x = 1. - p->x;
    p->x -= 0.5;
  } 
#line 491
{
    if (area) {
      p->y /= area*n.y;
      p->y = ( (p->y) < 0. ? 0. : (p->y) > 1. ? 1. : (p->y));
    }
    else
      p->y = 0.;
    if (m.y < 0.) p->y = 1. - p->y;
    p->y -= 0.5;
  } 
#line 491
{
    if (area) {
      p->z /= area*n.z;
      p->z = ( (p->z) < 0. ? 0. : (p->z) > 1. ? 1. : (p->z));
    }
    else
      p->z = 0.;
    if (m.z < 0.) p->z = 1. - p->z;
    p->z -= 0.5;
  }

  return area*sqrt (1./(( (n.x)*(n.x))*( (n.y)*(n.y))) +
      1./(( (n.x)*(n.x))*( (n.z)*(n.z))) +
      1./(( (n.z)*(n.z))*( (n.y)*(n.y))))/6.;
}






void line_center (coord m, double alpha, double a, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }    
#line 518
if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }

  p->z = 0.;
  if (alpha <= 0.) {
    p->x = p->y = -0.5;
    return;
  }

  if (alpha >= n.x + n.y) {
    p->x = p->y = 0.;
    return;
  }

  
    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = ( (int)((m.y) > 0 ? 1 : -1))*(a/2. - 0.5);
      return;
    }    
#line 535
if (n.y < 1e-4) {
      p->y = 0.;
      p->x = ( (int)((m.x) > 0 ? 1 : -1))*(a/2. - 0.5);
      return;
    }

  p->x = p->y = ( alpha*alpha*alpha);

   {
    double b = alpha - n.x;
    if (b > 0.) {
      p->x -= ( b*b)*(alpha + 2.*n.x);
      p->y -= ( b*b*b);
    }
  } 
#line 543
{
    double b = alpha - n.y;
    if (b > 0.) {
      p->y -= ( b*b)*(alpha + 2.*n.y);
      p->x -= ( b*b*b);
    }
  }

   {
    p->x /= 6.*( (n.x)*(n.x))*n.y*a;
    p->x = ( (int)((m.x) > 0 ? 1 : -1))*(p->x - 0.5);
  } 
#line 551
{
    p->y /= 6.*( (n.y)*(n.y))*n.x*a;
    p->y = ( (int)((m.y) > 0 ? 1 : -1))*(p->y - 0.5);
  }
}
#line 564 "/home/pwachara/basilisk/src/geometry.h"
void plane_center (coord m, double alpha, double a, coord * p)
{
  
    if (fabs (m.x) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.y;
      ((double *)&n)[1] = m.z;
      line_center (n, alpha, a, &q);
      p->x = 0.;
      p->y = ((double *)&q)[0];
      p->z = ((double *)&q)[1];
      return;
    }    
#line 567
if (fabs (m.y) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.z;
      ((double *)&n)[1] = m.x;
      line_center (n, alpha, a, &q);
      p->y = 0.;
      p->z = ((double *)&q)[0];
      p->x = ((double *)&q)[1];
      return;
    }    
#line 567
if (fabs (m.z) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.x;
      ((double *)&n)[1] = m.y;
      line_center (n, alpha, a, &q);
      p->z = 0.;
      p->x = ((double *)&q)[0];
      p->y = ((double *)&q)[1];
      return;
    }

  alpha += (m.x + m.y + m.z)/2.;
  coord n = m;
  
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }    
#line 581
if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }    
#line 581
if (n.z < 0.) {
      alpha -= n.z;
      n.z = - n.z;
    }

  if (alpha <= 0. || a == 0.) {
    p->x = p->y = p->z = -0.5;
    return;
  }

  if (alpha >= n.x + n.y + n.z || a == 1.) {
    p->x = p->y = p->z = 0.;
    return;
  }

  p->x = p->y = p->z = ( (( alpha*alpha))*(( alpha*alpha)));
   {
    double b = alpha - n.x;
    if (b > 0.) {
      p->x -= ( b*b*b)*(3.*n.x + alpha);
      p->y -= ( (( b*b))*(( b*b)));
      p->z -= ( (( b*b))*(( b*b)));
    }
  } 
#line 597
{
    double b = alpha - n.y;
    if (b > 0.) {
      p->y -= ( b*b*b)*(3.*n.y + alpha);
      p->z -= ( (( b*b))*(( b*b)));
      p->x -= ( (( b*b))*(( b*b)));
    }
  } 
#line 597
{
    double b = alpha - n.z;
    if (b > 0.) {
      p->z -= ( b*b*b)*(3.*n.z + alpha);
      p->x -= ( (( b*b))*(( b*b)));
      p->y -= ( (( b*b))*(( b*b)));
    }
  }

  double amax = alpha - (n.x + n.y + n.z);
   {
    double b = amax + n.z;
    if (b > 0.) {
      p->x += ( b*b*b)*(3.*n.x + alpha - n.y);
      p->y += ( b*b*b)*(3.*n.y + alpha - n.x);
      p->z += ( (( b*b))*(( b*b)));
    }
  } 
#line 607
{
    double b = amax + n.x;
    if (b > 0.) {
      p->y += ( b*b*b)*(3.*n.y + alpha - n.z);
      p->z += ( b*b*b)*(3.*n.z + alpha - n.y);
      p->x += ( (( b*b))*(( b*b)));
    }
  } 
#line 607
{
    double b = amax + n.y;
    if (b > 0.) {
      p->z += ( b*b*b)*(3.*n.z + alpha - n.x);
      p->x += ( b*b*b)*(3.*n.x + alpha - n.z);
      p->y += ( (( b*b))*(( b*b)));
    }
  }

  double b = 24.*n.x*n.y*n.z*a;
   {
    p->x /= b*n.x;
    p->x = ( (int)((m.x) > 0 ? 1 : -1))*(p->x - 0.5);
  } 
#line 617
{
    p->y /= b*n.y;
    p->y = ( (int)((m.y) > 0 ? 1 : -1))*(p->y - 0.5);
  } 
#line 617
{
    p->z /= b*n.z;
    p->z = ( (int)((m.z) > 0 ? 1 : -1))*(p->z - 0.5);
  }
}
#line 13 "/home/pwachara/basilisk/src/fractions.h"







#line 1 "myc.h"
#line 1 "/home/pwachara/basilisk/src/myc.h"
#line 16 "/home/pwachara/basilisk/src/myc.h"
coord mycs (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 18 "/home/pwachara/basilisk/src/myc.h"
double m1,m2,m[4][3],t0,t1,t2;
  int cn;



  m1 = val(c,-1,0,-1) + val(c,-1,0,1) + val(c,-1,-1,0) + val(c,-1,1,0) +
       val(c,-1,0,0);
  m2 = val(c,1,0,-1) + val(c,1,0,1) + val(c,1,-1,0) + val(c,1,1,0) +
       val(c,1,0,0);
  m[0][0] = m1 > m2 ? 1. : -1.;

  m1 = val(c,-1,-1,0)+ val(c,1,-1,0)+ val(c,0,-1,0);
  m2 = val(c,-1,1,0)+ val(c,1,1,0)+ val(c,0,1,0);
  m[0][1] = 0.5*(m1-m2);

  m1 = val(c,-1,0,-1)+ val(c,1,0,-1)+ val(c,0,0,-1);
  m2 = val(c,-1,0,1)+ val(c,1,0,1)+ val(c,0,0,1);
  m[0][2] = 0.5*(m1-m2);



  m1 = val(c,-1,-1,0) + val(c,-1,1,0) + val(c,-1,0,0);
  m2 = val(c,1,-1,0) + val(c,1,1,0) + val(c,1,0,0);
  m[1][0] = 0.5*(m1-m2);

  m1 = val(c,0,-1,-1) + val(c,0,-1,1) + val(c,1,-1,0) + val(c,-1,-1,0) +
       val(c,0,-1,0);
  m2 = val(c,0,1,-1) + val(c,0,1,1) + val(c,1,1,0) + val(c,-1,1,0) +
       val(c,0,1,0);
  m[1][1] = m1 > m2 ? 1. : -1.;

  m1 = val(c,0,-1,-1)+ val(c,0,0,-1)+ val(c,0,1,-1);
  m2 = val(c,0,-1,1)+ val(c,0,0,1)+ val(c,0,1,1);
  m[1][2] = 0.5*(m1-m2);




  m1 = val(c,-1,0,-1)+ val(c,-1,0,1)+ val(c,-1,0,0);
  m2 = val(c,1,0,-1)+ val(c,1,0,1)+ val(c,1,0,0);
  m[2][0] = 0.5*(m1-m2);

  m1 = val(c,0,-1,-1)+ val(c,0,-1,1)+ val(c,0,-1,0);
  m2 = val(c,0,1,-1)+ val(c,0,1,1)+ val(c,0,1,0);
  m[2][1] = 0.5*(m1-m2);

  m1 = val(c,-1,0,-1) + val(c,1,0,-1) + val(c,0,-1,-1) + val(c,0,1,-1) +
       val(c,0,0,-1);
  m2 = val(c,-1,0,1) + val(c,1,0,1) + val(c,0,-1,1) + val(c,0,1,1) +
       val(c,0,0,1);
  m[2][2] = m1 > m2 ? 1. : -1.;


  t0 = fabs(m[0][0]) + fabs(m[0][1]) + fabs(m[0][2]);
  m[0][0] /= t0;
  m[0][1] /= t0;
  m[0][2] /= t0;

  t0 = fabs(m[1][0]) + fabs(m[1][1]) + fabs(m[1][2]);
  m[1][0] /= t0;
  m[1][1] /= t0;
  m[1][2] /= t0;

  t0 = fabs(m[2][0]) + fabs(m[2][1]) + fabs(m[2][2]);
  m[2][0] /= t0;
  m[2][1] /= t0;
  m[2][2] /= t0;


  t0 = fabs(m[0][0]);
  t1 = fabs(m[1][1]);
  t2 = fabs(m[2][2]);

  cn = 0;
  if (t1 > t0) {
    t0 = t1;
    cn = 1;
  }
  if (t2 > t0)
    cn = 2;


  m1 = val(c,-1,-1,-1) + val(c,-1,1,-1) + val(c,-1,-1,1) + val(c,-1,1,1) +
       2.*(val(c,-1,-1,0) + val(c,-1,1,0) + val(c,-1,0,-1) + val(c,-1,0,1)) +
       4.*val(c,-1,0,0);
  m2 = val(c,1,-1,-1) + val(c,1,1,-1) + val(c,1,-1,1) + val(c,1,1,1) +
       2.*(val(c,1,-1,0) + val(c,1,1,0) + val(c,1,0,-1) + val(c,1,0,1)) +
       4.*val(c,1,0,0);
  m[3][0] = m1 - m2;

  m1 = val(c,-1,-1,-1) + val(c,-1,-1,1) + val(c,1,-1,-1) + val(c,1,-1,1) +
       2.*( val(c,-1,-1,0) + val(c,1,-1,0) + val(c,0,-1,-1) + val(c,0,-1,1)) +
       4.*val(c,0,-1,0);
  m2 = val(c,-1,1,-1) + val(c,-1,1,1) + val(c,1,1,-1) + val(c,1,1,1) +
       2.*(val(c,-1,1,0) + val(c,1,1,0) + val(c,0,1,-1) + val(c,0,1,1)) +
       4.*val(c,0,1,0);
  m[3][1] = m1 - m2;

  m1 = val(c,-1,-1,-1) + val(c,-1,1,-1) + val(c,1,-1,-1) + val(c,1,1,-1) +
       2.*(val(c,-1,0,-1) + val(c,1,0,-1) + val(c,0,-1,-1) + val(c,0,1,-1)) +
       4.*val(c,0,0,-1);
  m2 = val(c,-1,-1,1) + val(c,-1,1,1) + val(c,1,-1,1) + val(c,1,1,1) +
       2.*(val(c,-1,0,1) + val(c,1,0,1) + val(c,0,-1,1) + val(c,0,1,1)) +
       4.*val(c,0,0,1);
  m[3][2] = m1 - m2;


  t0 = fabs(m[3][0]) + fabs(m[3][1]) + fabs(m[3][2]);
  if (t0 < 1e-30) {
    coord mxyz = {1., 0., 0.};
    return mxyz;
  }

  m[3][0] /= t0;
  m[3][1] /= t0;
  m[3][2] /= t0;


  t0 = fabs (m[3][0]);
  t1 = fabs (m[3][1]);
  t2 = fabs (m[3][2]);
  if (t1 > t0)
    t0 = t1;
  if (t2 > t0)
    t0 = t2;

  if (fabs(m[cn][cn]) > t0)
    cn = 3;


  coord mxyz = {m[cn][0], m[cn][1], m[cn][2]};
  return mxyz;
}
#line 13 "/home/pwachara/basilisk/src/fractions.h"







#line 1 "myc.h"
#line 1 "/home/pwachara/basilisk/src/myc.h"
#line 16 "/home/pwachara/basilisk/src/myc.h"
static void _stencil_mycs (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;   
#line 23 "/home/pwachara/basilisk/src/myc.h"
_stencil_val(c,-1,0,-1); _stencil_val(c,-1,0,1); _stencil_val(c,-1,-1,0); _stencil_val(c,-1,1,0);
       _stencil_val(c,-1,0,0);       
_stencil_val(c,1,0,-1); _stencil_val(c,1,0,1); _stencil_val(c,1,-1,0); _stencil_val(c,1,1,0);
       _stencil_val(c,1,0,0);  
      
          
_stencil_val(c,-1,-1,0); _stencil_val(c,1,-1,0); _stencil_val(c,0,-1,0);    
_stencil_val(c,-1,1,0); _stencil_val(c,1,1,0); _stencil_val(c,0,1,0);
   
     
_stencil_val(c,-1,0,-1); _stencil_val(c,1,0,-1); _stencil_val(c,0,0,-1);    
_stencil_val(c,-1,0,1); _stencil_val(c,1,0,1); _stencil_val(c,0,0,1);
   
     


_stencil_val(c,-1,-1,0); _stencil_val(c,-1,1,0); _stencil_val(c,-1,0,0);      
_stencil_val(c,1,-1,0); _stencil_val(c,1,1,0); _stencil_val(c,1,0,0);
     
     
_stencil_val(c,0,-1,-1); _stencil_val(c,0,-1,1); _stencil_val(c,1,-1,0); _stencil_val(c,-1,-1,0);
       _stencil_val(c,0,-1,0);        
_stencil_val(c,0,1,-1); _stencil_val(c,0,1,1); _stencil_val(c,1,1,0); _stencil_val(c,-1,1,0);
       _stencil_val(c,0,1,0);
       
           
_stencil_val(c,0,-1,-1); _stencil_val(c,0,0,-1); _stencil_val(c,0,1,-1);    
_stencil_val(c,0,-1,1); _stencil_val(c,0,0,1); _stencil_val(c,0,1,1);
   
     



_stencil_val(c,-1,0,-1); _stencil_val(c,-1,0,1); _stencil_val(c,-1,0,0);    
_stencil_val(c,1,0,-1); _stencil_val(c,1,0,1); _stencil_val(c,1,0,0);
   
     
_stencil_val(c,0,-1,-1); _stencil_val(c,0,-1,1); _stencil_val(c,0,-1,0);    
_stencil_val(c,0,1,-1); _stencil_val(c,0,1,1); _stencil_val(c,0,1,0);
   
     
_stencil_val(c,-1,0,-1); _stencil_val(c,1,0,-1); _stencil_val(c,0,-1,-1); _stencil_val(c,0,1,-1);
       _stencil_val(c,0,0,-1);        
_stencil_val(c,-1,0,1); _stencil_val(c,1,0,1); _stencil_val(c,0,-1,1); _stencil_val(c,0,1,1);
       _stencil_val(c,0,0,1);    
       
          


       
    
    
    

        
    
    
    

        
    
    
    


    
   
   

    
      
     
   
      
     
      

_stencil_val(c,-1,-1,-1); _stencil_val(c,-1,1,-1); _stencil_val(c,-1,-1,1); _stencil_val(c,-1,1,1);
_stencil_val(c,-1,-1,0); _stencil_val(c,-1,1,0); _stencil_val(c,-1,0,-1); _stencil_val(c,-1,0,1);
_stencil_val(c,-1,0,0);        
_stencil_val(c,1,-1,-1); _stencil_val(c,1,1,-1); _stencil_val(c,1,-1,1); _stencil_val(c,1,1,1);
_stencil_val(c,1,-1,0); _stencil_val(c,1,1,0); _stencil_val(c,1,0,-1); _stencil_val(c,1,0,1);
_stencil_val(c,1,0,0);       


_stencil_val(c,-1,-1,-1); _stencil_val(c,-1,-1,1); _stencil_val(c,1,-1,-1); _stencil_val(c,1,-1,1); 
_stencil_val(c,-1,-1,0); _stencil_val(c,1,-1,0); _stencil_val(c,0,-1,-1); _stencil_val(c,0,-1,1);
_stencil_val(c,0,-1,0);        
_stencil_val(c,-1,1,-1); _stencil_val(c,-1,1,1); _stencil_val(c,1,1,-1); _stencil_val(c,1,1,1);
_stencil_val(c,-1,1,0); _stencil_val(c,1,1,0); _stencil_val(c,0,1,-1); _stencil_val(c,0,1,1);
_stencil_val(c,0,1,0);       


_stencil_val(c,-1,-1,-1); _stencil_val(c,-1,1,-1); _stencil_val(c,1,-1,-1); _stencil_val(c,1,1,-1);
_stencil_val(c,-1,0,-1); _stencil_val(c,1,0,-1); _stencil_val(c,0,-1,-1); _stencil_val(c,0,1,-1);
_stencil_val(c,0,0,-1);        
_stencil_val(c,-1,-1,1); _stencil_val(c,-1,1,1); _stencil_val(c,1,-1,1); _stencil_val(c,1,1,1);
_stencil_val(c,-1,0,1); _stencil_val(c,1,0,1); _stencil_val(c,0,-1,1); _stencil_val(c,0,1,1);
_stencil_val(c,0,0,1);  
#line 149
return ;
}
#line 21 "/home/pwachara/basilisk/src/fractions.h"
#line 120 "/home/pwachara/basilisk/src/fractions.h"
     
void fractions (scalar Phi, scalar c,
  vector s, double val)
{tracing("fractions","/home/pwachara/basilisk/src/fractions.h",121);

  vector   as=(s).x.i>0?(s):new_face_vector("as");
#line 134 "/home/pwachara/basilisk/src/fractions.h"
  vector  p=new_vector("p");
  
  
  
#line 146 "/home/pwachara/basilisk/src/fractions.h"
  
#line 314 "/home/pwachara/basilisk/src/grid/stencils.h"
{
    static int _first = 1.;
    ForeachData _loop = {
      .fname = "/home/pwachara/basilisk/src/grid/multigrid.h", .line = 927, .first = _first
    };
    if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
 _attribute[s.i].input = _attribute[s.i].output = _attribute[s.i].nowarning = false;
 _attribute[s.i].width = 0;
      }
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point); 
#line 335
{
    _loop.vertex = true; 
#line 927 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
    struct { int x, y, z; } _a = {point.i, point.j, point.k};
    
      if (_a.x < point.n.x + 2) 
#line 146 "/home/pwachara/basilisk/src/fractions.h"
{





_stencil_val(Phi,0,0,0);_stencil_val(Phi,1,0,0);{ {






_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,0,0); _stencil_val(Phi,1,0,0);      _stencil_val_a(p.x,0,0,0);
_stencil_val(Phi,0,0,0);
 { _stencil_val(p.x,0,0,0);_stencil_val_a(p.x,0,0,0);   }    
}
      








{_stencil_val(Phi,0,0,0); _stencil_val(Phi,1,0,0);_stencil_val_a(p.x,0,0,0);       }}





           
#line 171 "/home/pwachara/basilisk/src/fractions.h"
    
  
}
      
#line 930 "/home/pwachara/basilisk/src/grid/multigrid.h"
if (_a.y < point.n.y + 2) 
#line 146 "/home/pwachara/basilisk/src/fractions.h"
{





_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,1,0);{ {






_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,0,0); _stencil_val(Phi,0,1,0);      _stencil_val_a(p.y,0,0,0);
_stencil_val(Phi,0,0,0);
 { _stencil_val(p.y,0,0,0);_stencil_val_a(p.y,0,0,0);   }    
}
      








{_stencil_val(Phi,0,0,0); _stencil_val(Phi,0,1,0);_stencil_val_a(p.y,0,0,0);       }}





           
#line 171 "/home/pwachara/basilisk/src/fractions.h"
    
  
}
      
#line 930 "/home/pwachara/basilisk/src/grid/multigrid.h"
if (_a.z < point.n.z + 2) 
#line 146 "/home/pwachara/basilisk/src/fractions.h"
{





_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,0,1);{ {






_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,0,0); _stencil_val(Phi,0,0,1);      _stencil_val_a(p.z,0,0,0);
_stencil_val(Phi,0,0,0);
 { _stencil_val(p.z,0,0,0);_stencil_val_a(p.z,0,0,0);   }    
}
      








{_stencil_val(Phi,0,0,0); _stencil_val(Phi,0,0,1);_stencil_val_a(p.z,0,0,0);       }}





           
#line 171 "/home/pwachara/basilisk/src/fractions.h"
    
  
}
  
#line 932 "/home/pwachara/basilisk/src/grid/multigrid.h"
}  
#line 338 "/home/pwachara/basilisk/src/grid/stencils.h"
}    
#line 328
check_stencil (&_loop);
    boundary_stencil (&_loop);
    _first = 0;
  } 
#line 912 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
#line 300
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
      for (_k = 2; _k <= point.n.x + 2; _k++) {
 point.i = _k;

 for (point.j = 2; point.j <= point.n.y + 2; point.j++)

   for (point.k = 2; point.k <= point.n.z + 2; point.k++) 
#line 913
{
    int ig = -1; NOT_UNUSED(ig);

    int jg = -1; NOT_UNUSED(jg);


    int kg = -1; NOT_UNUSED(kg); 







{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 927
{
    struct { int x, y, z; } _a = {point.i, point.j, point.k};
    
      if (_a.x < point.n.x + 2) 
#line 146 "/home/pwachara/basilisk/src/fractions.h"
{





    if ((val(Phi,0,0,0) - val)*(val(Phi,1,0,0) - val) < 0.) {






      val(p.x,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,1,0,0));
      if (val(Phi,0,0,0) < val)
 val(p.x,0,0,0) = 1. - val(p.x,0,0,0);
    }
#line 171 "/home/pwachara/basilisk/src/fractions.h"
    else
      val(p.x,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,1,0,0) > val);
  }
      
#line 930 "/home/pwachara/basilisk/src/grid/multigrid.h"
if (_a.y < point.n.y + 2) 
#line 146 "/home/pwachara/basilisk/src/fractions.h"
{





    if ((val(Phi,0,0,0) - val)*(val(Phi,0,1,0) - val) < 0.) {






      val(p.y,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,0,1,0));
      if (val(Phi,0,0,0) < val)
 val(p.y,0,0,0) = 1. - val(p.y,0,0,0);
    }
#line 171 "/home/pwachara/basilisk/src/fractions.h"
    else
      val(p.y,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,0,1,0) > val);
  }
      
#line 930 "/home/pwachara/basilisk/src/grid/multigrid.h"
if (_a.z < point.n.z + 2) 
#line 146 "/home/pwachara/basilisk/src/fractions.h"
{





    if ((val(Phi,0,0,0) - val)*(val(Phi,0,0,1) - val) < 0.) {






      val(p.z,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,0,0,1));
      if (val(Phi,0,0,0) < val)
 val(p.z,0,0,0) = 1. - val(p.z,0,0,0);
    }
#line 171 "/home/pwachara/basilisk/src/fractions.h"
    else
      val(p.z,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,0,0,1) > val);
  }
  
#line 932 "/home/pwachara/basilisk/src/grid/multigrid.h"
}}  
#line 922
}      
#line 317
}
  }
}

#line 923
}
#line 190 "/home/pwachara/basilisk/src/fractions.h"
  
    _attribute[p.x.i].dirty = false;    _attribute[p.y.i].dirty = false;    _attribute[p.z.i].dirty = false;

  scalar s_x = as.x, s_y = as.y, s_z = as.z;
  
  
#line 314 "/home/pwachara/basilisk/src/grid/stencils.h"
{
    static int _first = 1.;
    ForeachData _loop = {
      .fname = "/home/pwachara/basilisk/src/fractions.h", .line = 194, .first = _first
    };
    if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
 _attribute[s.i].input = _attribute[s.i].output = _attribute[s.i].nowarning = false;
 _attribute[s.i].width = 0;
      }
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point);
#line 194 "/home/pwachara/basilisk/src/fractions.h"
{ 
#line 379 "/home/pwachara/basilisk/src/grid/stencils.h"
_loop.face |= (1 << 2);  
#line 199 "/home/pwachara/basilisk/src/fractions.h"
{    
#line 231 "/home/pwachara/basilisk/src/fractions.h"
    
    
     { 
_stencil_val(p.y,0,0,0); _stencil_val(p.y,1,0,0);    

} 
#line 233
{ 
_stencil_val(p.x,0,0,0); _stencil_val(p.x,0,1,0);    

}





{
      { _stencil_val(p.x,0,0,0);_stencil_val_a(s_z,0,0,0); } 
{      





      
   






      
      for (int i = 0; i <= 1; i++)
 {
   {_stencil_val(p.x,0,i,0); _stencil_val(p.x,0,i,0); {              
     _stencil_val(p.x,0,i,0);_stencil_val(Phi,0,i,0); 
          
     
   }      }   
#line 261
{_stencil_val(p.y,i,0,0); _stencil_val(p.y,i,0,0); {              
     _stencil_val(p.y,i,0,0);_stencil_val(Phi,i,0,0); 
          
     
   }      }}








{
 {_stencil_val(p.x,0,0,0);_stencil_val(p.y,0,0,0);_stencil_val(p.x,0,0,0);_stencil_val(p.y,0,0,0);_stencil_val_a(s_z,0,0,0);         }
{
 {_stencil_val_a(s_z,0,0,0);     } 
{

_stencil_val(p.x,0,0,0); _stencil_val(p.x,0,1,0); _stencil_val(p.y,0,0,0); _stencil_val(p.y,1,0,0); _stencil_val_a(s_z,0,0,0);       



      }}}
#line 274 "/home/pwachara/basilisk/src/fractions.h"
         
          
      
    







}}  
} 
#line 377 "/home/pwachara/basilisk/src/grid/stencils.h"
_loop.face |= (1 << 0);  
#line 199 "/home/pwachara/basilisk/src/fractions.h"
{    
#line 231 "/home/pwachara/basilisk/src/fractions.h"
    
    
     { 
_stencil_val(p.z,0,0,0); _stencil_val(p.z,0,1,0);    

} 
#line 233
{ 
_stencil_val(p.y,0,0,0); _stencil_val(p.y,0,0,1);    

}





{
      { _stencil_val(p.y,0,0,0);_stencil_val_a(s_x,0,0,0); } 
{      





      
   






      
      for (int i = 0; i <= 1; i++)
 {
   {_stencil_val(p.y,0,0,i); _stencil_val(p.y,0,0,i); {              
     _stencil_val(p.y,0,0,i);_stencil_val(Phi,0,0,i); 
          
     
   }      }   
#line 261
{_stencil_val(p.z,0,i,0); _stencil_val(p.z,0,i,0); {              
     _stencil_val(p.z,0,i,0);_stencil_val(Phi,0,i,0); 
          
     
   }      }}








{
 {_stencil_val(p.y,0,0,0);_stencil_val(p.z,0,0,0);_stencil_val(p.y,0,0,0);_stencil_val(p.z,0,0,0);_stencil_val_a(s_x,0,0,0);         }
{
 {_stencil_val_a(s_x,0,0,0);     } 
{

_stencil_val(p.y,0,0,0); _stencil_val(p.y,0,0,1); _stencil_val(p.z,0,0,0); _stencil_val(p.z,0,1,0); _stencil_val_a(s_x,0,0,0);       



      }}}
#line 274 "/home/pwachara/basilisk/src/fractions.h"
         
          
      
    







}}  
} 
#line 378 "/home/pwachara/basilisk/src/grid/stencils.h"
_loop.face |= (1 << 1);  
#line 199 "/home/pwachara/basilisk/src/fractions.h"
{    
#line 231 "/home/pwachara/basilisk/src/fractions.h"
    
    
     { 
_stencil_val(p.x,0,0,0); _stencil_val(p.x,0,0,1);    

} 
#line 233
{ 
_stencil_val(p.z,0,0,0); _stencil_val(p.z,1,0,0);    

}





{
      { _stencil_val(p.z,0,0,0);_stencil_val_a(s_y,0,0,0); } 
{      





      
   






      
      for (int i = 0; i <= 1; i++)
 {
   {_stencil_val(p.z,i,0,0); _stencil_val(p.z,i,0,0); {              
     _stencil_val(p.z,i,0,0);_stencil_val(Phi,i,0,0); 
          
     
   }      }   
#line 261
{_stencil_val(p.x,0,0,i); _stencil_val(p.x,0,0,i); {              
     _stencil_val(p.x,0,0,i);_stencil_val(Phi,0,0,i); 
          
     
   }      }}








{
 {_stencil_val(p.z,0,0,0);_stencil_val(p.x,0,0,0);_stencil_val(p.z,0,0,0);_stencil_val(p.x,0,0,0);_stencil_val_a(s_y,0,0,0);         }
{
 {_stencil_val_a(s_y,0,0,0);     } 
{

_stencil_val(p.z,0,0,0); _stencil_val(p.z,1,0,0); _stencil_val(p.x,0,0,0); _stencil_val(p.x,0,0,1); _stencil_val_a(s_y,0,0,0);       



      }}}
#line 274 "/home/pwachara/basilisk/src/fractions.h"
         
          
      
    







}}  
}}

    
#line 328 "/home/pwachara/basilisk/src/grid/stencils.h"
check_stencil (&_loop);
    boundary_stencil (&_loop);
    _first = 0;
  }
#line 300 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
      for (_k = 2; _k <= point.n.x + 2; _k++) {
 point.i = _k;

 for (point.j = 2; point.j <= point.n.y + 2; point.j++)

   for (point.k = 2; point.k <= point.n.z + 2; point.k++)
#line 194 "/home/pwachara/basilisk/src/fractions.h"
{
  
#line 400 "/home/pwachara/basilisk/src/grid/multigrid.h"
if (point.i < point.n.x + 2 && point.j < point.n.y + 2) {
    int kg = -1; NOT_UNUSED(kg);  
#line 199 "/home/pwachara/basilisk/src/fractions.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 199 "/home/pwachara/basilisk/src/fractions.h"
{
#line 231 "/home/pwachara/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
     {
      n.x = val(p.y,0,0,0) - val(p.y,1,0,0);
      nn += fabs(n.x);
    } 
#line 233
{
      n.y = val(p.x,0,0,0) - val(p.x,0,1,0);
      nn += fabs(n.y);
    }





    if (nn == 0.)
      val(s_z,0,0,0) = val(p.x,0,0,0);
    else {





      
 n.x /= nn; n.y /= nn;






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
   if (val(p.x,0,i,0) > 0. && val(p.x,0,i,0) < 1.) {
     double a = ( (int)((val(Phi,0,i,0) - val) > 0 ? 1 : -1))*(val(p.x,0,i,0) - 0.5);
     alpha += n.x*a + n.y*(i - 0.5);
     ni++;
   }   
#line 261
if (val(p.y,i,0,0) > 0. && val(p.y,i,0,0) < 1.) {
     double a = ( (int)((val(Phi,i,0,0) - val) > 0 ? 1 : -1))*(val(p.y,i,0,0) - 0.5);
     alpha += n.y*a + n.x*(i - 0.5);
     ni++;
   }}
#line 274 "/home/pwachara/basilisk/src/fractions.h"
      if (ni == 0)
 val(s_z,0,0,0) = ( (val(p.x,0,0,0)) > (val(p.y,0,0,0)) ? (val(p.x,0,0,0)) : (val(p.y,0,0,0)));
      else if (ni != 4)
 val(s_z,0,0,0) = line_area (n.x, n.y, alpha/ni);
      else {

 val(s_z,0,0,0) = (val(p.x,0,0,0) + val(p.x,0,1,0) + val(p.y,0,0,0) + val(p.y,1,0,0) > 2.);



      }
    }
  }}
  
#line 403 "/home/pwachara/basilisk/src/grid/multigrid.h"
}  
#line 386
if (point.j < point.n.y + 2 && point.k < point.n.z + 2) {
    int ig = -1; NOT_UNUSED(ig);  
#line 199 "/home/pwachara/basilisk/src/fractions.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_y = Delta;

  double Delta_z = Delta;


  double Delta_x = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_y);

  NOT_UNUSED(Delta_z);


  NOT_UNUSED(Delta_x);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 199 "/home/pwachara/basilisk/src/fractions.h"
{
#line 231 "/home/pwachara/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
     {
      n.y = val(p.z,0,0,0) - val(p.z,0,1,0);
      nn += fabs(n.y);
    } 
#line 233
{
      n.z = val(p.y,0,0,0) - val(p.y,0,0,1);
      nn += fabs(n.z);
    }





    if (nn == 0.)
      val(s_x,0,0,0) = val(p.y,0,0,0);
    else {





      
 n.y /= nn; n.z /= nn;






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
   if (val(p.y,0,0,i) > 0. && val(p.y,0,0,i) < 1.) {
     double a = ( (int)((val(Phi,0,0,i) - val) > 0 ? 1 : -1))*(val(p.y,0,0,i) - 0.5);
     alpha += n.y*a + n.z*(i - 0.5);
     ni++;
   }   
#line 261
if (val(p.z,0,i,0) > 0. && val(p.z,0,i,0) < 1.) {
     double a = ( (int)((val(Phi,0,i,0) - val) > 0 ? 1 : -1))*(val(p.z,0,i,0) - 0.5);
     alpha += n.z*a + n.y*(i - 0.5);
     ni++;
   }}
#line 274 "/home/pwachara/basilisk/src/fractions.h"
      if (ni == 0)
 val(s_x,0,0,0) = ( (val(p.y,0,0,0)) > (val(p.z,0,0,0)) ? (val(p.y,0,0,0)) : (val(p.z,0,0,0)));
      else if (ni != 4)
 val(s_x,0,0,0) = line_area (n.y, n.z, alpha/ni);
      else {

 val(s_x,0,0,0) = (val(p.y,0,0,0) + val(p.y,0,0,1) + val(p.z,0,0,0) + val(p.z,0,1,0) > 2.);



      }
    }
  }}
  
#line 389 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
  


if (point.i < point.n.x + 2 && point.k < point.n.z + 2) {
    int jg = -1; NOT_UNUSED(jg);  
#line 199 "/home/pwachara/basilisk/src/fractions.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_z = Delta;

  double Delta_x = Delta;


  double Delta_y = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_z);

  NOT_UNUSED(Delta_x);


  NOT_UNUSED(Delta_y);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 199 "/home/pwachara/basilisk/src/fractions.h"
{
#line 231 "/home/pwachara/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
     {
      n.z = val(p.x,0,0,0) - val(p.x,0,0,1);
      nn += fabs(n.z);
    } 
#line 233
{
      n.x = val(p.z,0,0,0) - val(p.z,1,0,0);
      nn += fabs(n.x);
    }





    if (nn == 0.)
      val(s_y,0,0,0) = val(p.z,0,0,0);
    else {





      
 n.z /= nn; n.x /= nn;






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
   if (val(p.z,i,0,0) > 0. && val(p.z,i,0,0) < 1.) {
     double a = ( (int)((val(Phi,i,0,0) - val) > 0 ? 1 : -1))*(val(p.z,i,0,0) - 0.5);
     alpha += n.z*a + n.x*(i - 0.5);
     ni++;
   }   
#line 261
if (val(p.x,0,0,i) > 0. && val(p.x,0,0,i) < 1.) {
     double a = ( (int)((val(Phi,0,0,i) - val) > 0 ? 1 : -1))*(val(p.x,0,0,i) - 0.5);
     alpha += n.x*a + n.z*(i - 0.5);
     ni++;
   }}
#line 274 "/home/pwachara/basilisk/src/fractions.h"
      if (ni == 0)
 val(s_y,0,0,0) = ( (val(p.z,0,0,0)) > (val(p.x,0,0,0)) ? (val(p.z,0,0,0)) : (val(p.x,0,0,0)));
      else if (ni != 4)
 val(s_y,0,0,0) = line_area (n.z, n.x, alpha/ni);
      else {

 val(s_y,0,0,0) = (val(p.z,0,0,0) + val(p.z,1,0,0) + val(p.x,0,0,0) + val(p.x,0,0,1) > 2.);



      }
    }
  }}
  
#line 396 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
#line 286 "/home/pwachara/basilisk/src/fractions.h"
}
      
#line 317 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
  }
}  
#line 314 "/home/pwachara/basilisk/src/grid/stencils.h"
{
    static int _first = 1.;
    ForeachData _loop = {
      .fname = "/home/pwachara/basilisk/src/fractions.h", .line = 294, .first = _first
    };
    if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
 _attribute[s.i].input = _attribute[s.i].output = _attribute[s.i].nowarning = false;
 _attribute[s.i].width = 0;
      }
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point); 
#line 294 "/home/pwachara/basilisk/src/fractions.h"
{    




    
    
     { 
_stencil_val(as.x,0,0,0); _stencil_val(as.x,1,0,0);    

} 
#line 301
{ 
_stencil_val(as.y,0,0,0); _stencil_val(as.y,0,1,0);    

} 
#line 301
{ 
_stencil_val(as.z,0,0,0); _stencil_val(as.z,0,0,1);    

}
{
      { _stencil_val(as.x,0,0,0);_stencil_val_a(c,0,0,0); } 
{      
      
   






      
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   {
     {_stencil_val(p.x,0,i,j); _stencil_val(p.x,0,i,j); {              
       _stencil_val(p.x,0,i,j);_stencil_val(Phi,0,i,j); 
                
       
     }      }     
#line 320
{_stencil_val(p.y,j,0,i); _stencil_val(p.y,j,0,i); {              
       _stencil_val(p.y,j,0,i);_stencil_val(Phi,j,0,i); 
                
       
     }      }     
#line 320
{_stencil_val(p.z,i,j,0); _stencil_val(p.z,i,j,0); {              
       _stencil_val(p.z,i,j,0);_stencil_val(Phi,i,j,0); 
                
       
     }      }}




{
 { _stencil_val(as.x,0,0,0);_stencil_val_a(c,0,0,0); }
{
 {_stencil_val_a(c,0,0,0);  }
 
{_stencil_val_a(c,0,0,0);    }}}    
}}  
}    
#line 328 "/home/pwachara/basilisk/src/grid/stencils.h"
check_stencil (&_loop);
    boundary_stencil (&_loop);
    _first = 0;
  } 
#line 265 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
      for (_k = 2; _k < point.n.x + 2; _k++) {
 point.i = _k;

 for (point.j = 2; point.j < point.n.y + 2; point.j++)

   for (point.k = 2; point.k < point.n.z + 2; point.k++) 
#line 294 "/home/pwachara/basilisk/src/fractions.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 294 "/home/pwachara/basilisk/src/fractions.h"
{




    coord n;
    double nn = 0.;
     {
      n.x = val(as.x,0,0,0) - val(as.x,1,0,0);
      nn += fabs(n.x);
    } 
#line 301
{
      n.y = val(as.y,0,0,0) - val(as.y,0,1,0);
      nn += fabs(n.y);
    } 
#line 301
{
      n.z = val(as.z,0,0,0) - val(as.z,0,0,1);
      nn += fabs(n.z);
    }
    if (nn == 0.)
      val(c,0,0,0) = val(as.x,0,0,0);
    else {
      
 n.x /= nn; n.y /= nn; n.z /= nn;






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   {
     if (val(p.x,0,i,j) > 0. && val(p.x,0,i,j) < 1.) {
       double a = ( (int)((val(Phi,0,i,j) - val) > 0 ? 1 : -1))*(val(p.x,0,i,j) - 0.5);
       alpha += n.x*a + n.y*(i - 0.5) + n.z*(j - 0.5);
       ni++;
     }     
#line 320
if (val(p.y,j,0,i) > 0. && val(p.y,j,0,i) < 1.) {
       double a = ( (int)((val(Phi,j,0,i) - val) > 0 ? 1 : -1))*(val(p.y,j,0,i) - 0.5);
       alpha += n.y*a + n.z*(i - 0.5) + n.x*(j - 0.5);
       ni++;
     }     
#line 320
if (val(p.z,i,j,0) > 0. && val(p.z,i,j,0) < 1.) {
       double a = ( (int)((val(Phi,i,j,0) - val) > 0 ? 1 : -1))*(val(p.z,i,j,0) - 0.5);
       alpha += n.z*a + n.x*(i - 0.5) + n.y*(j - 0.5);
       ni++;
     }}




      if (ni == 0)
 val(c,0,0,0) = val(as.x,0,0,0);
      else if (ni < 3 || ni > 6)
 val(c,0,0,0) = 0.;
      else
 val(c,0,0,0) = plane_volume (n, alpha/ni);
    }
  }}      
#line 282 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
  }
}
#line 336 "/home/pwachara/basilisk/src/fractions.h"
delete((scalar*)((vector[]){p,{{-1},{-1},{-1}}}));
#line 351 "/home/pwachara/basilisk/src/fractions.h"
end_tracing("fractions","/home/pwachara/basilisk/src/fractions.h",351);}





void macro_fraction (scalar f, double func)
{
  {
    scalar  phi=new_vertex_scalar("phi");    
#line 314 "/home/pwachara/basilisk/src/grid/stencils.h"
{
    static int _first = 1.;
    ForeachData _loop = {
      .fname = "/home/pwachara/basilisk/src/fractions.h", .line = 361, .first = _first
    };
    if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
 _attribute[s.i].input = _attribute[s.i].output = _attribute[s.i].nowarning = false;
 _attribute[s.i].width = 0;
      }
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point); 
#line 335
{
    _loop.vertex = true;
      
#line 362 "/home/pwachara/basilisk/src/fractions.h"
{_stencil_val_a(phi,0,0,0);  }  
#line 338 "/home/pwachara/basilisk/src/grid/stencils.h"
}    
#line 328
check_stencil (&_loop);
    boundary_stencil (&_loop);
    _first = 0;
  } 
#line 912 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
#line 300
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
      for (_k = 2; _k <= point.n.x + 2; _k++) {
 point.i = _k;

 for (point.j = 2; point.j <= point.n.y + 2; point.j++)

   for (point.k = 2; point.k <= point.n.z + 2; point.k++) 
#line 913
{
    int ig = -1; NOT_UNUSED(ig);

    int jg = -1; NOT_UNUSED(jg);


    int kg = -1; NOT_UNUSED(kg);      
#line 362 "/home/pwachara/basilisk/src/fractions.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 362 "/home/pwachara/basilisk/src/fractions.h"
val(phi,0,0,0) = func;}
  
#line 922 "/home/pwachara/basilisk/src/grid/multigrid.h"
}      
#line 317
}
  }
}

#line 923
}    
#line 363 "/home/pwachara/basilisk/src/fractions.h"
fractions (phi, f
#line 121
,
(  vector) {0}, 0.
#line 363
);delete((scalar*)((scalar[]){phi,{-1}}));
  }
}

void macro_solid (scalar cs, vector fs, double func)
{
  {
    scalar  phi=new_vertex_scalar("phi");    
#line 314 "/home/pwachara/basilisk/src/grid/stencils.h"
{
    static int _first = 1.;
    ForeachData _loop = {
      .fname = "/home/pwachara/basilisk/src/fractions.h", .line = 371, .first = _first
    };
    if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
 _attribute[s.i].input = _attribute[s.i].output = _attribute[s.i].nowarning = false;
 _attribute[s.i].width = 0;
      }
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point); 
#line 335
{
    _loop.vertex = true;
      
#line 372 "/home/pwachara/basilisk/src/fractions.h"
{_stencil_val_a(phi,0,0,0);  }  
#line 338 "/home/pwachara/basilisk/src/grid/stencils.h"
}    
#line 328
check_stencil (&_loop);
    boundary_stencil (&_loop);
    _first = 0;
  } 
#line 912 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
#line 300
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
      for (_k = 2; _k <= point.n.x + 2; _k++) {
 point.i = _k;

 for (point.j = 2; point.j <= point.n.y + 2; point.j++)

   for (point.k = 2; point.k <= point.n.z + 2; point.k++) 
#line 913
{
    int ig = -1; NOT_UNUSED(ig);

    int jg = -1; NOT_UNUSED(jg);


    int kg = -1; NOT_UNUSED(kg);      
#line 372 "/home/pwachara/basilisk/src/fractions.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 372 "/home/pwachara/basilisk/src/fractions.h"
val(phi,0,0,0) = func;}
  
#line 922 "/home/pwachara/basilisk/src/grid/multigrid.h"
}      
#line 317
}
  }
}

#line 923
}    
#line 373 "/home/pwachara/basilisk/src/fractions.h"
fractions (phi, cs, fs
#line 122
, 0.
#line 373
);delete((scalar*)((scalar[]){phi,{-1}}));
  }
}
#line 401 "/home/pwachara/basilisk/src/fractions.h"
coord youngs_normal (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
  
#line 403 "/home/pwachara/basilisk/src/fractions.h"
coord n;
  double nn = 0.;
  if (!(3 == 2)) qassert ("/home/pwachara/basilisk/src/fractions.h", 405, "dimension == 2");
   {
    n.x = (val(c,-1,1,0) + 2.*val(c,-1,0,0) + val(c,-1,-1,0) -
    val(c,+1,1,0) - 2.*val(c,+1,0,0) - val(c,+1,-1,0));
    nn += fabs(n.x);
  } 
#line 406
{
    n.y = (val(c,0,-1,1) + 2.*val(c,0,-1,0) + val(c,0,-1,-1) -
    val(c,0,+1,1) - 2.*val(c,0,+1,0) - val(c,0,+1,-1));
    nn += fabs(n.y);
  } 
#line 406
{
    n.z = (val(c,1,0,-1) + 2.*val(c,0,0,-1) + val(c,-1,0,-1) -
    val(c,1,0,+1) - 2.*val(c,0,0,+1) - val(c,-1,0,+1));
    nn += fabs(n.z);
  }

  if (nn > 0.)
    {
      n.x /= nn;      n.y /= nn;      n.z /= nn;}
  else
    n.x = 1.;
  return n;
}





coord facet_normal (Point point, scalar c, vector s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
  
#line 426 "/home/pwachara/basilisk/src/fractions.h"
if (s.x.i >= 0) {
    coord n;
    double nn = 0.;
     {
      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
      nn += fabs(n.x);
    } 
#line 429
{
      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
      nn += fabs(n.y);
    } 
#line 429
{
      n.z = val(s.z,0,0,0) - val(s.z,0,0,1);
      nn += fabs(n.z);
    }
    if (nn > 0.)
      {
 n.x /= nn; n.y /= nn; n.z /= nn;}
    else
      {
 n.x = 1./3; n.y = 1./3; n.z = 1./3;}
    return n;
  }
  return mycs (point, c);
}
#line 424
static void _stencil_facet_normal (Point point, scalar c, vector s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
  
#line 426 "/home/pwachara/basilisk/src/fractions.h"
if (s.x.i >= 0) {    
    
    
     { 
_stencil_val(s.x,0,0,0); _stencil_val(s.x,1,0,0);    

} 
#line 429
{ 
_stencil_val(s.y,0,0,0); _stencil_val(s.y,0,1,0);    

} 
#line 429
{ 
_stencil_val(s.z,0,0,0); _stencil_val(s.z,0,0,1);    

}
      
   
       
    
       
  
    return ;
  } 
_stencil_mycs (point, c);  return;
}
#line 451 "/home/pwachara/basilisk/src/fractions.h"
     
void reconstruction (const scalar c, vector n, scalar alpha)
{tracing("reconstruction","/home/pwachara/basilisk/src/fractions.h",452);  
#line 314 "/home/pwachara/basilisk/src/grid/stencils.h"
{
    static int _first = 1.;
    ForeachData _loop = {
      .fname = "/home/pwachara/basilisk/src/fractions.h", .line = 454, .first = _first
    };
    if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
 _attribute[s.i].input = _attribute[s.i].output = _attribute[s.i].nowarning = false;
 _attribute[s.i].width = 0;
      }
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point); 
#line 454 "/home/pwachara/basilisk/src/fractions.h"
{





_stencil_val(c,0,0,0); _stencil_val(c,0,0,0);{ {
      _stencil_val_a(alpha,0,0,0);  
      
 {_stencil_val_a(n.x,0,0,0);  } {_stencil_val_a(n.y,0,0,0);  } {_stencil_val_a(n.z,0,0,0);  }
    } 
{  






       _stencil_mycs (point, c);
      
 {_stencil_val_a(n.x,0,0,0);  } {_stencil_val_a(n.y,0,0,0);  } {_stencil_val_a(n.z,0,0,0);  }
_stencil_val(c,0,0,0);      _stencil_val_a(alpha,0,0,0);    
    }}  
}    
#line 328 "/home/pwachara/basilisk/src/grid/stencils.h"
check_stencil (&_loop);
    boundary_stencil (&_loop);
    _first = 0;
  } 
#line 265 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
      for (_k = 2; _k < point.n.x + 2; _k++) {
 point.i = _k;

 for (point.j = 2; point.j < point.n.y + 2; point.j++)

   for (point.k = 2; point.k < point.n.z + 2; point.k++) 
#line 454 "/home/pwachara/basilisk/src/fractions.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 454 "/home/pwachara/basilisk/src/fractions.h"
{





    if (val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) {
      val(alpha,0,0,0) = 0.;
      
 val(n.x,0,0,0) = 0.; val(n.y,0,0,0) = 0.; val(n.z,0,0,0) = 0.;
    }
    else {






      coord m = mycs (point, c);
      
 val(n.x,0,0,0) = m.x; val(n.y,0,0,0) = m.y; val(n.z,0,0,0) = m.z;
      val(alpha,0,0,0) = plane_alpha (val(c,0,0,0), m);
    }
  }}      
#line 282 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
  }
}
#line 495 "/home/pwachara/basilisk/src/fractions.h"
end_tracing("reconstruction","/home/pwachara/basilisk/src/fractions.h",495);}
#line 515 "/home/pwachara/basilisk/src/fractions.h"
     
void output_facets (scalar c, FILE * fp, vector s)
{tracing("output_facets","/home/pwachara/basilisk/src/fractions.h",516);
{  
#line 314 "/home/pwachara/basilisk/src/grid/stencils.h"
{
    static int _first = 1.;
    ForeachData _loop = {
      .fname = "/home/pwachara/basilisk/src/fractions.h", .line = 518, .first = _first
    };
    if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
 _attribute[s.i].input = _attribute[s.i].output = _attribute[s.i].nowarning = false;
 _attribute[s.i].width = 0;
      }
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point);
    
#line 519 "/home/pwachara/basilisk/src/fractions.h"
{_stencil_val(c,0,0,0); _stencil_val(c,0,0,0); {  
       _stencil_facet_normal (point, c, s);     
      _stencil_val(c,0,0,0);        
#line 531 "/home/pwachara/basilisk/src/fractions.h"
      
                 
              
     
 
        
 

    }        }    
#line 328 "/home/pwachara/basilisk/src/grid/stencils.h"
check_stencil (&_loop);
    boundary_stencil (&_loop);
    _first = 0;
  } 
#line 265 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
      for (_k = 2; _k < point.n.x + 2; _k++) {
 point.i = _k;

 for (point.j = 2; point.j < point.n.y + 2; point.j++)

   for (point.k = 2; point.k < point.n.z + 2; point.k++)
    
#line 519 "/home/pwachara/basilisk/src/fractions.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 519 "/home/pwachara/basilisk/src/fractions.h"
if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      double alpha = plane_alpha (val(c,0,0,0), n);
#line 531 "/home/pwachara/basilisk/src/fractions.h"
      coord v[12];
      int m = facets (n, alpha, v, 1.);
      for (int i = 0; i < m; i++)
 fprintf (fp, "%g %g %g\n",
   x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
      if (m > 0)
 fputc ('\n', fp);

    }}      
#line 282 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
  }
}
#line 539 "/home/pwachara/basilisk/src/fractions.h"
}

  fflush (fp);
end_tracing("output_facets","/home/pwachara/basilisk/src/fractions.h",542);}







     
double interface_area (scalar c)
{tracing("interface_area","/home/pwachara/basilisk/src/fractions.h",551);
  double area = 0.;  
#line 314 "/home/pwachara/basilisk/src/grid/stencils.h"
{
    static int _first = 1.;
    ForeachData _loop = {
      .fname = "/home/pwachara/basilisk/src/fractions.h", .line = 554, .first = _first
    };
    if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
 _attribute[s.i].input = _attribute[s.i].output = _attribute[s.i].nowarning = false;
 _attribute[s.i].width = 0;
      }
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point);
    
#line 555 "/home/pwachara/basilisk/src/fractions.h"
{_stencil_val(c,0,0,0); _stencil_val(c,0,0,0); {   
       _stencil_mycs (point, c);     
      _stencil_val(c,0,0,0); 
             
    }        }    
#line 328 "/home/pwachara/basilisk/src/grid/stencils.h"
check_stencil (&_loop);
    boundary_stencil (&_loop);
    _first = 0;
  } 
#line 265 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
  OMP_PARALLEL (reduction(+:area)) {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
      for (_k = 2; _k < point.n.x + 2; _k++) {
 point.i = _k;

 for (point.j = 2; point.j < point.n.y + 2; point.j++)

   for (point.k = 2; point.k < point.n.z + 2; point.k++)
    
#line 555 "/home/pwachara/basilisk/src/fractions.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 555 "/home/pwachara/basilisk/src/fractions.h"
if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = mycs (point, c), p;
      double alpha = plane_alpha (val(c,0,0,0), n);
      area += pow(Delta, 3 - 1)*plane_area_center (n, alpha, &p);
    }}      
#line 282 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
  }
}
#line 559 "/home/pwachara/basilisk/src/fractions.h"
{mpi_all_reduce_array(&area,MPI_DOUBLE,MPI_SUM,1);}
  {end_tracing("interface_area","/home/pwachara/basilisk/src/fractions.h",560);return area;}
end_tracing("interface_area","/home/pwachara/basilisk/src/fractions.h",561);}
#line 13 "/home/pwachara/basilisk/src/embed.h"






scalar  cs={0};
vector  fs={{1},{2},{3}};

double (* metric_embed_factor) (Point, coord) = NULL;
#line 95 "/home/pwachara/basilisk/src/embed.h"

static inline coord embed_face_barycentre_z (Point point, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 99 "/home/pwachara/basilisk/src/embed.h"
coord n1 = {0};
  double nn = 0.;
  scalar f = fs.z;
   {
    n1.x = (val(f,-1,-1,i) + 2.*val(f,-1,0,i) + val(f,-1,1,i) -
     val(f,+1,-1,i) - 2.*val(f,+1,0,i) - val(f,+1,1,i));
    nn += fabs(n1.x);
  } 
#line 102
{
    n1.y = (val(f,-1,-1,i) + 2.*val(f,0,-1,i) + val(f,1,-1,i) -
     val(f,-1,+1,i) - 2.*val(f,0,+1,i) - val(f,1,+1,i));
    nn += fabs(n1.y);
  }
  if (!nn)
    return (coord){0.,0.,0.};
  
    n1.x /= nn;    n1.y /= nn;

  coord n, p1, p;
  ((double *)&n)[0] = n1.x, ((double *)&n)[1] = n1.y;
  double alpha = line_alpha (val(f,0,0,i), n);
  line_center (n, alpha, val(f,0,0,i), &p1);
  p.x = ((double *)&p1)[0], p.y = ((double *)&p1)[1], p.z = 0.;
  return p;
}
#line 96
static inline coord embed_face_barycentre_x (Point point, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_y = Delta;

  double Delta_z = Delta;


  double Delta_x = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_y);

  NOT_UNUSED(Delta_z);


  NOT_UNUSED(Delta_x);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 99 "/home/pwachara/basilisk/src/embed.h"
coord n1 = {0};
  double nn = 0.;
  scalar f = fs.x;
   {
    n1.y = (val(f,i,-1,-1) + 2.*val(f,i,-1,0) + val(f,i,-1,1) -
     val(f,i,+1,-1) - 2.*val(f,i,+1,0) - val(f,i,+1,1));
    nn += fabs(n1.y);
  } 
#line 102
{
    n1.z = (val(f,i,-1,-1) + 2.*val(f,i,0,-1) + val(f,i,1,-1) -
     val(f,i,-1,+1) - 2.*val(f,i,0,+1) - val(f,i,1,+1));
    nn += fabs(n1.z);
  }
  if (!nn)
    return (coord){0.,0.,0.};
  
    n1.y /= nn;    n1.z /= nn;

  coord n, p1, p;
  ((double *)&n)[0] = n1.y, ((double *)&n)[1] = n1.z;
  double alpha = line_alpha (val(f,i,0,0), n);
  line_center (n, alpha, val(f,i,0,0), &p1);
  p.y = ((double *)&p1)[0], p.z = ((double *)&p1)[1], p.x = 0.;
  return p;
}
#line 96
static inline coord embed_face_barycentre_y (Point point, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_z = Delta;

  double Delta_x = Delta;


  double Delta_y = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_z);

  NOT_UNUSED(Delta_x);


  NOT_UNUSED(Delta_y);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 99 "/home/pwachara/basilisk/src/embed.h"
coord n1 = {0};
  double nn = 0.;
  scalar f = fs.y;
   {
    n1.z = (val(f,-1,i,-1) + 2.*val(f,0,i,-1) + val(f,1,i,-1) -
     val(f,-1,i,+1) - 2.*val(f,0,i,+1) - val(f,1,i,+1));
    nn += fabs(n1.z);
  } 
#line 102
{
    n1.x = (val(f,-1,i,-1) + 2.*val(f,-1,i,0) + val(f,-1,i,1) -
     val(f,+1,i,-1) - 2.*val(f,+1,i,0) - val(f,+1,i,1));
    nn += fabs(n1.x);
  }
  if (!nn)
    return (coord){0.,0.,0.};
  
    n1.z /= nn;    n1.x /= nn;

  coord n, p1, p;
  ((double *)&n)[0] = n1.z, ((double *)&n)[1] = n1.x;
  double alpha = line_alpha (val(f,0,i,0), n);
  line_center (n, alpha, val(f,0,i,0), &p1);
  p.z = ((double *)&p1)[0], p.x = ((double *)&p1)[1], p.y = 0.;
  return p;
}
#line 95 "/home/pwachara/basilisk/src/embed.h"

static void _stencil_embed_face_barycentre_z (Point point, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 101 "/home/pwachara/basilisk/src/embed.h"
scalar f = fs.z;
   {
_stencil_val(f,-1,-1,i);_stencil_val(f,-1,0,i); _stencil_val(f,-1,1,i);
     _stencil_val(f,+1,-1,i);_stencil_val(f,+1,0,i); _stencil_val(f,+1,1,i);  

} 
#line 102
{
_stencil_val(f,-1,-1,i);_stencil_val(f,0,-1,i); _stencil_val(f,1,-1,i);
     _stencil_val(f,-1,+1,i);_stencil_val(f,0,+1,i); _stencil_val(f,1,+1,i);  

}         
    
   
  
      

    
       
  _stencil_val(f,0,0,i); 
_stencil_val(f,0,0,i);  

return ;
}
#line 96
static void _stencil_embed_face_barycentre_x (Point point, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_y = Delta;

  double Delta_z = Delta;


  double Delta_x = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_y);

  NOT_UNUSED(Delta_z);


  NOT_UNUSED(Delta_x);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 101 "/home/pwachara/basilisk/src/embed.h"
scalar f = fs.x;
   {
_stencil_val(f,i,-1,-1);_stencil_val(f,i,-1,0); _stencil_val(f,i,-1,1);
     _stencil_val(f,i,+1,-1);_stencil_val(f,i,+1,0); _stencil_val(f,i,+1,1);  

} 
#line 102
{
_stencil_val(f,i,-1,-1);_stencil_val(f,i,0,-1); _stencil_val(f,i,1,-1);
     _stencil_val(f,i,-1,+1);_stencil_val(f,i,0,+1); _stencil_val(f,i,1,+1);  

}         
    
   
  
      

    
       
  _stencil_val(f,i,0,0); 
_stencil_val(f,i,0,0);  

return ;
}
#line 96
static void _stencil_embed_face_barycentre_y (Point point, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_z = Delta;

  double Delta_x = Delta;


  double Delta_y = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_z);

  NOT_UNUSED(Delta_x);


  NOT_UNUSED(Delta_y);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 101 "/home/pwachara/basilisk/src/embed.h"
scalar f = fs.y;
   {
_stencil_val(f,-1,i,-1);_stencil_val(f,0,i,-1); _stencil_val(f,1,i,-1);
     _stencil_val(f,-1,i,+1);_stencil_val(f,0,i,+1); _stencil_val(f,1,i,+1);  

} 
#line 102
{
_stencil_val(f,-1,i,-1);_stencil_val(f,-1,i,0); _stencil_val(f,-1,i,1);
     _stencil_val(f,+1,i,-1);_stencil_val(f,+1,i,0); _stencil_val(f,+1,i,1);  

}         
    
   
  
      

    
       
  _stencil_val(f,0,i,0); 
_stencil_val(f,0,i,0);  

return ;
}
#line 138 "/home/pwachara/basilisk/src/embed.h"

static inline double embed_face_gradient_x (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 141 "/home/pwachara/basilisk/src/embed.h"
if (!(val(cs,i,0,0) && val(cs,i-1,0,0))) qassert ("/home/pwachara/basilisk/src/embed.h", 141, "cs[i] && cs[i-1]");
  coord p = embed_face_barycentre_x (point, i);

  int j = ( (int)((p.y) > 0 ? 1 : -1)), k = ( (int)((p.z) > 0 ? 1 : -1));
  if ((val(fs.x,i,j,k) > 0.5 && (val(fs.x,i,j,0) > 0.5 || val(fs.x,i,0,k) > 0.5) && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(fs.y,i,j + (j < 0),k) && val(fs.y,i-1,j + (j < 0),k) && val(fs.z,i,0,k + (k < 0)) && val(fs.z,i-1,0,k + (k < 0)) && val(fs.z,i,j,k + (k < 0)) && val(fs.z,i-1,j,k + (k < 0)) && val(cs,i-1,j,0) && val(cs,i-1,0,k) && val(cs,i-1,j,k) && val(cs,i,j,0) && val(cs,i,0,k) && val(cs,i,j,k))) {
    p.y = fabs(p.y), p.z = fabs(p.z);
    return (((val(a,i,0,0) - val(a,i-1,0,0))*(1. - p.y) +
      (val(a,i,j,0) - val(a,i-1,j,0))*p.y)*(1. - p.z) +
     ((val(a,i,0,k) - val(a,i-1,0,k))*(1. - p.y) +
      (val(a,i,j,k) - val(a,i-1,j,k))*p.y)*p.z)/Delta;
  }
  return (val(a,i,0,0) - val(a,i-1,0,0))/Delta;
}
#line 139
static inline double embed_face_gradient_y (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_y = Delta;

  double Delta_z = Delta;


  double Delta_x = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_y);

  NOT_UNUSED(Delta_z);


  NOT_UNUSED(Delta_x);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 141 "/home/pwachara/basilisk/src/embed.h"
if (!(val(cs,0,i,0) && val(cs,0,i-1,0))) qassert ("/home/pwachara/basilisk/src/embed.h", 141, "cs[i] && cs[i-1]");
  coord p = embed_face_barycentre_y (point, i);

  int j = ( (int)((p.z) > 0 ? 1 : -1)), k = ( (int)((p.x) > 0 ? 1 : -1));
  if ((val(fs.y,k,i,j) > 0.5 && (val(fs.y,0,i,j) > 0.5 || val(fs.y,k,i,0) > 0.5) && val(fs.z,0,i,j + (j < 0)) && val(fs.z,0,i-1,j + (j < 0)) && val(fs.z,k,i,j + (j < 0)) && val(fs.z,k,i-1,j + (j < 0)) && val(fs.x,k + (k < 0),i,0) && val(fs.x,k + (k < 0),i-1,0) && val(fs.x,k + (k < 0),i,j) && val(fs.x,k + (k < 0),i-1,j) && val(cs,0,i-1,j) && val(cs,k,i-1,0) && val(cs,k,i-1,j) && val(cs,0,i,j) && val(cs,k,i,0) && val(cs,k,i,j))) {
    p.z = fabs(p.z), p.x = fabs(p.x);
    return (((val(a,0,i,0) - val(a,0,i-1,0))*(1. - p.z) +
      (val(a,0,i,j) - val(a,0,i-1,j))*p.z)*(1. - p.x) +
     ((val(a,k,i,0) - val(a,k,i-1,0))*(1. - p.z) +
      (val(a,k,i,j) - val(a,k,i-1,j))*p.z)*p.x)/Delta;
  }
  return (val(a,0,i,0) - val(a,0,i-1,0))/Delta;
}
#line 139
static inline double embed_face_gradient_z (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_z = Delta;

  double Delta_x = Delta;


  double Delta_y = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_z);

  NOT_UNUSED(Delta_x);


  NOT_UNUSED(Delta_y);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 141 "/home/pwachara/basilisk/src/embed.h"
if (!(val(cs,0,0,i) && val(cs,0,0,i-1))) qassert ("/home/pwachara/basilisk/src/embed.h", 141, "cs[i] && cs[i-1]");
  coord p = embed_face_barycentre_z (point, i);

  int j = ( (int)((p.x) > 0 ? 1 : -1)), k = ( (int)((p.y) > 0 ? 1 : -1));
  if ((val(fs.z,j,k,i) > 0.5 && (val(fs.z,j,0,i) > 0.5 || val(fs.z,0,k,i) > 0.5) && val(fs.x,j + (j < 0),0,i) && val(fs.x,j + (j < 0),0,i-1) && val(fs.x,j + (j < 0),k,i) && val(fs.x,j + (j < 0),k,i-1) && val(fs.y,0,k + (k < 0),i) && val(fs.y,0,k + (k < 0),i-1) && val(fs.y,j,k + (k < 0),i) && val(fs.y,j,k + (k < 0),i-1) && val(cs,j,0,i-1) && val(cs,0,k,i-1) && val(cs,j,k,i-1) && val(cs,j,0,i) && val(cs,0,k,i) && val(cs,j,k,i))) {
    p.x = fabs(p.x), p.y = fabs(p.y);
    return (((val(a,0,0,i) - val(a,0,0,i-1))*(1. - p.x) +
      (val(a,j,0,i) - val(a,j,0,i-1))*p.x)*(1. - p.y) +
     ((val(a,0,k,i) - val(a,0,k,i-1))*(1. - p.x) +
      (val(a,j,k,i) - val(a,j,k,i-1))*p.x)*p.y)/Delta;
  }
  return (val(a,0,0,i) - val(a,0,0,i-1))/Delta;
}
#line 138 "/home/pwachara/basilisk/src/embed.h"

static void _stencil_embed_face_gradient_x (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 141 "/home/pwachara/basilisk/src/embed.h"
_stencil_val(cs,i,0,0); _stencil_val(cs,i-1,0,0);   
_stencil_embed_face_barycentre_x (point, i);                    

  
_stencil_val(fs.x,i,o_stencil,o_stencil);_stencil_val(fs.x,i,o_stencil,0); _stencil_val(fs.x,i,0,o_stencil); _stencil_val(fs.y,i,o_stencil    ,0); _stencil_val(fs.y,i-1,o_stencil    ,0); _stencil_val(fs.y,i,o_stencil,    o_stencil); _stencil_val(fs.y,i-1,o_stencil,    o_stencil); _stencil_val(fs.z,i,0,o_stencil    ); _stencil_val(fs.z,i-1,0,o_stencil    ); _stencil_val(fs.z,i,o_stencil,o_stencil    ); _stencil_val(fs.z,i-1,o_stencil,o_stencil    ); _stencil_val(cs,i-1,o_stencil,0); _stencil_val(cs,i-1,0,o_stencil); _stencil_val(cs,i-1,o_stencil,o_stencil); _stencil_val(cs,i,o_stencil,0); _stencil_val(cs,i,0,o_stencil); _stencil_val(cs,i,o_stencil,o_stencil); {
         
_stencil_val(a,i,0,0); _stencil_val(a,i-1,0,0);
_stencil_val(a,i,o_stencil,0); _stencil_val(a,i-1,o_stencil,0);
_stencil_val(a,i,0,o_stencil); _stencil_val(a,i-1,0,o_stencil);
_stencil_val(a,i,o_stencil,o_stencil); _stencil_val(a,i-1,o_stencil,o_stencil);  
}
_stencil_val(a,i,0,0); _stencil_val(a,i-1,0,0);  return  ;
}
#line 139
static void _stencil_embed_face_gradient_y (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_y = Delta;

  double Delta_z = Delta;


  double Delta_x = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_y);

  NOT_UNUSED(Delta_z);


  NOT_UNUSED(Delta_x);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 141 "/home/pwachara/basilisk/src/embed.h"
_stencil_val(cs,0,i,0); _stencil_val(cs,0,i-1,0);   
_stencil_embed_face_barycentre_y (point, i);                    

  
_stencil_val(fs.y,o_stencil,i,o_stencil);_stencil_val(fs.y,0,i,o_stencil); _stencil_val(fs.y,o_stencil,i,0); _stencil_val(fs.z,0,i    ,o_stencil); _stencil_val(fs.z,0,i-1    ,o_stencil); _stencil_val(fs.z,    o_stencil,i,o_stencil); _stencil_val(fs.z,    o_stencil,i-1,o_stencil); _stencil_val(fs.x,o_stencil,i,0    ); _stencil_val(fs.x,o_stencil,i-1,0    ); _stencil_val(fs.x,o_stencil,i,o_stencil    ); _stencil_val(fs.x,o_stencil,i-1,o_stencil    ); _stencil_val(cs,0,i-1,o_stencil); _stencil_val(cs,o_stencil,i-1,0); _stencil_val(cs,o_stencil,i-1,o_stencil); _stencil_val(cs,0,i,o_stencil); _stencil_val(cs,o_stencil,i,0); _stencil_val(cs,o_stencil,i,o_stencil); {
         
_stencil_val(a,0,i,0); _stencil_val(a,0,i-1,0);
_stencil_val(a,0,i,o_stencil); _stencil_val(a,0,i-1,o_stencil);
_stencil_val(a,o_stencil,i,0); _stencil_val(a,o_stencil,i-1,0);
_stencil_val(a,o_stencil,i,o_stencil); _stencil_val(a,o_stencil,i-1,o_stencil);  
}
_stencil_val(a,0,i,0); _stencil_val(a,0,i-1,0);  return  ;
}
#line 139
static void _stencil_embed_face_gradient_z (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_z = Delta;

  double Delta_x = Delta;


  double Delta_y = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_z);

  NOT_UNUSED(Delta_x);


  NOT_UNUSED(Delta_y);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 141 "/home/pwachara/basilisk/src/embed.h"
_stencil_val(cs,0,0,i); _stencil_val(cs,0,0,i-1);   
_stencil_embed_face_barycentre_z (point, i);                    

  
_stencil_val(fs.z,o_stencil,o_stencil,i);_stencil_val(fs.z,o_stencil,0,i); _stencil_val(fs.z,0,o_stencil,i); _stencil_val(fs.x,o_stencil,0    ,i); _stencil_val(fs.x,o_stencil,0    ,i-1); _stencil_val(fs.x,o_stencil,    o_stencil,i); _stencil_val(fs.x,o_stencil,    o_stencil,i-1); _stencil_val(fs.y,0,o_stencil,i    ); _stencil_val(fs.y,0,o_stencil,i-1    ); _stencil_val(fs.y,o_stencil,o_stencil,i    ); _stencil_val(fs.y,o_stencil,o_stencil,i-1    ); _stencil_val(cs,o_stencil,0,i-1); _stencil_val(cs,0,o_stencil,i-1); _stencil_val(cs,o_stencil,o_stencil,i-1); _stencil_val(cs,o_stencil,0,i); _stencil_val(cs,0,o_stencil,i); _stencil_val(cs,o_stencil,o_stencil,i); {
         
_stencil_val(a,0,0,i); _stencil_val(a,0,0,i-1);
_stencil_val(a,o_stencil,0,i); _stencil_val(a,o_stencil,0,i-1);
_stencil_val(a,0,o_stencil,i); _stencil_val(a,0,o_stencil,i-1);
_stencil_val(a,o_stencil,o_stencil,i); _stencil_val(a,o_stencil,o_stencil,i-1);  
}
_stencil_val(a,0,0,i); _stencil_val(a,0,0,i-1);  return  ;
}


static inline double embed_face_value_x (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 158 "/home/pwachara/basilisk/src/embed.h"
coord p = embed_face_barycentre_x (point, i);

  int j = ( (int)((p.y) > 0 ? 1 : -1)), k = ( (int)((p.z) > 0 ? 1 : -1));
  if ((val(fs.x,i,j,k) > 0.5 && (val(fs.x,i,j,0) > 0.5 || val(fs.x,i,0,k) > 0.5) && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(fs.y,i,j + (j < 0),k) && val(fs.y,i-1,j + (j < 0),k) && val(fs.z,i,0,k + (k < 0)) && val(fs.z,i-1,0,k + (k < 0)) && val(fs.z,i,j,k + (k < 0)) && val(fs.z,i-1,j,k + (k < 0)) && val(cs,i-1,j,0) && val(cs,i-1,0,k) && val(cs,i-1,j,k) && val(cs,i,j,0) && val(cs,i,0,k) && val(cs,i,j,k))) {
    p.y = fabs(p.y), p.z = fabs(p.z);
    return ((((val(a,i,0,0)*(1.5 + val(cs,i,0,0)) + val(a,i-1,0,0)*(1.5 + val(cs,i-1,0,0)))/ (val(cs,i,0,0) + val(cs,i-1,0,0) + 3.))*(1. - p.y) + ((val(a,i,j,0)*(1.5 + val(cs,i,j,0)) + val(a,i-1,j,0)*(1.5 + val(cs,i-1,j,0)))/ (val(cs,i,j,0) + val(cs,i-1,j,0) + 3.))*p.y)*(1. - p.z) +
     (((val(a,i,0,k)*(1.5 + val(cs,i,0,k)) + val(a,i-1,0,k)*(1.5 + val(cs,i-1,0,k)))/ (val(cs,i,0,k) + val(cs,i-1,0,k) + 3.))*(1. - p.y) + ((val(a,i,j,k)*(1.5 + val(cs,i,j,k)) + val(a,i-1,j,k)*(1.5 + val(cs,i-1,j,k)))/ (val(cs,i,j,k) + val(cs,i-1,j,k) + 3.))*p.y)*p.z);
  }
  return ((val(a,i,0,0)*(1.5 + val(cs,i,0,0)) + val(a,i-1,0,0)*(1.5 + val(cs,i-1,0,0)))/ (val(cs,i,0,0) + val(cs,i-1,0,0) + 3.));
}
#line 156
static inline double embed_face_value_y (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_y = Delta;

  double Delta_z = Delta;


  double Delta_x = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_y);

  NOT_UNUSED(Delta_z);


  NOT_UNUSED(Delta_x);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 158 "/home/pwachara/basilisk/src/embed.h"
coord p = embed_face_barycentre_y (point, i);

  int j = ( (int)((p.z) > 0 ? 1 : -1)), k = ( (int)((p.x) > 0 ? 1 : -1));
  if ((val(fs.y,k,i,j) > 0.5 && (val(fs.y,0,i,j) > 0.5 || val(fs.y,k,i,0) > 0.5) && val(fs.z,0,i,j + (j < 0)) && val(fs.z,0,i-1,j + (j < 0)) && val(fs.z,k,i,j + (j < 0)) && val(fs.z,k,i-1,j + (j < 0)) && val(fs.x,k + (k < 0),i,0) && val(fs.x,k + (k < 0),i-1,0) && val(fs.x,k + (k < 0),i,j) && val(fs.x,k + (k < 0),i-1,j) && val(cs,0,i-1,j) && val(cs,k,i-1,0) && val(cs,k,i-1,j) && val(cs,0,i,j) && val(cs,k,i,0) && val(cs,k,i,j))) {
    p.z = fabs(p.z), p.x = fabs(p.x);
    return ((((val(a,0,i,0)*(1.5 + val(cs,0,i,0)) + val(a,0,i-1,0)*(1.5 + val(cs,0,i-1,0)))/ (val(cs,0,i,0) + val(cs,0,i-1,0) + 3.))*(1. - p.z) + ((val(a,0,i,j)*(1.5 + val(cs,0,i,j)) + val(a,0,i-1,j)*(1.5 + val(cs,0,i-1,j)))/ (val(cs,0,i,j) + val(cs,0,i-1,j) + 3.))*p.z)*(1. - p.x) +
     (((val(a,k,i,0)*(1.5 + val(cs,k,i,0)) + val(a,k,i-1,0)*(1.5 + val(cs,k,i-1,0)))/ (val(cs,k,i,0) + val(cs,k,i-1,0) + 3.))*(1. - p.z) + ((val(a,k,i,j)*(1.5 + val(cs,k,i,j)) + val(a,k,i-1,j)*(1.5 + val(cs,k,i-1,j)))/ (val(cs,k,i,j) + val(cs,k,i-1,j) + 3.))*p.z)*p.x);
  }
  return ((val(a,0,i,0)*(1.5 + val(cs,0,i,0)) + val(a,0,i-1,0)*(1.5 + val(cs,0,i-1,0)))/ (val(cs,0,i,0) + val(cs,0,i-1,0) + 3.));
}
#line 156
static inline double embed_face_value_z (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_z = Delta;

  double Delta_x = Delta;


  double Delta_y = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_z);

  NOT_UNUSED(Delta_x);


  NOT_UNUSED(Delta_y);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;  
#line 158 "/home/pwachara/basilisk/src/embed.h"
coord p = embed_face_barycentre_z (point, i);

  int j = ( (int)((p.x) > 0 ? 1 : -1)), k = ( (int)((p.y) > 0 ? 1 : -1));
  if ((val(fs.z,j,k,i) > 0.5 && (val(fs.z,j,0,i) > 0.5 || val(fs.z,0,k,i) > 0.5) && val(fs.x,j + (j < 0),0,i) && val(fs.x,j + (j < 0),0,i-1) && val(fs.x,j + (j < 0),k,i) && val(fs.x,j + (j < 0),k,i-1) && val(fs.y,0,k + (k < 0),i) && val(fs.y,0,k + (k < 0),i-1) && val(fs.y,j,k + (k < 0),i) && val(fs.y,j,k + (k < 0),i-1) && val(cs,j,0,i-1) && val(cs,0,k,i-1) && val(cs,j,k,i-1) && val(cs,j,0,i) && val(cs,0,k,i) && val(cs,j,k,i))) {
    p.x = fabs(p.x), p.y = fabs(p.y);
    return ((((val(a,0,0,i)*(1.5 + val(cs,0,0,i)) + val(a,0,0,i-1)*(1.5 + val(cs,0,0,i-1)))/ (val(cs,0,0,i) + val(cs,0,0,i-1) + 3.))*(1. - p.x) + ((val(a,j,0,i)*(1.5 + val(cs,j,0,i)) + val(a,j,0,i-1)*(1.5 + val(cs,j,0,i-1)))/ (val(cs,j,0,i) + val(cs,j,0,i-1) + 3.))*p.x)*(1. - p.y) +
     (((val(a,0,k,i)*(1.5 + val(cs,0,k,i)) + val(a,0,k,i-1)*(1.5 + val(cs,0,k,i-1)))/ (val(cs,0,k,i) + val(cs,0,k,i-1) + 3.))*(1. - p.x) + ((val(a,j,k,i)*(1.5 + val(cs,j,k,i)) + val(a,j,k,i-1)*(1.5 + val(cs,j,k,i-1)))/ (val(cs,j,k,i) + val(cs,j,k,i-1) + 3.))*p.x)*p.y);
  }
  return ((val(a,0,0,i)*(1.5 + val(cs,0,0,i)) + val(a,0,0,i-1)*(1.5 + val(cs,0,0,i-1)))/ (val(cs,0,0,i) + val(cs,0,0,i-1) + 3.));
}
#line 156
static void _stencil_embed_face_value_x (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;   
#line 158 "/home/pwachara/basilisk/src/embed.h"
_stencil_embed_face_barycentre_x (point, i);                    

  
_stencil_val(fs.x,i,o_stencil,o_stencil);_stencil_val(fs.x,i,o_stencil,0); _stencil_val(fs.x,i,0,o_stencil); _stencil_val(fs.y,i,o_stencil    ,0); _stencil_val(fs.y,i-1,o_stencil    ,0); _stencil_val(fs.y,i,o_stencil,    o_stencil); _stencil_val(fs.y,i-1,o_stencil,    o_stencil); _stencil_val(fs.z,i,0,o_stencil    ); _stencil_val(fs.z,i-1,0,o_stencil    ); _stencil_val(fs.z,i,o_stencil,o_stencil    ); _stencil_val(fs.z,i-1,o_stencil,o_stencil    ); _stencil_val(cs,i-1,o_stencil,0); _stencil_val(cs,i-1,0,o_stencil); _stencil_val(cs,i-1,o_stencil,o_stencil); _stencil_val(cs,i,o_stencil,0); _stencil_val(cs,i,0,o_stencil); _stencil_val(cs,i,o_stencil,o_stencil); {
         
_stencil_val(a,i,0,0); _stencil_val(cs,i,0,0); _stencil_val(a,i-1,0,0); _stencil_val(cs,i-1,0,0);_stencil_val(cs,i,0,0); _stencil_val(cs,i-1,0,0);_stencil_val(a,i,o_stencil,0); _stencil_val(cs,i,o_stencil,0); _stencil_val(a,i-1,o_stencil,0); _stencil_val(cs,i-1,o_stencil,0);_stencil_val(cs,i,o_stencil,0); _stencil_val(cs,i-1,o_stencil,0);
_stencil_val(a,i,0,o_stencil); _stencil_val(cs,i,0,o_stencil); _stencil_val(a,i-1,0,o_stencil); _stencil_val(cs,i-1,0,o_stencil);_stencil_val(cs,i,0,o_stencil); _stencil_val(cs,i-1,0,o_stencil);_stencil_val(a,i,o_stencil,o_stencil); _stencil_val(cs,i,o_stencil,o_stencil); _stencil_val(a,i-1,o_stencil,o_stencil); _stencil_val(cs,i-1,o_stencil,o_stencil);_stencil_val(cs,i,o_stencil,o_stencil); _stencil_val(cs,i-1,o_stencil,o_stencil);  
}
_stencil_val(a,i,0,0); _stencil_val(cs,i,0,0); _stencil_val(a,i-1,0,0); _stencil_val(cs,i-1,0,0);_stencil_val(cs,i,0,0); _stencil_val(cs,i-1,0,0);  return        ;
}
#line 156
static void _stencil_embed_face_value_y (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_y = Delta;

  double Delta_z = Delta;


  double Delta_x = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_y);

  NOT_UNUSED(Delta_z);


  NOT_UNUSED(Delta_x);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;   
#line 158 "/home/pwachara/basilisk/src/embed.h"
_stencil_embed_face_barycentre_y (point, i);                    

  
_stencil_val(fs.y,o_stencil,i,o_stencil);_stencil_val(fs.y,0,i,o_stencil); _stencil_val(fs.y,o_stencil,i,0); _stencil_val(fs.z,0,i    ,o_stencil); _stencil_val(fs.z,0,i-1    ,o_stencil); _stencil_val(fs.z,    o_stencil,i,o_stencil); _stencil_val(fs.z,    o_stencil,i-1,o_stencil); _stencil_val(fs.x,o_stencil,i,0    ); _stencil_val(fs.x,o_stencil,i-1,0    ); _stencil_val(fs.x,o_stencil,i,o_stencil    ); _stencil_val(fs.x,o_stencil,i-1,o_stencil    ); _stencil_val(cs,0,i-1,o_stencil); _stencil_val(cs,o_stencil,i-1,0); _stencil_val(cs,o_stencil,i-1,o_stencil); _stencil_val(cs,0,i,o_stencil); _stencil_val(cs,o_stencil,i,0); _stencil_val(cs,o_stencil,i,o_stencil); {
         
_stencil_val(a,0,i,0); _stencil_val(cs,0,i,0); _stencil_val(a,0,i-1,0); _stencil_val(cs,0,i-1,0);_stencil_val(cs,0,i,0); _stencil_val(cs,0,i-1,0);_stencil_val(a,0,i,o_stencil); _stencil_val(cs,0,i,o_stencil); _stencil_val(a,0,i-1,o_stencil); _stencil_val(cs,0,i-1,o_stencil);_stencil_val(cs,0,i,o_stencil); _stencil_val(cs,0,i-1,o_stencil);
_stencil_val(a,o_stencil,i,0); _stencil_val(cs,o_stencil,i,0); _stencil_val(a,o_stencil,i-1,0); _stencil_val(cs,o_stencil,i-1,0);_stencil_val(cs,o_stencil,i,0); _stencil_val(cs,o_stencil,i-1,0);_stencil_val(a,o_stencil,i,o_stencil); _stencil_val(cs,o_stencil,i,o_stencil); _stencil_val(a,o_stencil,i-1,o_stencil); _stencil_val(cs,o_stencil,i-1,o_stencil);_stencil_val(cs,o_stencil,i,o_stencil); _stencil_val(cs,o_stencil,i-1,o_stencil);  
}
_stencil_val(a,0,i,0); _stencil_val(cs,0,i,0); _stencil_val(a,0,i-1,0); _stencil_val(cs,0,i-1,0);_stencil_val(cs,0,i,0); _stencil_val(cs,0,i-1,0);  return        ;
}
#line 156
static void _stencil_embed_face_value_z (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_z = Delta;

  double Delta_x = Delta;


  double Delta_y = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_z);

  NOT_UNUSED(Delta_x);


  NOT_UNUSED(Delta_y);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;   
#line 158 "/home/pwachara/basilisk/src/embed.h"
_stencil_embed_face_barycentre_z (point, i);                    

  
_stencil_val(fs.z,o_stencil,o_stencil,i);_stencil_val(fs.z,o_stencil,0,i); _stencil_val(fs.z,0,o_stencil,i); _stencil_val(fs.x,o_stencil,0    ,i); _stencil_val(fs.x,o_stencil,0    ,i-1); _stencil_val(fs.x,o_stencil,    o_stencil,i); _stencil_val(fs.x,o_stencil,    o_stencil,i-1); _stencil_val(fs.y,0,o_stencil,i    ); _stencil_val(fs.y,0,o_stencil,i-1    ); _stencil_val(fs.y,o_stencil,o_stencil,i    ); _stencil_val(fs.y,o_stencil,o_stencil,i-1    ); _stencil_val(cs,o_stencil,0,i-1); _stencil_val(cs,0,o_stencil,i-1); _stencil_val(cs,o_stencil,o_stencil,i-1); _stencil_val(cs,o_stencil,0,i); _stencil_val(cs,0,o_stencil,i); _stencil_val(cs,o_stencil,o_stencil,i); {
         
_stencil_val(a,0,0,i); _stencil_val(cs,0,0,i); _stencil_val(a,0,0,i-1); _stencil_val(cs,0,0,i-1);_stencil_val(cs,0,0,i); _stencil_val(cs,0,0,i-1);_stencil_val(a,o_stencil,0,i); _stencil_val(cs,o_stencil,0,i); _stencil_val(a,o_stencil,0,i-1); _stencil_val(cs,o_stencil,0,i-1);_stencil_val(cs,o_stencil,0,i); _stencil_val(cs,o_stencil,0,i-1);
_stencil_val(a,0,o_stencil,i); _stencil_val(cs,0,o_stencil,i); _stencil_val(a,0,o_stencil,i-1); _stencil_val(cs,0,o_stencil,i-1);_stencil_val(cs,0,o_stencil,i); _stencil_val(cs,0,o_stencil,i-1);_stencil_val(a,o_stencil,o_stencil,i); _stencil_val(cs,o_stencil,o_stencil,i); _stencil_val(a,o_stencil,o_stencil,i-1); _stencil_val(cs,o_stencil,o_stencil,i-1);_stencil_val(cs,o_stencil,o_stencil,i); _stencil_val(cs,o_stencil,o_stencil,i-1);  
}
_stencil_val(a,0,0,i); _stencil_val(cs,0,0,i); _stencil_val(a,0,0,i-1); _stencil_val(cs,0,0,i-1);_stencil_val(cs,0,0,i); _stencil_val(cs,0,0,i-1);  return        ;
}
#line 177 "/home/pwachara/basilisk/src/embed.h"

#line 222 "/home/pwachara/basilisk/src/embed.h"
static inline
double embed_geometry (Point point, coord * p, coord * n)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
  
#line 225 "/home/pwachara/basilisk/src/embed.h"
*n = facet_normal (point, cs, fs);
  double alpha = plane_alpha (val(cs,0,0,0), *n);
  double area = plane_area_center (*n, alpha, p);
  normalize (n);
  return area;
}
#line 222 "/home/pwachara/basilisk/src/embed.h"
static void 
_stencil_embed_geometry (Point point,_stencil_undefined * p,_stencil_undefined * n)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2; 
#line 225 "/home/pwachara/basilisk/src/embed.h"
_stencil_facet_normal (point, cs, fs);  
_stencil_val(cs,0,0,0);      
   
  
  return ;
}





static inline
double embed_area_center (Point point, double * x1, double * y1, double * z1)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
  
#line 239 "/home/pwachara/basilisk/src/embed.h"
double area = 0.;
  if (val(cs,0,0,0) > 0. && val(cs,0,0,0) < 1.) {
    coord n, p;
    area = embed_geometry (point, &p, &n);
    *x1 += p.x*Delta, *y1 += p.y*Delta, *z1 += p.z*Delta;
  }
  return area;
}
#line 255 "/home/pwachara/basilisk/src/embed.h"
double embed_interpolate (Point point, scalar s, coord p)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
  
#line 257 "/home/pwachara/basilisk/src/embed.h"
if (!(3 == 2)) qassert ("/home/pwachara/basilisk/src/embed.h", 257, "dimension == 2");
  int i = ( (int)((p.x) > 0 ? 1 : -1)), j = ( (int)((p.y) > 0 ? 1 : -1));
  if (val(cs,i,0,0) && val(cs,0,j,0) && val(cs,i,j,0))

    return ((val(s,0,0,0)*(1. - fabs(p.x)) + val(s,i,0,0)*fabs(p.x))*(1. - fabs(p.y)) +
     (val(s,0,j,0)*(1. - fabs(p.x)) + val(s,i,j,0)*fabs(p.x))*fabs(p.y));
  else {


    double val = val(s,0,0,0);
     {
      int i = ( (int)((p.x) > 0 ? 1 : -1));
      if (val(cs,i,0,0))
 val += fabs(p.x)*(val(s,i,0,0) - val(s,0,0,0));
      else if (val(cs,-i,0,0))
 val += fabs(p.x)*(val(s,0,0,0) - val(s,-i,0,0));
    } 
#line 267
{
      int i = ( (int)((p.y) > 0 ? 1 : -1));
      if (val(cs,0,i,0))
 val += fabs(p.y)*(val(s,0,i,0) - val(s,0,0,0));
      else if (val(cs,0,-i,0))
 val += fabs(p.y)*(val(s,0,0,0) - val(s,0,-i,0));
    } 
#line 267
{
      int i = ( (int)((p.z) > 0 ? 1 : -1));
      if (val(cs,0,0,i))
 val += fabs(p.z)*(val(s,0,0,i) - val(s,0,0,0));
      else if (val(cs,0,0,-i))
 val += fabs(p.z)*(val(s,0,0,0) - val(s,0,0,-i));
    }
    return val;
  }
}
#line 255 "/home/pwachara/basilisk/src/embed.h"
static void _stencil_embed_interpolate (Point point, scalar s,_stencil_undefined * p)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;                       
      
  
#line 259 "/home/pwachara/basilisk/src/embed.h"
_stencil_val(cs,o_stencil,0,0); _stencil_val(cs,0,o_stencil,0); _stencil_val(cs,o_stencil,o_stencil,0);{

    {_stencil_val(s,0,0,0);_stencil_val(s, o_stencil,0,0);
_stencil_val(s,0,o_stencil,0); _stencil_val(s,o_stencil,o_stencil,0);      } 
{  


     _stencil_val(s,0,0,0);
     {          
      
_stencil_val(cs,o_stencil,0,0);{
 {_stencil_val(s,o_stencil,0,0); _stencil_val(s,0,0,0);   } 
{_stencil_val(cs,o_stencil,0,0);
 {_stencil_val(s,0,0,0);_stencil_val(s, o_stencil,0,0);   } }}    
} 
#line 267
{          
      
_stencil_val(cs,0,o_stencil,0);{
 {_stencil_val(s,0,o_stencil,0); _stencil_val(s,0,0,0);   } 
{_stencil_val(cs,0,o_stencil,0);
 {_stencil_val(s,0,0,0);_stencil_val(s,0, o_stencil,0);   } }}    
} 
#line 267
{          
      
_stencil_val(cs,0,0,o_stencil);{
 {_stencil_val(s,0,0,o_stencil); _stencil_val(s,0,0,0);   } 
{_stencil_val(cs,0,0,o_stencil);
 {_stencil_val(s,0,0,0);_stencil_val(s,0,0, o_stencil);   } }}    
} 
    
  }}
}
#line 285 "/home/pwachara/basilisk/src/embed.h"
struct Cleanup {
  scalar c;
  vector s;
  double smin;
  bool opposite;
};

     
int fractions_cleanup (scalar c, vector s,
         double smin, bool opposite)
{tracing("fractions_cleanup","/home/pwachara/basilisk/src/embed.h",293);







  int changed = 1, schanged = 0, i;
  for (i = 0; i < 100 && changed; i++) {
  
  
#line 314 "/home/pwachara/basilisk/src/grid/stencils.h"
{
    static int _first = 1.;
    ForeachData _loop = {
      .fname = "/home/pwachara/basilisk/src/embed.h", .line = 309, .first = _first
    };
    if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
 _attribute[s.i].input = _attribute[s.i].output = _attribute[s.i].nowarning = false;
 _attribute[s.i].width = 0;
      }
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point);
#line 309 "/home/pwachara/basilisk/src/embed.h"
{ 
#line 377 "/home/pwachara/basilisk/src/grid/stencils.h"
_loop.face |= (1 << 0);      
#line 310 "/home/pwachara/basilisk/src/embed.h"
{_stencil_val(s.x,0,0,0);_stencil_val(c,0,0,0);_stencil_val(c,-1,0,0); _stencil_val(s.x,0,0,0);
 {_stencil_val_a(s.x,0,0,0);  }        } 
#line 378 "/home/pwachara/basilisk/src/grid/stencils.h"
_loop.face |= (1 << 1);      
#line 310 "/home/pwachara/basilisk/src/embed.h"
{_stencil_val(s.y,0,0,0);_stencil_val(c,0,0,0);_stencil_val(c,0,-1,0); _stencil_val(s.y,0,0,0);
 {_stencil_val_a(s.y,0,0,0);  }        } 
#line 379 "/home/pwachara/basilisk/src/grid/stencils.h"
_loop.face |= (1 << 2);      
#line 310 "/home/pwachara/basilisk/src/embed.h"
{_stencil_val(s.z,0,0,0);_stencil_val(c,0,0,0);_stencil_val(c,0,0,-1); _stencil_val(s.z,0,0,0);
 {_stencil_val_a(s.z,0,0,0);  }        }}

    
#line 328 "/home/pwachara/basilisk/src/grid/stencils.h"
check_stencil (&_loop);
    boundary_stencil (&_loop);
    _first = 0;
  }
#line 300 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
  OMP_PARALLEL () {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
      for (_k = 2; _k <= point.n.x + 2; _k++) {
 point.i = _k;

 for (point.j = 2; point.j <= point.n.y + 2; point.j++)

   for (point.k = 2; point.k <= point.n.z + 2; point.k++)
#line 309 "/home/pwachara/basilisk/src/embed.h"
{
  
#line 386 "/home/pwachara/basilisk/src/grid/multigrid.h"
if (point.j < point.n.y + 2 && point.k < point.n.z + 2) {
    int ig = -1; NOT_UNUSED(ig);      
#line 310 "/home/pwachara/basilisk/src/embed.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 310 "/home/pwachara/basilisk/src/embed.h"
if (val(s.x,0,0,0) && ((!val(c,0,0,0) || !val(c,-1,0,0)) || val(s.x,0,0,0) < smin))
 val(s.x,0,0,0) = 0.;}
  
#line 389 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
  


if (point.i < point.n.x + 2 && point.k < point.n.z + 2) {
    int jg = -1; NOT_UNUSED(jg);      
#line 310 "/home/pwachara/basilisk/src/embed.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_y = Delta;

  double Delta_z = Delta;


  double Delta_x = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_y);

  NOT_UNUSED(Delta_z);


  NOT_UNUSED(Delta_x);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 310 "/home/pwachara/basilisk/src/embed.h"
if (val(s.y,0,0,0) && ((!val(c,0,0,0) || !val(c,0,-1,0)) || val(s.y,0,0,0) < smin))
 val(s.y,0,0,0) = 0.;}
  
#line 396 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
  


if (point.i < point.n.x + 2 && point.j < point.n.y + 2) {
    int kg = -1; NOT_UNUSED(kg);      
#line 310 "/home/pwachara/basilisk/src/embed.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_z = Delta;

  double Delta_x = Delta;


  double Delta_y = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_z);

  NOT_UNUSED(Delta_x);


  NOT_UNUSED(Delta_y);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 310 "/home/pwachara/basilisk/src/embed.h"
if (val(s.z,0,0,0) && ((!val(c,0,0,0) || !val(c,0,0,-1)) || val(s.z,0,0,0) < smin))
 val(s.z,0,0,0) = 0.;}
  
#line 403 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
#line 311 "/home/pwachara/basilisk/src/embed.h"
}
      
#line 317 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
  }
}    
#line 313 "/home/pwachara/basilisk/src/embed.h"
changed = 0;    
#line 314 "/home/pwachara/basilisk/src/grid/stencils.h"
{
    static int _first = 1.;
    ForeachData _loop = {
      .fname = "/home/pwachara/basilisk/src/embed.h", .line = 314, .first = _first
    };
    if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
 _attribute[s.i].input = _attribute[s.i].output = _attribute[s.i].nowarning = false;
 _attribute[s.i].width = 0;
      }
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point);      
#line 315 "/home/pwachara/basilisk/src/embed.h"
{_stencil_val(c,0,0,0); _stencil_val(c,0,0,0); {   
 
  {
   for (int i = 0; i <= 1; i++)
     {_stencil_val(s.x,i,0,0);
          } 









_stencil_val(s.x,0,0,0); _stencil_val(s.x,1,0,0);
     {_stencil_val_a(c,0,0,0);   }
#line 330 "/home/pwachara/basilisk/src/embed.h"
          
 
} 
#line 317
{
   for (int i = 0; i <= 1; i++)
     {_stencil_val(s.y,0,i,0);
          } 









_stencil_val(s.y,0,0,0); _stencil_val(s.y,0,1,0);
     {_stencil_val_a(c,0,0,0);   }
#line 330 "/home/pwachara/basilisk/src/embed.h"
          
 
} 
#line 317
{
   for (int i = 0; i <= 1; i++)
     {_stencil_val(s.z,0,0,i);
          } 









_stencil_val(s.z,0,0,0); _stencil_val(s.z,0,0,1);
     {_stencil_val_a(c,0,0,0);   }
#line 330 "/home/pwachara/basilisk/src/embed.h"
          
 
}
   







{_stencil_val_a(c,0,0,0);   }      
}      }    
#line 328 "/home/pwachara/basilisk/src/grid/stencils.h"
check_stencil (&_loop);
    boundary_stencil (&_loop);
    _first = 0;
  } 
#line 265 "/home/pwachara/basilisk/src/grid/multigrid.h"
{
  OMP_PARALLEL (reduction(+:changed)) {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    point.n.x = point.n.y = point.n.z = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
      for (_k = 2; _k < point.n.x + 2; _k++) {
 point.i = _k;

 for (point.j = 2; point.j < point.n.y + 2; point.j++)

   for (point.k = 2; point.k < point.n.z + 2; point.k++)
      
#line 315 "/home/pwachara/basilisk/src/embed.h"
{  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 315 "/home/pwachara/basilisk/src/embed.h"
if (val(c,0,0,0) > 0. && val(c,0,0,0) < 1.) {
 int n = 0;
  {
   for (int i = 0; i <= 1; i++)
     if (val(s.x,i,0,0) > 0.)
       n++;
#line 330 "/home/pwachara/basilisk/src/embed.h"
   if (opposite && val(s.x,0,0,0) == 0. && val(s.x,1,0,0) == 0.)
     val(c,0,0,0) = 0., changed++;
 } 
#line 317
{
   for (int i = 0; i <= 1; i++)
     if (val(s.y,0,i,0) > 0.)
       n++;
#line 330 "/home/pwachara/basilisk/src/embed.h"
   if (opposite && val(s.y,0,0,0) == 0. && val(s.y,0,1,0) == 0.)
     val(c,0,0,0) = 0., changed++;
 } 
#line 317
{
   for (int i = 0; i <= 1; i++)
     if (val(s.z,0,0,i) > 0.)
       n++;
#line 330 "/home/pwachara/basilisk/src/embed.h"
   if (opposite && val(s.z,0,0,0) == 0. && val(s.z,0,0,1) == 0.)
     val(c,0,0,0) = 0., changed++;
 }







 if (n < 3)
   val(c,0,0,0) = 0., changed++;
      }}      
#line 282 "/home/pwachara/basilisk/src/grid/multigrid.h"
}
  }
}
#line 342 "/home/pwachara/basilisk/src/embed.h"
{mpi_all_reduce_array(&changed,MPI_INT,MPI_SUM,1);}

    schanged += changed;
  }
  if (changed)
    fprintf (ferr, "src/embed.h:%d: warning: fractions_cleanup() did not converge after "
      "%d iterations\n", 348, i);
  {end_tracing("fractions_cleanup","/home/pwachara/basilisk/src/embed.h",349);return schanged;}
end_tracing("fractions_cleanup","/home/pwachara/basilisk/src/embed.h",350);}
#line 374 "/home/pwachara/basilisk/src/embed.h"

static inline double dirichlet_gradient_x (Point point, scalar s, scalar cs,
        coord n, coord p, double bc,
        double * coef)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
  
    
#line 380 "/home/pwachara/basilisk/src/embed.h"
n.x = - n.x;    n.y = - n.y;    n.z = - n.z;
  double d[2], v[2] = {1e30f,1e30f};
  bool defined = true;
  
    if (defined && !val(fs.x,(n.x > 0.),0,0))
      defined = false;    
#line 384
if (defined && !val(fs.y,0,(n.y > 0.),0))
      defined = false;    
#line 384
if (defined && !val(fs.z,0,0,(n.z > 0.)))
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*( (int)((n.x) > 0 ? 1 : -1));
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;





      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = val(fs.x,i + (i < 0),j,k);
      for (int m = -1; m <= 1 && defined; m++)
 if (!val(fs.y,i,j,k+m) || !val(fs.y,i,j+1,k+m) ||
     !val(fs.z,i,j+m,k) || !val(fs.z,i,j+m,k+1) ||
     !val(cs,i,j+m,k-1) || !val(cs,i,j+m,k) || !val(cs,i,j+m,k+1))
   defined = false;
      if (defined)

 v[l] =
   (((((((val(s,i,j-1,k-1)))*((y1) - 1.) + ((val(s,i,j+1,k-1)))*((y1) + 1.))*(y1)/2. - ((val(s,i,j,k-1)))*((y1) - 1.)*((y1) + 1.)))*((z) - 1.) + (((((val(s,i,j-1,k+1)))*((y1) - 1.) + ((val(s,i,j+1,k+1)))*((y1) + 1.))*(y1)/2. - ((val(s,i,j,k+1)))*((y1) - 1.)*((y1) + 1.)))*((z) + 1.))*(z)/2. - (((((val(s,i,j-1,k)))*((y1) - 1.) + ((val(s,i,j+1,k)))*((y1) + 1.))*(y1)/2. - ((val(s,i,j,k)))*((y1) - 1.)*((y1) + 1.)))*((z) - 1.)*((z) + 1.))





                                                  ;

      else
 break;
    }
  if (v[0] == 1e30f) {





    d[0] = ( 1e-3 > (fabs(p.x/n.x)) ? 1e-3 : (fabs(p.x/n.x)));
    *coef = - 1./(d[0]*Delta);
    return bc/(d[0]*Delta);
  }





  *coef = 0.;
  if (v[1] != 1e30f)
    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
  return (bc - v[0])/(d[0]*Delta);
}
#line 375
static inline double dirichlet_gradient_y (Point point, scalar s, scalar cs,
        coord n, coord p, double bc,
        double * coef)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_y = Delta;

  double Delta_z = Delta;


  double Delta_x = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_y);

  NOT_UNUSED(Delta_z);


  NOT_UNUSED(Delta_x);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
  
    
#line 380 "/home/pwachara/basilisk/src/embed.h"
n.y = - n.y;    n.z = - n.z;    n.x = - n.x;
  double d[2], v[2] = {1e30f,1e30f};
  bool defined = true;
  
    if (defined && !val(fs.y,0,(n.y > 0.),0))
      defined = false;    
#line 384
if (defined && !val(fs.z,0,0,(n.z > 0.)))
      defined = false;    
#line 384
if (defined && !val(fs.x,(n.x > 0.),0,0))
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*( (int)((n.y) > 0 ? 1 : -1));
      d[l] = (i - p.y)/n.y;
      double y1 = p.z + d[l]*n.z;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;





      double z = p.x + d[l]*n.x;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = val(fs.y,k,i + (i < 0),j);
      for (int m = -1; m <= 1 && defined; m++)
 if (!val(fs.z,k+m,i,j) || !val(fs.z,k+m,i,j+1) ||
     !val(fs.x,k,i,j+m) || !val(fs.x,k+1,i,j+m) ||
     !val(cs,k-1,i,j+m) || !val(cs,k,i,j+m) || !val(cs,k+1,i,j+m))
   defined = false;
      if (defined)

 v[l] =
   (((((((val(s,k-1,i,j-1)))*((y1) - 1.) + ((val(s,k-1,i,j+1)))*((y1) + 1.))*(y1)/2. - ((val(s,k-1,i,j)))*((y1) - 1.)*((y1) + 1.)))*((z) - 1.) + (((((val(s,k+1,i,j-1)))*((y1) - 1.) + ((val(s,k+1,i,j+1)))*((y1) + 1.))*(y1)/2. - ((val(s,k+1,i,j)))*((y1) - 1.)*((y1) + 1.)))*((z) + 1.))*(z)/2. - (((((val(s,k,i,j-1)))*((y1) - 1.) + ((val(s,k,i,j+1)))*((y1) + 1.))*(y1)/2. - ((val(s,k,i,j)))*((y1) - 1.)*((y1) + 1.)))*((z) - 1.)*((z) + 1.))





                                                  ;

      else
 break;
    }
  if (v[0] == 1e30f) {





    d[0] = ( 1e-3 > (fabs(p.y/n.y)) ? 1e-3 : (fabs(p.y/n.y)));
    *coef = - 1./(d[0]*Delta);
    return bc/(d[0]*Delta);
  }





  *coef = 0.;
  if (v[1] != 1e30f)
    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
  return (bc - v[0])/(d[0]*Delta);
}
#line 375
static inline double dirichlet_gradient_z (Point point, scalar s, scalar cs,
        coord n, coord p, double bc,
        double * coef)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_z = Delta;

  double Delta_x = Delta;


  double Delta_y = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_z);

  NOT_UNUSED(Delta_x);


  NOT_UNUSED(Delta_y);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
  
    
#line 380 "/home/pwachara/basilisk/src/embed.h"
n.z = - n.z;    n.x = - n.x;    n.y = - n.y;
  double d[2], v[2] = {1e30f,1e30f};
  bool defined = true;
  
    if (defined && !val(fs.z,0,0,(n.z > 0.)))
      defined = false;    
#line 384
if (defined && !val(fs.x,(n.x > 0.),0,0))
      defined = false;    
#line 384
if (defined && !val(fs.y,0,(n.y > 0.),0))
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*( (int)((n.z) > 0 ? 1 : -1));
      d[l] = (i - p.z)/n.z;
      double y1 = p.x + d[l]*n.x;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;





      double z = p.y + d[l]*n.y;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = val(fs.z,j,k,i + (i < 0));
      for (int m = -1; m <= 1 && defined; m++)
 if (!val(fs.x,j,k+m,i) || !val(fs.x,j+1,k+m,i) ||
     !val(fs.y,j+m,k,i) || !val(fs.y,j+m,k+1,i) ||
     !val(cs,j+m,k-1,i) || !val(cs,j+m,k,i) || !val(cs,j+m,k+1,i))
   defined = false;
      if (defined)

 v[l] =
   (((((((val(s,j-1,k-1,i)))*((y1) - 1.) + ((val(s,j+1,k-1,i)))*((y1) + 1.))*(y1)/2. - ((val(s,j,k-1,i)))*((y1) - 1.)*((y1) + 1.)))*((z) - 1.) + (((((val(s,j-1,k+1,i)))*((y1) - 1.) + ((val(s,j+1,k+1,i)))*((y1) + 1.))*(y1)/2. - ((val(s,j,k+1,i)))*((y1) - 1.)*((y1) + 1.)))*((z) + 1.))*(z)/2. - (((((val(s,j-1,k,i)))*((y1) - 1.) + ((val(s,j+1,k,i)))*((y1) + 1.))*(y1)/2. - ((val(s,j,k,i)))*((y1) - 1.)*((y1) + 1.)))*((z) - 1.)*((z) + 1.))





                                                  ;

      else
 break;
    }
  if (v[0] == 1e30f) {





    d[0] = ( 1e-3 > (fabs(p.z/n.z)) ? 1e-3 : (fabs(p.z/n.z)));
    *coef = - 1./(d[0]*Delta);
    return bc/(d[0]*Delta);
  }





  *coef = 0.;
  if (v[1] != 1e30f)
    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
  return (bc - v[0])/(d[0]*Delta);
}
#line 374 "/home/pwachara/basilisk/src/embed.h"

static void _stencil_dirichlet_gradient_x (Point point, scalar s, scalar cs,
_stencil_undefined * n,_stencil_undefined * p,_stencil_undefined * bc,
_stencil_undefined * coef)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;       
  
       
  
  
  
    
#line 384 "/home/pwachara/basilisk/src/embed.h"
{_stencil_val(fs.x,o_stencil,0,0  ); 
          }    
#line 384
{_stencil_val(fs.y,0,o_stencil,0  ); 
          }    
#line 384
{_stencil_val(fs.z,0,0,o_stencil  ); 
          }
    
for (int l = 0; l <= 1; l++) {                                                      
       
         
      
      
        





      
      
        
       _stencil_val(fs.x,    o_stencil,o_stencil,o_stencil);   
             
 {_stencil_val(fs.y,o_stencil,o_stencil,o_stencil);_stencil_val(fs.y,o_stencil,o_stencil,o_stencil);
_stencil_val(fs.z,o_stencil,o_stencil,o_stencil);_stencil_val(fs.z,o_stencil,o_stencil,o_stencil);
_stencil_val(cs,o_stencil,o_stencil,o_stencil);_stencil_val(cs,o_stencil,o_stencil,o_stencil);_stencil_val(cs,o_stencil,o_stencil,o_stencil);         
}
{

 {
_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);                                                  





}
 

}    
}         
     
   
     





   
     
  







return   ;
}
#line 375
static void _stencil_dirichlet_gradient_y (Point point, scalar s, scalar cs,
_stencil_undefined * n,_stencil_undefined * p,_stencil_undefined * bc,
_stencil_undefined * coef)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_y = Delta;

  double Delta_z = Delta;


  double Delta_x = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_y);

  NOT_UNUSED(Delta_z);


  NOT_UNUSED(Delta_x);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;       
  
       
  
  
  
    
#line 384 "/home/pwachara/basilisk/src/embed.h"
{_stencil_val(fs.y,0,o_stencil,0  ); 
          }    
#line 384
{_stencil_val(fs.z,0,0,o_stencil  ); 
          }    
#line 384
{_stencil_val(fs.x,o_stencil,0,0  ); 
          }
    
for (int l = 0; l <= 1; l++) {                                                      
       
         
      
      
        





      
      
        
       _stencil_val(fs.y,o_stencil,    o_stencil,o_stencil);   
             
 {_stencil_val(fs.z,o_stencil,o_stencil,o_stencil);_stencil_val(fs.z,o_stencil,o_stencil,o_stencil);
_stencil_val(fs.x,o_stencil,o_stencil,o_stencil);_stencil_val(fs.x,o_stencil,o_stencil,o_stencil);
_stencil_val(cs,o_stencil,o_stencil,o_stencil);_stencil_val(cs,o_stencil,o_stencil,o_stencil);_stencil_val(cs,o_stencil,o_stencil,o_stencil);         
}
{

 {
_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);                                                  





}
 

}    
}         
     
   
     





   
     
  







return   ;
}
#line 375
static void _stencil_dirichlet_gradient_z (Point point, scalar s, scalar cs,
_stencil_undefined * n,_stencil_undefined * p,_stencil_undefined * bc,
_stencil_undefined * coef)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_z = Delta;

  double Delta_x = Delta;


  double Delta_y = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_z);

  NOT_UNUSED(Delta_x);


  NOT_UNUSED(Delta_y);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;       
  
       
  
  
  
    
#line 384 "/home/pwachara/basilisk/src/embed.h"
{_stencil_val(fs.z,0,0,o_stencil  ); 
          }    
#line 384
{_stencil_val(fs.x,o_stencil,0,0  ); 
          }    
#line 384
{_stencil_val(fs.y,0,o_stencil,0  ); 
          }
    
for (int l = 0; l <= 1; l++) {                                                      
       
         
      
      
        





      
      
        
       _stencil_val(fs.z,o_stencil,o_stencil,    o_stencil);   
             
 {_stencil_val(fs.x,o_stencil,o_stencil,o_stencil);_stencil_val(fs.x,o_stencil,o_stencil,o_stencil);
_stencil_val(fs.y,o_stencil,o_stencil,o_stencil);_stencil_val(fs.y,o_stencil,o_stencil,o_stencil);
_stencil_val(cs,o_stencil,o_stencil,o_stencil);_stencil_val(cs,o_stencil,o_stencil,o_stencil);_stencil_val(cs,o_stencil,o_stencil,o_stencil);         
}
{

 {
_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);_stencil_val(s,o_stencil,o_stencil,o_stencil);                                                  





}
 

}    
}         
     
   
     





   
     
  







return   ;
}

double dirichlet_gradient (Point point, scalar s, scalar cs,
      coord n, coord p, double bc, double * coef)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;





  
#line 450 "/home/pwachara/basilisk/src/embed.h"
if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return dirichlet_gradient_x (point, s, cs, n, p, bc, coef);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return dirichlet_gradient_y (point, s, cs, n, p, bc, coef);
  return dirichlet_gradient_z (point, s, cs, n, p, bc, coef);

  return 1e30f;
}
#line 442
static void _stencil_dirichlet_gradient (Point point, scalar s, scalar cs,
_stencil_undefined * n,_stencil_undefined * p,_stencil_undefined * bc,_stencil_undefined * coef)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
#line 450 "/home/pwachara/basilisk/src/embed.h"
{ {
      
{ _stencil_dirichlet_gradient_x (point, s, cs,NULL ,NULL ,NULL ,NULL );}  
}
    
{ _stencil_dirichlet_gradient_y (point, s, cs,NULL ,NULL ,NULL ,NULL );}} 
_stencil_dirichlet_gradient_z (point, s, cs,NULL ,NULL ,NULL ,NULL );  return;

  return ;
}

bid embed;
#line 470 "/home/pwachara/basilisk/src/embed.h"
static inline
coord embed_gradient (Point point, vector u, coord p, coord n)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2;
  
#line 473 "/home/pwachara/basilisk/src/embed.h"
coord dudn;
   {
    bool dirichlet = false;
    double vb = _attribute[u.x.i].boundary[embed] (point, point, u.x, &dirichlet);
    if (dirichlet) {
      double val;
      dudn.x = dirichlet_gradient (point, u.x, cs, n, p, vb, &val);
    }
    else
      dudn.x = vb;
    if (dudn.x == 1e30f)
      dudn.x = 0.;
  } 
#line 474
{
    bool dirichlet = false;
    double vb = _attribute[u.y.i].boundary[embed] (point, point, u.y, &dirichlet);
    if (dirichlet) {
      double val;
      dudn.y = dirichlet_gradient (point, u.y, cs, n, p, vb, &val);
    }
    else
      dudn.y = vb;
    if (dudn.y == 1e30f)
      dudn.y = 0.;
  } 
#line 474
{
    bool dirichlet = false;
    double vb = _attribute[u.z.i].boundary[embed] (point, point, u.z, &dirichlet);
    if (dirichlet) {
      double val;
      dudn.z = dirichlet_gradient (point, u.z, cs, n, p, vb, &val);
    }
    else
      dudn.z = vb;
    if (dudn.z == 1e30f)
      dudn.z = 0.;
  }
  return dudn;
}
#line 470 "/home/pwachara/basilisk/src/embed.h"
static void 
_stencil_embed_gradient (Point point, vector u,_stencil_undefined * p,_stencil_undefined * n)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);  
#line 3 "/home/pwachara/basilisk/src/grid/variables.h"
double Delta = L0*(1./((1 << point.level)*Dimensions_scale));
  double Delta_x = Delta;

  double Delta_y = Delta;


  double Delta_z = Delta;


  double x = ((ig + 1)/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)))*Delta + X0; NOT_UNUSED(x);

  double y = ((jg + 1)/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)))*Delta + Y0;



  NOT_UNUSED(y);

  double z = ((kg + 1)/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)))*Delta + Z0;



  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);


  NOT_UNUSED(Delta_z);


  ;
  
#line 206 "/home/pwachara/basilisk/src/grid/multigrid.h"
int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + 2)%2) - 1,
    2*((point.j + 2)%2) - 1,
    2*((point.k + 2)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;
  parent.j = (point.j + 2)/2;
  parent.k = (point.k + 2)/2; 
  
   
#line 474 "/home/pwachara/basilisk/src/embed.h"
{     
    
    default_stencil ( point,((scalar[]){ u.x,{-1}}) ); 
{ 
       
_stencil_dirichlet_gradient (point, u.x, cs,NULL ,NULL ,NULL ,NULL );    
}  




} 
#line 474
{     
    
    default_stencil ( point,((scalar[]){ u.y,{-1}}) ); 
{ 
       
_stencil_dirichlet_gradient (point, u.y, cs,NULL ,NULL ,NULL ,NULL );    
}  




} 
#line 474
{     
    
    default_stencil ( point,((scalar[]){ u.z,{-1}}) ); 
{ 
       
_stencil_dirichlet_gradient (point, u.z, cs,NULL ,NULL ,NULL ,NULL );    
}  




}
  return ;
}
#line 508 "/home/pwachara/basilisk/src/embed.h"
     
void embed_force (scalar p, vector u, vecto