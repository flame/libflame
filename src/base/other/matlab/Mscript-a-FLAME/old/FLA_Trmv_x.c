/* Mex file created by interface.pl for the FLAME C function FLA_Trmvsx_external() */

#include "mex.h"
#include "FLA_Matlab2C.h"

extern double FLA_Clock();

#define NRHS 8
#define NINT 3

#define NOBJ (NRHS-NINT)


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int attr[NINT];
  FLA_Obj obj[NOBJ];
  double *dtime;

  FLA_Init();

  /* Check if the number of arguments supplied is correct */
  FLA_M2C_CheckNumArgs(NRHS, nrhs);

  /* Convert Matlab arguments into the appropriate FLAME C arguments */
  FLA_M2C_ConvertArgs(NRHS, prhs, NINT, attr, obj);

  /* If an extra argument is supplied, collect timing informaion in it. */
  if (nrhs == NRHS+1)
    dtime = FLA_M2C_ConvertDoublePtr(prhs[NRHS]);

  /* Now call the C FLAME function, timing it if the extra argument is given. */

  if (nrhs == NRHS+1)
    *dtime = FLA_Clock();

  FLA_Trmvsx_external(attr[0], attr[1], attr[2], obj[0], obj[1], obj[2], obj[3], obj[4]);

  if (nrhs == NRHS+1)
    *dtime = FLA_Clock() - *dtime;

  FLA_Finalize();

}

