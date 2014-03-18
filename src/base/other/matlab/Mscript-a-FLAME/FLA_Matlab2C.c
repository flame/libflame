/* ====================================================================================== */
/*
   C file containing functions to convert things from Matlab FLAME to C FLAME.
*/
/* ====================================================================================== */


#include "mex.h"
#include "FLA_Matlab2C.h"
#include <string.h>
#include <stdio.h>


/* Check if the number of arguments with which the function is called is correct */
/* rrhs: required number of arguments */
/* srhs: supplied number of arguments */
void FLA_M2C_CheckNumArgs(int rrhs, int srhs) {
  char err[200];
  
  /* rrhs is the required number of arguments.  If an additional optional argument is
     specified, the routine is timed and the timing information returned in that argument. */
  if (srhs != rrhs && srhs != rrhs+1) {
    sprintf(err, "Wrong number of arguments: %d arguments supplied; %d or %d arguments required", srhs, rrhs, rrhs+1);
    mexErrMsgTxt(err);
  }
  
}


/* Convert FLAME attributes specified as a string in MATLAB to the numerical values 
   required when calling FLAME C functions. */
int FLA_M2C_ConvertAttribute(const mxArray *mobj) {
  char *attr;
  int value;

  if (!mxIsChar(mobj)) {
    mexErrMsgTxt("Char array expected.");
  }
  else {
    /* Convert to C style string */
    attr = mxArrayToString(mobj);

    if (strcmp(attr, "FLA_NO_TRANSPOSE") == 0)
      value = FLA_NO_TRANSPOSE;
    else if (strcmp(attr, "FLA_TRANSPOSE") == 0)
      value = FLA_TRANSPOSE;
    else if (strcmp(attr, "FLA_CONJUGATE") == 0)
      value = FLA_CONJUGATE;
    else if (strcmp(attr, "FLA_NO_CONJUGATE") == 0)
      value = FLA_NO_CONJUGATE;
    else if (strcmp(attr, "FLA_CONJ_TRANSPOSE") == 0)
      value = FLA_CONJ_TRANSPOSE;
    else if (strcmp(attr, "FLA_LOWER_TRIANGULAR") == 0)
      value = FLA_LOWER_TRIANGULAR;
    else if (strcmp(attr, "FLA_UPPER_TRIANGULAR") == 0)
      value = FLA_UPPER_TRIANGULAR;
    else if (strcmp(attr, "FLA_NONUNIT_DIAG") == 0)
      value = FLA_NONUNIT_DIAG;
    else if (strcmp(attr, "FLA_UNIT_DIAG") == 0)
      value = FLA_UNIT_DIAG;
    else if (strcmp(attr, "FLA_ZERO_DIAG") == 0)
      value = FLA_ZERO_DIAG;
    else if (strcmp(attr, "FLA_LEFT") == 0)
      value = FLA_LEFT;
    else if (strcmp(attr, "FLA_RIGHT") == 0)
      value = FLA_RIGHT;
    else
      value = -1;

    mxFree(attr);
    
    return value;
  }

}
    

/* Convert Matlab array to C FLAME object */
FLA_Obj FLA_M2C_ConvertMxArray(const mxArray *mobj) {
  int type;
  FLA_Obj fobj;

  if (!mxIsNumeric(mobj)) {
    mexErrMsgTxt("Numeric data expected.");
  }
  else {
    /* Determine data type of the matlab array. */
    if (mxIsDouble(mobj))
      if (mxIsComplex(mobj))
	type = FLA_DOUBLE_COMPLEX;
      else
	type = FLA_DOUBLE;
    else if (mxIsSingle(mobj))
      if (mxIsComplex(mobj))
	type = FLA_COMPLEX;
      else
	type = FLA_FLOAT;
    else if (mxIsInt32(mobj))
      type = FLA_INT;
    else
      mexErrMsgTxt("Data type not supported.");

    fobj.datatype = type;
    fobj.m = mxGetM(mobj);
    fobj.n = mxGetN(mobj);
    fobj.buffer = mxGetPr(mobj);
    fobj.ldim = fobj.m;

    return fobj;
  }

}


/* Convert Matlab array to double */
double *FLA_M2C_ConvertDoublePtr(const mxArray *mobj) {
  double *var = 0;

  if (!mxIsDouble(mobj)) {
    mexErrMsgTxt("Double precision variable expected.");
  }
  else if (mxGetM(mobj) != 1 || mxGetN(mobj) != 1) {
    mexErrMsgTxt("Scalar variable expected.");
  }
  else {
    var = mxGetPr(mobj);
  }

  return var;
}


/* Convert the attribute and FLA_Obj arguments */
void FLA_M2C_ConvertArgs(int nrhs, const mxArray *prhs[], int nint, 
			 int attr[], FLA_Obj obj[]) {
  char err[200];
  int i, j;

  /* First convert all the attributes */
  /* Check if the argument is a string before calling the conversion function */
  for (i=0; i<nint; i++) {
    if (mxIsChar(prhs[i])) {
      attr[i] = FLA_M2C_ConvertAttribute(prhs[i]);
    }
    else {
      sprintf(err, "The %dth argument must be a FLAME attribute string.", i+1);
      mexErrMsgTxt(err);
    }
  }

  /* The remaining arguments are converted to FLA_Obj */
  /* Check if the argument is numeric before calling the conversion function */
  for (i=i; i<nrhs; i++) {
    if (mxIsNumeric(prhs[i])) {
      obj[i-nint] = FLA_M2C_ConvertMxArray(prhs[i]);
    }
    else {
      sprintf(err, "The %dth argument must be numeric.", i+1);
      mexErrMsgTxt(err);
    }
  }

} 
