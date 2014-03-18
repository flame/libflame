#include "FLAME.h"

void FLA_M2C_CheckNumArgs(int rrhs, int srhs);
int FLA_M2C_ConvertAttribute(const mxArray *mobj);
FLA_Obj FLA_M2C_ConvertMxArray(const mxArray *mobj);
double *FLA_M2C_ConvertDoublePtr(const mxArray *mobj);
void FLA_M2C_ConvertArgs(int nrhs, const mxArray *prhs[], int nint, 
			 int attr[], FLA_Obj obj[]);
