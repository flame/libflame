/*
    Copyright (c) 2021 Advanced Micro Devices, Inc. All rights reserved.
    Aug 4, 2021
*/

#include "FLAME.h"
#include "FLA_f2c.h" 
void multaat(double* AA, const double*  A, int lda);
void multata(double* AA, const double*  A, int lda);
//Helper functions to calculate A*Atranspose and Atranspose*A
void multaat(double* AA, const double*  A, int lda)
{
        AA[0] = A[0] * A[0] + A[1] * A[1];
        AA[1] = A[0] * A[lda] + A[1] * A[lda + 1];
        AA[lda] = A[lda] * A[0] + A[lda + 1] * A[1];
        AA[lda + 1] = A[lda] * A[lda] + A[lda + 1] * A[lda + 1];
}

void multata(double * AA, const double * A, int lda)
{
        AA[0] = A[0] * A[0] + A[lda] * A[lda];
        AA[1] = A[0] * A[1] + A[lda] * A[lda + 1];
        AA[lda] = A[1] * A[0] + A[lda + 1] * A[lda];
        AA[lda + 1] = A[1] * A[1] + A[lda + 1] * A[lda + 1];
}
//Sign calculator (C) to calculate Vtranspose
double signs(double val)
{
        return val > 0 ? 1 : (val == 0 ? 0 : -1);
}
        

/* Subroutine */ int dgesvd2x2(char *jobu, char *jobvt, integer *m, integer *n, 
	doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
	ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, 
	integer *info)
{
   int i__2;       
   double AAT[4], ATA[4], phi, cphi, sphi, theta, ctheta, stheta, s11, s22, temp, susum, sudiff, s1, s2;    
   *info = 0;
 #if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "dgesvd inputs: jobu %c, jobvt %c, m %d, n %d, lda %d, ldu %d, ldvt %d\n", *jobu, *jobvt, *m, *n, *lda, *ldu, *ldvt);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
 #endif
   // Error checking , here m, n sizes will always be valid because of specific size
   if (! (*jobu == 'A' || *jobu == 'S' || *jobu == 'O' || *jobu == 'N')) {
       *info = -1;
   } else if (! (*jobvt == 'A' || *jobvt == 'S' || *jobvt == 'O' || *jobvt == 'N')) { 
       *info = -2;   
   } else if (*lda < 2){
       *info = -3;
   } else if (*ldu < 1 || ( (*jobu == 'A' && *ldu < 2)  || (*jobu == 'S' && *ldu < 2))) {
       *info = -9;
   } else if (*ldvt < 1 || ( (*jobvt == 'A' && *ldvt < 2)  || (*jobvt == 'S' && *ldvt < 2))) {
       *info = -11;
   }
   //special case which is not allowed
   if(*jobu == 'O' && *jobvt == 'O') {
         *info = -2;
   } 
   if (*info == 0) {
      if (*lwork == -1) {    //lwork -1 means workspace assignment though not required
        if ( *jobvt == 'O' || *jobu == 'O'){
           *lwork=10;
           work[0] = (doublereal) *lwork;
           return 0;
        }
        else { 
           *lwork=4;
           work[0] = (doublereal) *lwork;
           return 0;
        }
      }
   }
   //Print error
   if (*info != 0) {
      i__2 = -(*info);
      xerbla_("DGESVD", &i__2);
      return 0;
   } 
   //To account for column major(taking transpose)
   temp = a[1]; 
   a[1] = a[*lda];
   a[*lda] = temp;
   // U calculation (to be read as column major)
   multaat( (double*) AAT, (double*)a, *lda);
   phi= 0.5 * atan2(AAT[1] + AAT[*lda], AAT[0] - AAT[*lda + 1]);
   cphi = cos(phi);
   sphi= sin(phi);
  	
   if (*jobu != 'N') { //If JobU=N , U is not calculated
      //assignment for U matrix 
      u[0] = cphi, u[1] = sphi, u[*ldu] = -sphi, u[*ldu + 1] = cphi;
   }
   //VT matrix calculation	
   multata(ATA,(double*)a, *lda);
   if (*jobvt != 'N') {
      theta = 0.5 * atan2(ATA[1] + ATA[*lda], ATA[0] - ATA[*lda + 1]);
      ctheta = cos(theta);
      stheta = sin(theta);
      //s11 and s22 are only used to calculate the signs lost due to squaring for VT
      s11 = ( a[0] * cphi + a[*lda] * sphi) * ctheta + ( a[1] * cphi + a[*lda + 1] * sphi) * stheta;
      s22 = ( a[0]*sphi - a[*lda] * cphi ) * stheta + (-a[1] * sphi + a[*lda+1] * cphi) * ctheta;
      //VT assignment
      vt[0] = signs(s11) * ctheta, vt[1] = -signs(s22) * stheta, vt[*ldvt] = signs(s11) * stheta, vt[*ldvt + 1] = signs(s22) * ctheta;	
   }
   susum = AAT[0] + AAT[*lda + 1];
   sudiff = AAT[0] - AAT[*lda + 1];
   s1 = susum;
   s2 = sqrt(sudiff * sudiff + 4 * AAT[1] * AAT[*lda]);
   s[0] = sqrt( (s1 + s2) * 0.5);
   s[1] = sqrt( (s1 - s2) * 0.5);
   //Output A matrix contains sigma in its diagonal and rest values as 0 ( exception when jobu or jobv is O then A is overwritten)
   a[0] = s[0], a[1]=0, a[*lda] = 0, a[*lda + 1] = s[1];
   //handling the case when either jobu or jobvt is O	
   if (*jobu == 'O') {
      a[0] = u[0], a[1]= u[1], a[*lda] = u[*ldu], a[*lda + 1] = u[*ldu + 1];
      //U is reset		
      u[0] = 0, u[1] = 0, u[*ldu] = 0, u[*ldu + 1] = 0;
   } else if (*jobvt == 'O') {  
      a[0]= vt[0], a[1] = vt[1], a[*lda] = vt[*ldvt], a[*lda + 1] = vt[*ldvt + 1];
      //VT is reset
      vt[0] = 0, vt[1] = 0, vt[*ldvt] = 0, vt[*ldvt + 1] = 0;
   }
return 0;
}
