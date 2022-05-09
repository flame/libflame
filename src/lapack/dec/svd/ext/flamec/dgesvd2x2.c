/*
    Copyright (c) 2021 Advanced Micro Devices, Inc. All rights reserved.
    Aug 4, 2021
*/

#include "FLAME.h"
#include "FLA_f2c.h" 
        

/* Subroutine */ int dgesvd2x2(char *jobu, char *jobvt, integer *m, integer *n, 
	doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
	ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, 
	integer *info)
{
   int i__2;       
   double AAT[4], ATA[4], tempu[4], tempvt[4], phi, cphi, sphi, theta, ctheta, stheta, s11, s22, temp, susum, sudiff, s1, s2;    
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
           *lwork = 10;
           work[0] = (doublereal) *lwork;
           return 0;
        }
        else { 
           *lwork = 4;
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


   //S calculation
   s[0] = (sqrt(pow(a[0] - a[*lda+1], 2) + pow(a[1] + a[*lda], 2)) + sqrt(pow(a[0] + a[*lda+1], 2) + pow(a[1] - a[*lda], 2))) / 2;
   s[1] = fabs(s[0] - sqrt(pow(a[0] - a[*lda+1], 2) + pow(a[1] + a[*lda], 2)));  
   
  //V calculation
   tempvt[2] = (s[0] > s[1]) ? sin((atan2(2 * (a[0] * a[1] + a[*lda] * a[*lda+1]), a[0] * a[0] - a[1] * a[1] + a[*lda] * a[*lda] - a[*lda+1] * a[*lda+1])) / 2) : 0;   // atan2 =0 will never occur as in that case s[0]=s[1]
   tempvt[0] = sqrt(1 - tempvt[2] * tempvt[2]);
   tempvt[1] = -tempvt[2];
   tempvt[3] = tempvt[0];
  
  //U calculation
   tempu[0] = (s[0] != 0) ? (a[0] * tempvt[0] + a[1] * tempvt[2]) / s[0] : 1;
   tempu[1] = (s[0] != 0) ? (a[*lda] * tempvt[0] + a[*lda+1] * tempvt[2]) / s[0] : 0;
   tempu[2] = (s[1] != 0) ? (a[0] * tempvt[1] + a[1] * tempvt[3]) / s[1] : -tempu[1];
   tempu[3] = (s[1] != 0) ? (a[*lda] * tempvt[1] + a[*lda+1] * tempvt[3]) / s[1] : tempu[0];

   if (*jobvt == 'A' || *jobvt == 'S') {
       vt[0] = tempvt[0], vt[1] = tempvt[1], vt[*ldvt] = tempvt[2], vt[*ldvt+1] = tempvt[3];
   }
   // U calculation (to be read as column major)
      //assignment for U matrix 
   if (*jobu == 'A' || *jobu == 'S') {  
       u[0] = tempu[0], u[1] = tempu[1], u[*ldu] = tempu[2], u[*ldu + 1] = tempu[3];
   }
   //Output A matrix contains sigma in its diagonal and rest values as 0 ( exception when jobu or jobv is O then A is overwritten)
   a[0] = s[0], a[1]=0, a[*lda] = 0, a[*lda + 1] = s[1];
   //handling the case when either jobu or jobvt is O	
   if (*jobu == 'O') {
       a[0] = tempu[0], a[1]= tempu[1], a[*lda] = tempu[*ldu], a[*lda + 1] = tempu[*ldu + 1];
   } else if (*jobvt == 'O') {  
       a[0]= tempvt[0], a[1] = tempvt[1], a[*lda] = tempvt[*ldvt], a[*lda + 1] = tempvt[*ldvt + 1];
   }
   return 0;
}
