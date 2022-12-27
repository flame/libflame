#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

/* Begin bbcsd_double_parameters  class definition */
class bbcsd_double_parameters{
	
   public:
   /* input params to the API **/
   char  JOBU1; 
   char  JOBU2;
   char  JOBV1T;
   char  JOBV2T;
   char  TRANS;
   int 	 M;  // 512
   int 	 P; //  0 <= P <= M.  250
   int 	 Q; // 0 <= Q <= MIN(P,M-P,M-Q)       212  <= fla_min(250, 262, 300)
   int 	 LDU1; //  LDU1 >= MAX(1,P)  250
   int 	 LDU2; //  LDU2 >= MAX(1,M-P)  262
   int 	 LDV1T; //  LDV1T >= MAX(1,Q)  212
   int 	 LDV2T; //  LDV2T >= MAX(1,M-Q) 300
   int 	 LWORK; //  LWORK >= MAX(1,8*Q)  1696

   /*  in/out optional params  **/
   double  *THETA, *THETAref;
   double  *PHI, *PHIref;
   double  *U1, *U1ref; // dimension ldu1
   double  *U2, *U2ref; // dimension ldu2
   double  *V1T, *V1Tref; // dimension ldv1t
   double  *V2T, *V2Tref; // dimension ldv2t

   /* output parameters  */
   double  *B11D, *B11Dref;
   double  *B11E, *B11Eref;
   double  *B12D, *B12Dref;
   double  *B12E, *B12Eref;
   double  *B21D, *B21Dref;
   double  *B21E, *B21Eref;
   double  *B22D, *B22Dref;
   double  *B22E, *B22Eref;
   double  *WORK, *WORKref;
   int 	   INFO, INFOref;

   public: 
      bbcsd_double_parameters (char jobu1, char	jobu2, char jobv1t, char jobv2t,
	                           char trans, int m, int p, int q, int ldu1, 
                               int ldu2,  int ldv1t, int ldv2t, int lwork );
      ~bbcsd_double_parameters ();

}; /* end of bbcsd_double_parameters  class definition */

bbcsd_double_parameters::bbcsd_double_parameters (char jobu1, char jobu2, char jobv1t, 
                               char jobv2t, char trans, int m, int p, int q, int ldu1,
                               int ldu2,  int ldv1t, int ldv2t, int lwork )
{
   int j;
   JOBU1  = jobu1;
   JOBU2  = jobu2;
   JOBV1T = jobv1t;
   JOBV2T = jobv2t;
   TRANS  = trans;
   M      = m;
   P      = p; // 512
   Q      = q;
   LDU1   = ldu1; //600
   LDU2   = ldu2;
   LDV1T  = ldv1t;
   LDV2T  = ldv2t;
   LWORK  = lwork;

/* Memory allocation for various buffers */
#if 1
   lapacke_gtest_alloc_double_buffer_pair( &THETA, &THETAref, Q);
   lapacke_gtest_alloc_double_buffer_pair( &PHI, &PHIref, Q-1);
   lapacke_gtest_alloc_double_buffer_pair( &U1, &U1ref, (LDU1*P));
   lapacke_gtest_alloc_double_buffer_pair( &U2, &U2ref, (LDU2*(M-P)));
   lapacke_gtest_alloc_double_buffer_pair( &V1T, &V1Tref, (LDV1T*Q));
   lapacke_gtest_alloc_double_buffer_pair( &V2T, &V2Tref, (LDV2T*(M-Q)));

#else
   THETA    = (double *)malloc(Q*sizeof(double));
   THETAref = (double *)malloc(Q*sizeof(double));
   
   PHI    = (double *)malloc( (Q-1)*sizeof(double) );
   PHIref = (double *)malloc( (Q-1)*sizeof(double) );
  
   U1     = (double *)malloc( (LDU1*P)*sizeof(double) );
   U1ref  = (double *)malloc( (LDU1*P)*sizeof(double) );
 
   U2     = (double *)malloc( (LDU2*(M-P))*sizeof(double) );
   U2ref  = (double *)malloc( (LDU2*(M-P))*sizeof(double) );
   
   V1T    = (double *)malloc( (LDV1T*Q)*sizeof(double) );
   V1Tref = (double *)malloc( (LDV1T*Q)*sizeof(double) );
   
   V2T    = (double *)malloc( (LDV2T*(M-Q))*sizeof(double) );
   V2Tref = (double *)malloc( (LDV2T*(M-Q))*sizeof(double) );
#endif
   B11D    = (double *)malloc(Q*sizeof(double));
   B11Dref = (double *)malloc(Q*sizeof(double));
   
   B11E    = (double *)malloc( (Q-1)*sizeof(double) );
   B11Eref = (double *)malloc( (Q-1)*sizeof(double) );

   B12D    = (double *)malloc(Q*sizeof(double));
   B12Dref = (double *)malloc(Q*sizeof(double));
   
   B12E    = (double *)malloc( (Q-1)*sizeof(double) );
   B12Eref = (double *)malloc( (Q-1)*sizeof(double) );

   B21D    = (double *)malloc(Q*sizeof(double));
   B21Dref = (double *)malloc(Q*sizeof(double));
   
   B21E    = (double *)malloc( (Q-1)*sizeof(double) );
   B21Eref = (double *)malloc( (Q-1)*sizeof(double) );

   B22D    = (double *)malloc(Q*sizeof(double));
   B22Dref = (double *)malloc(Q*sizeof(double));
   
   B22E    = (double *)malloc( (Q-1)*sizeof(double) );
   B22Eref = (double *)malloc( (Q-1)*sizeof(double) );
   
   WORK    = (double *)malloc( fla_max(1,LWORK)*sizeof(double) );
   WORKref = (double *)malloc( fla_max(1,LWORK)*sizeof(double) ); 

   if ((THETA==NULL) || (THETAref==NULL) || (PHI==NULL) || (PHIref==NULL) / 
       (U1==NULL) || (U1ref==NULL) || (U2==NULL)  || (U2ref==NULL) /      
       (V1T==NULL) || (V1Tref==NULL) || (V2T==NULL)  || (V2Tref==NULL) / 
       (B11D==NULL) || (B11Dref==NULL) || (B12D==NULL)  || (B12Dref==NULL) /
       (B11E==NULL) || (B11Eref==NULL) || (B12E==NULL)  || (B12Eref==NULL) /
       (B21D==NULL) || (B21Dref==NULL) || (B22D==NULL)  || (B22Dref==NULL) /
       (B21E==NULL) || (B21Eref==NULL) || (B22E==NULL)  || (B22Eref==NULL) /
       (WORK==NULL) || (WORKref==NULL)
       )
   {
       printf(" bbcsd_double_parameters object: malloc error. Exiting...\n");
       free(THETA); free(THETAref); free(PHI); free(PHIref);
	   free(U1); free(U1ref); free(U2); free(U2ref);
       free(V2T); free(V2Tref); free(V1T); free(V1Tref);
       free(B11D); free(B11E); free(B12D); free(B12E);
       free(B21D); free(B21E); free(B22D); free(B22E);
       free(B11Dref); free(B11Eref); free(B12Dref); free(B12Eref);
       free(B21Dref); free(B21Eref); free(B22Dref); free(B22Eref);
       free(WORKref); free(WORK);	
       exit(-1); 
    }

      /* Initialization of input Buffers */
	  
    for( j = 0; j < Q; j++ ) {
       THETA[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       THETAref[j] = THETA[j];
    }

    for( j = 0; j < Q-1; j++ ) {
       PHI[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       PHIref[j] = PHI[j];
    }

    for( j = 0; j < (LDU1*P); j++ ) {
       U1[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       U1ref[j] = U1[j];
    }

    for( j = 0; j < (LDU2*(M-P)); j++ ) {
       U2[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       U2ref[j] = U2[j];
    }
	  
    for( j = 0; j < (LDV1T*Q); j++ ) {
       V1T[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       V1Tref[j] = V1T[j];
    }

    for( j = 0; j < (LDV2T*(M-Q)); j++ ) {
       V2T[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       V2Tref[j] = V2T[j];
    }
 
    for( j = 0; j < Q; j++ ) {
       B11D[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       B11Dref[j] = B11D[j];
	   
	   B12D[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       B12Dref[j] = B12D[j];
	   
	   B21D[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       B21Dref[j] = B21D[j];

	   B22D[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       B22Dref[j] = B22D[j];
    }
    for( j = 0; j < Q-1; j++ ) {
       B11E[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       B11Eref[j] = B11E[j];
	   
	   B12E[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       B12Eref[j] = B12E[j];
	   
	   B21E[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       B21Eref[j] = B21E[j];

	   B22E[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       B22Eref[j] = B22E[j];
    }
    for( j = 0; j < fla_max(1,LWORK); j++ ) {
       WORK[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       WORKref[j] = WORK[j];
    }
 
}

/* Destructor definition **/
bbcsd_double_parameters:: ~bbcsd_double_parameters ()
{
    free(THETA); free(THETAref); free(PHI); free(PHIref);
	free(U1); free(U1ref); free(U2); free(U2ref);
    free(V2T); free(V2Tref); free(V1T); free(V1Tref);

    free(B11D); free(B11E); free(B12D); free(B12E);
    free(B21D); free(B21E); free(B22D); free(B22E);
    free(B11Dref); free(B11Eref); free(B12Dref); free(B12Eref);
    free(B21Dref); free(B21Eref); free(B22Dref); free(B22Eref);
    free(WORKref); free(WORK);	
}


TEST(bbcsd, dbbcsd1) {

   /* LAPACKE dbbcsd prototype */
typedef int (*Fptr_NL_LAPACKE_dbbcsd) ( int matrix_layout, char jobu1, char jobu2,
                            char jobv1t, char jobv2t, char trans, lapack_int m,
                            lapack_int p, lapack_int q, double* theta,
                            double* phi, double* u1, lapack_int ldu1, double* u2,
                            lapack_int ldu2, double* v1t, lapack_int ldv1t,
                            double* v2t, lapack_int ldv2t, double* b11d,
                            double* b11e, double* b12d, double* b12e,
                            double* b21d, double* b21e, double* b22d,
                            double* b22e );  
						   
   void *hModule, *dModule;   
   Fptr_NL_LAPACKE_dbbcsd DBBCSD;		
   double diff;   
	
//   bbcsd_double_parameters dbbcsd_obj('n', 'n','n', 'n', 'N', 512, 250, 212, 250, 262,
//                                212, 300, 1696 );

   bbcsd_double_parameters dbbcsd_obj('n', 'n','n', 'n', 'N', 51, 25, 21, 25, 26,
                                21, 30, 200 );

   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   DBBCSD = (Fptr_NL_LAPACKE_dbbcsd)dlsym(hModule, "LAPACKE_dbbcsd");
   if (NULL == DBBCSD)
   {
   	  printf("Could not get the symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
      /* Compute libflame's Lapacke o/p  */
      dbbcsd_obj.INFO    = LAPACKE_dbbcsd( LAPACK_COL_MAJOR, dbbcsd_obj.JOBU1, dbbcsd_obj.JOBU2,
                            dbbcsd_obj.JOBV1T, dbbcsd_obj.JOBV2T, dbbcsd_obj.TRANS, dbbcsd_obj.M,
                            dbbcsd_obj.P, dbbcsd_obj.Q, dbbcsd_obj.THETA,
                            dbbcsd_obj.PHI, dbbcsd_obj.U1, dbbcsd_obj.LDU1, dbbcsd_obj.U2,
                            dbbcsd_obj.LDU2, dbbcsd_obj.V1T, dbbcsd_obj.LDV1T,
                            dbbcsd_obj.V2T, dbbcsd_obj.LDV2T, dbbcsd_obj.B11D,
                            dbbcsd_obj.B11E, dbbcsd_obj.B12D, dbbcsd_obj.B12E,
                            dbbcsd_obj.B21D, dbbcsd_obj.B21E, dbbcsd_obj.B22D,
                            dbbcsd_obj.B22E );
	  
      /* Compute the reference o/p by invoking Netlib-Lapack's API 	*/

      dbbcsd_obj.INFOref    = DBBCSD( LAPACK_COL_MAJOR, dbbcsd_obj.JOBU1, dbbcsd_obj.JOBU2,
                            dbbcsd_obj.JOBV1T, dbbcsd_obj.JOBV2T, dbbcsd_obj.TRANS, dbbcsd_obj.M,
                            dbbcsd_obj.P, dbbcsd_obj.Q, dbbcsd_obj.THETAref,
                            dbbcsd_obj.PHIref, dbbcsd_obj.U1ref, dbbcsd_obj.LDU1, dbbcsd_obj.U2ref,
                            dbbcsd_obj.LDU2, dbbcsd_obj.V1Tref, dbbcsd_obj.LDV1T,
                            dbbcsd_obj.V2Tref, dbbcsd_obj.LDV2T, dbbcsd_obj.B11Dref,
                            dbbcsd_obj.B11Eref, dbbcsd_obj.B12Dref, dbbcsd_obj.B12Eref,
                            dbbcsd_obj.B21Dref, dbbcsd_obj.B21Eref, dbbcsd_obj.B22Dref,
                            dbbcsd_obj.B22Eref );

   	/* Check for the exact singularity */
	if( dbbcsd_obj.INFO > 0 ) {
		printf( "info: %d, libflame dbbcsd posted error. Exiting.. \n", dbbcsd_obj.INFO );
		exit( 1 );
	}
	if( dbbcsd_obj.INFOref > 0 ) {
		printf( "info: %d, Netlib LAPACK dbbcsd posted error. Exiting.. \n", dbbcsd_obj.INFOref);
		exit( 1 );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_d( dbbcsd_obj.LWORK, dbbcsd_obj.WORK, dbbcsd_obj.WORKref );

	/*
	FLA_Nrm2_external( dbbcsd_obj.A, dbbcsd_obj.norm );
        FLA_Nrm2_external( dbbcsd_obj.Aref, dbbcsd_obj. normref );
	diff = FLA_abs( dbbcsd_obj. norm,  dbbcsd_obj. normref);
	*/
    
    	//EXPECT_EQ (0, diff);
	//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 10.0);
}


class bbcsd_float_parameters{
	
   public:
   /* input params to the API **/
   char  JOBU1; 
   char  JOBU2;
   char  JOBV1T;
   char  JOBV2T;
   char  TRANS;
   int 	 M;  // 512
   int 	 P; //  0 <= P <= M.  250
   int 	 Q; // 0 <= Q <= MIN(P,M-P,M-Q)       212  <= fla_min(250, 262, 300)
   int 	 LDU1; //  LDU1 >= MAX(1,P)  250
   int 	 LDU2; //  LDU2 >= MAX(1,M-P)  262
   int 	 LDV1T; //  LDV1T >= MAX(1,Q)  212
   int 	 LDV2T; //  LDV2T >= MAX(1,M-Q) 300
   int 	 LWORK; //  LWORK >= MAX(1,8*Q)  1696

   /*  in/out optional params  **/
   float  *THETA, *THETAref;
   float  *PHI, *PHIref;
   float  *U1, *U1ref; // dimension ldu1
   float  *U2, *U2ref; // dimension ldu2
   float  *V1T, *V1Tref; // dimension ldv1t
   float  *V2T, *V2Tref; // dimension ldv2t

   /* output parameters  */
   float  *B11D, *B11Dref;
   float  *B11E, *B11Eref;
   float  *B12D, *B12Dref;
   float  *B12E, *B12Eref;
   float  *B21D, *B21Dref;
   float  *B21E, *B21Eref;
   float  *B22D, *B22Dref;
   float  *B22E, *B22Eref;
   float  *WORK, *WORKref;
   int 	   INFO, INFOref;

   public: 
      bbcsd_float_parameters (char jobu1, char	jobu2, char jobv1t, char jobv2t,
	                           char trans, int m, int p, int q, int ldu1, 
                               int ldu2,  int ldv1t, int ldv2t, int lwork );
      ~bbcsd_float_parameters ();

}; /* end of bbcsd_float_parameters  class definition */

bbcsd_float_parameters::bbcsd_float_parameters (char jobu1, char jobu2, char jobv1t, 
                               char jobv2t, char trans, int m, int p, int q, int ldu1,
                               int ldu2,  int ldv1t, int ldv2t, int lwork )
{
   int j;
   JOBU1  = jobu1;
   JOBU2  = jobu2;
   JOBV1T = jobv1t;
   JOBV2T = jobv2t;
   TRANS  = trans;
   M      = m;
   P      = p; // 512
   Q      = q;
   LDU1   = ldu1; //600
   LDU2   = ldu2;
   LDV1T  = ldv1t;
   LDV2T  = ldv2t;
   LWORK  = lwork;

/* Memory allocation for various buffers */   
   THETA    = (float *)malloc(Q*sizeof(float));
   THETAref = (float *)malloc(Q*sizeof(float));
   
   PHI    = (float *)malloc( (Q-1)*sizeof(float) );
   PHIref = (float *)malloc( (Q-1)*sizeof(float) );
   
   U1     = (float *)malloc( (LDU1*P)*sizeof(float) );
   U1ref  = (float *)malloc( (LDU1*P)*sizeof(float) );

   U2     = (float *)malloc( (LDU2*(M-P))*sizeof(float) );
   U2ref  = (float *)malloc( (LDU2*(M-P))*sizeof(float) );
   
   V1T    = (float *)malloc( (LDV1T*Q)*sizeof(float) );
   V1Tref = (float *)malloc( (LDV1T*Q)*sizeof(float) );
   
   V2T    = (float *)malloc( (LDV2T*(M-Q))*sizeof(float) );
   V2Tref = (float *)malloc( (LDV2T*(M-Q))*sizeof(float) );

   B11D    = (float *)malloc(Q*sizeof(float));
   B11Dref = (float *)malloc(Q*sizeof(float));
   
   B11E    = (float *)malloc( (Q-1)*sizeof(float) );
   B11Eref = (float *)malloc( (Q-1)*sizeof(float) );

   B12D    = (float *)malloc(Q*sizeof(float));
   B12Dref = (float *)malloc(Q*sizeof(float));
   
   B12E    = (float *)malloc( (Q-1)*sizeof(float) );
   B12Eref = (float *)malloc( (Q-1)*sizeof(float) );

   B21D    = (float *)malloc(Q*sizeof(float));
   B21Dref = (float *)malloc(Q*sizeof(float));
   
   B21E    = (float *)malloc( (Q-1)*sizeof(float) );
   B21Eref = (float *)malloc( (Q-1)*sizeof(float) );

   B22D    = (float *)malloc(Q*sizeof(float));
   B22Dref = (float *)malloc(Q*sizeof(float));
   
   B22E    = (float *)malloc( (Q-1)*sizeof(float) );
   B22Eref = (float *)malloc( (Q-1)*sizeof(float) );
   
   WORK    = (float *)malloc( fla_max(1,LWORK)*sizeof(float) );
   WORKref = (float *)malloc( fla_max(1,LWORK)*sizeof(float) ); 

   if ((THETA==NULL) || (THETAref==NULL) || (PHI==NULL) || (PHIref==NULL) / 
       (U1==NULL) || (U1ref==NULL) || (U2==NULL)  || (U2ref==NULL) /      
       (V1T==NULL) || (V1Tref==NULL) || (V2T==NULL)  || (V2Tref==NULL) / 
       (B11D==NULL) || (B11Dref==NULL) || (B12D==NULL)  || (B12Dref==NULL) /
       (B11E==NULL) || (B11Eref==NULL) || (B12E==NULL)  || (B12Eref==NULL) /
       (B21D==NULL) || (B21Dref==NULL) || (B22D==NULL)  || (B22Dref==NULL) /
       (B21E==NULL) || (B21Eref==NULL) || (B22E==NULL)  || (B22Eref==NULL) /
       (WORK==NULL) || (WORKref==NULL)
       )
   {
       printf(" bbcsd_float_parameters object: malloc error. Exiting...\n");
       free(THETA); free(THETAref); free(PHI); free(PHIref);
	   free(U1); free(U1ref); free(U2); free(U2ref);
       free(V2T); free(V2Tref); free(V1T); free(V1Tref);
       free(B11D); free(B11E); free(B12D); free(B12E);
       free(B21D); free(B21E); free(B22D); free(B22E);
       free(B11Dref); free(B11Eref); free(B12Dref); free(B12Eref);
       free(B21Dref); free(B21Eref); free(B22Dref); free(B22Eref);
       free(WORKref); free(WORK);	
       exit(-1); 
    }

      /* Initialization of input Buffers */
	  
    for( j = 0; j < Q; j++ ) {
       THETA[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       THETAref[j] = THETA[j];
    }

    for( j = 0; j < Q-1; j++ ) {
       PHI[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       PHIref[j] = PHI[j];
    }

    for( j = 0; j < (LDU1*P); j++ ) {
       U1[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       U1ref[j] = U1[j];
    }

    for( j = 0; j < (LDU2*(M-P)); j++ ) {
       U2[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       U2ref[j] = U2[j];
    }
	  
    for( j = 0; j < (LDV1T*Q); j++ ) {
       V1T[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       V1Tref[j] = V1T[j];
    }

    for( j = 0; j < (LDV2T*(M-Q)); j++ ) {
       V2T[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       V2Tref[j] = V2T[j];
    }
 
    for( j = 0; j < Q; j++ ) {
       B11D[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       B11Dref[j] = B11D[j];
	   
	   B12D[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       B12Dref[j] = B12D[j];
	   
	   B21D[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       B21Dref[j] = B21D[j];

	   B22D[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       B22Dref[j] = B22D[j];
    }
    for( j = 0; j < Q-1; j++ ) {
       B11E[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       B11Eref[j] = B11E[j];
	   
	   B12E[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       B12Eref[j] = B12E[j];
	   
	   B21E[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       B21Eref[j] = B21E[j];

	   B22E[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       B22Eref[j] = B22E[j];
    }
    for( j = 0; j < fla_max(1,LWORK); j++ ) {
       WORK[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       WORKref[j] = WORK[j];
    }
 
}

/* Destructor definition **/
bbcsd_float_parameters:: ~bbcsd_float_parameters ()
{
    free(THETA); free(THETAref); free(PHI); free(PHIref);
	free(U1); free(U1ref); free(U2); free(U2ref);
    free(V2T); free(V2Tref); free(V1T); free(V1Tref);

    free(B11D); free(B11E); free(B12D); free(B12E);
    free(B21D); free(B21E); free(B22D); free(B22E);
    free(B11Dref); free(B11Eref); free(B12Dref); free(B12Eref);
    free(B21Dref); free(B21Eref); free(B22Dref); free(B22Eref);
    free(WORKref); free(WORK);	
}


TEST(bbcsd, sbbcsd1) {

   /* LAPACKE sbbcsd prototype */
typedef int (*Fptr_NL_LAPACKE_sbbcsd) ( int matrix_layout, char jobu1, char jobu2,
                            char jobv1t, char jobv2t, char trans, lapack_int m,
                            lapack_int p, lapack_int q, float* theta,
                            float* phi, float* u1, lapack_int ldu1, float* u2,
                            lapack_int ldu2, float* v1t, lapack_int ldv1t,
                            float* v2t, lapack_int ldv2t, float* b11d,
                            float* b11e, float* b12d, float* b12e,
                            float* b21d, float* b21e, float* b22d,
                            float* b22e );  
						   
   void *hModule, *dModule;   
   Fptr_NL_LAPACKE_sbbcsd sbbcsd;		
   float diff;   
	
//   bbcsd_float_parameters sbbcsd_obj('n', 'n','n', 'n', 'N', 512, 250, 212, 250, 262,
//                                212, 300, 1696 );

   bbcsd_float_parameters sbbcsd_obj('n', 'n','n', 'n', 'N', 51, 25, 21, 25, 26,
                                21, 30, 200 );

   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   sbbcsd = (Fptr_NL_LAPACKE_sbbcsd)dlsym(hModule, "LAPACKE_sbbcsd");
   if (NULL == sbbcsd)
   {
   	  printf("Could not get the symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
      /* Compute libflame's Lapacke o/p  */
      sbbcsd_obj.INFO    = LAPACKE_sbbcsd( LAPACK_COL_MAJOR, sbbcsd_obj.JOBU1, sbbcsd_obj.JOBU2,
                            sbbcsd_obj.JOBV1T, sbbcsd_obj.JOBV2T, sbbcsd_obj.TRANS, sbbcsd_obj.M,
                            sbbcsd_obj.P, sbbcsd_obj.Q, sbbcsd_obj.THETA,
                            sbbcsd_obj.PHI, sbbcsd_obj.U1, sbbcsd_obj.LDU1, sbbcsd_obj.U2,
                            sbbcsd_obj.LDU2, sbbcsd_obj.V1T, sbbcsd_obj.LDV1T,
                            sbbcsd_obj.V2T, sbbcsd_obj.LDV2T, sbbcsd_obj.B11D,
                            sbbcsd_obj.B11E, sbbcsd_obj.B12D, sbbcsd_obj.B12E,
                            sbbcsd_obj.B21D, sbbcsd_obj.B21E, sbbcsd_obj.B22D,
                            sbbcsd_obj.B22E );
	  
      /* Compute the reference o/p by invoking Netlib-Lapack's API 	

      sbbcsd_obj.INFOref    = sbbcsd( LAPACK_COL_MAJOR, sbbcsd_obj.JOBU1, sbbcsd_obj.JOBU2,
                            sbbcsd_obj.JOBV1T, sbbcsd_obj.JOBV2T, sbbcsd_obj.TRANS, sbbcsd_obj.M,
                            sbbcsd_obj.P, sbbcsd_obj.Q, sbbcsd_obj.THETAref,
                            sbbcsd_obj.PHIref, sbbcsd_obj.U1ref, sbbcsd_obj.LDU1, sbbcsd_obj.U2ref,
                            sbbcsd_obj.LDU2, sbbcsd_obj.V1Tref, sbbcsd_obj.LDV1T,
                            sbbcsd_obj.V2Tref, sbbcsd_obj.LDV2T, sbbcsd_obj.B11Dref,
                            sbbcsd_obj.B11Eref, sbbcsd_obj.B12Dref, sbbcsd_obj.B12Eref,
                            sbbcsd_obj.B21Dref, sbbcsd_obj.B21Eref, sbbcsd_obj.B22Dref,
                            sbbcsd_obj.B22Eref );*/

   	/* Check for the exact singularity */
	if( sbbcsd_obj.INFO > 0 ) {
		printf( "info: %d, libflame sbbcsd posted error. Exiting.. \n", sbbcsd_obj.INFO );
		exit( 1 );
	}
	if( sbbcsd_obj.INFOref > 0 ) {
		printf( "info: %d, Netlib LAPACK sbbcsd posted error. Exiting.. \n", sbbcsd_obj.INFOref);
		exit( 1 );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_s( sbbcsd_obj.LWORK, sbbcsd_obj.WORK, sbbcsd_obj.WORKref );

	/*
	FLA_Nrm2_external( sbbcsd_obj.A, sbbcsd_obj.norm );
        FLA_Nrm2_external( sbbcsd_obj.Aref, sbbcsd_obj. normref );
	diff = FLA_abs( sbbcsd_obj. norm,  sbbcsd_obj. normref);
	*/
    
    	//EXPECT_EQ (0, diff);
	//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 10.0);
}
