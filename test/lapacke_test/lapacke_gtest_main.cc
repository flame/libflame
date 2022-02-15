#include "gtest/gtest.h"
#include "lapacke_gtest_main.h"
#include "lapacke_gtest_aux.h"
#include "lapack.h"

#include <stdio.h>

/* Macros */
#define common_parameters_free() \
       free (A   ); \
       free (Aref); \
       free (ipiv   ); \
       free (ipivref); \
       free (b   ); \
       free (bref);

/* Global variables */
#define MAC2NAME(s) NAME(s)
#define NAME(s) #s
/*
The macros 'REF_LPKE_LIB', 'REF_BLAS_LIB' hold the paths of 
'Netlib refrence libraries (liblapacke.so, libblas.so)' respectively.
The above macros are provided in the build system
*/
const char *NETLIB_LAPACKE_LIB = MAC2NAME(REF_LPKE_LIB);
const char *NETLIB_BLAS_LIB = MAC2NAME(REF_BLAS_LIB);

Lin_solver_paramlist lin_solver_paramslist[NUM_SUB_TESTS];
Lin_driver_paramlist lin_driver_paramslist[NUM_SUB_TESTS];
EIG_paramlist eig_paramslist[NUM_SUB_TESTS];
EIG_Non_symmetric_paramlist_t eig_non_sym_paramslist[NUM_SUB_TESTS];
SVD_paramlist svd_paramslist[NUM_SUB_TESTS];


/* Functions definitions  */

/* This function reads parameters needed for Linear solver APIs 
   from the config settings file 'LIN_SLVR.dat' and saves in the 
   'lin_solver_paramslist' structure array   */
void Read_Lin_solver_param ( const char *file_name )
{
   FILE *fp;
   int i;
   char line[20];
   char *str;
   
   fp = fopen( file_name, "r");
   if (fp == NULL){
    printf("Error: Lin solver config file missing. Exiting.. \n");
    exit(-1);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].matrix_layout);
   }

   str = &line[0];
   fscanf(fp, "%s", str);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      lin_solver_paramslist[i].transr = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      lin_solver_paramslist[i].Uplo = *str;
   }
   
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].m);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].n);
   }
  
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].nrhs);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].lda);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].ldb);
   }
  
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].ldab);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].kl);
   }
  
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].ku);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].kd);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      lin_solver_paramslist[i].diag = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      lin_solver_paramslist[i].fact = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      lin_solver_paramslist[i].equed = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      lin_solver_paramslist[i].symm = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%f", &lin_solver_paramslist[i].solver_threhold);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      lin_solver_paramslist[i].equed_porfsx = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].n_err_bnds_porfsx);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].nparams_porfsx);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      lin_solver_paramslist[i].norm_gbcon = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].kl_gbcon);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].ku_gbcon);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].ldab_gbcon);
   }

   fclose(fp);
   
}

/* This function reads additional params needed for Linear solver driver 
   APIs from the config settings file 'LIN_DRVR.dat' and saves in the 
   'lin_driver_paramslist' structure array   */
void Read_Lin_driver_param ( const char *file_name )
{
   FILE *fp;
   int i;
   char line[20];
   char *str;
   fp = fopen( file_name, "r");
   str = &line[0];

   if (fp == NULL){
    printf("Error: Lin solver config file missing. Exiting.. \n");
    exit(-1);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &lin_driver_paramslist[i].nparams);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &lin_driver_paramslist[i].nerrbnds);
   }
   
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      lin_driver_paramslist[i].fact = *str;
   }
   
   fclose(fp);

}

/* This function reads parameters needed for Eigen APIs 
   from the config settings file 'EIG_PARAMS.dat' and saves in the 
   'eig_paramslist' structure array   */
void Read_EIG_params( const char *file_name )
{
   FILE *fp;
   int i;
   char line[20];
   char *str;
   
   fp = fopen( file_name, "r");
   if (fp == NULL){
    printf("Error: EIG params config file missing. Exiting.. \n");
    exit(-1);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &eig_paramslist[i].matrix_layout);
   }

   str = &line[0];
   fscanf(fp, "%s", str);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_paramslist[i].trans = *str;
   }
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_paramslist[i].uplo = *str;
   }
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_paramslist[i].job = *str;
   }
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_paramslist[i].vect = *str;
   }
   
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &eig_paramslist[i].m);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &eig_paramslist[i].n);
   }
   
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &eig_paramslist[i].p);
   }
  
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &eig_paramslist[i].nrhs);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &eig_paramslist[i].lda);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &eig_paramslist[i].ldb);
   }
  
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &eig_paramslist[i].nb);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &eig_paramslist[i].ldt);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &eig_paramslist[i].k);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &eig_paramslist[i].isgn);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_paramslist[i].compz = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &eig_paramslist[i].kb);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &eig_paramslist[i].itype);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_paramslist[i].vect_rd = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_paramslist[i].side = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_paramslist[i].job_seqr = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_paramslist[i].eigsrc = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_paramslist[i].initv = *str;
   }
   
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_paramslist[i].norm = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_paramslist[i].diag = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_paramslist[i].storev = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &eig_paramslist[i].tsize);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &eig_paramslist[i].threshold_value);
   }


   fclose(fp);
   
}

/* This function reads parameters needed for Non symmetric Eigen APIs 
   from the config settings file 'EIG_NSYM_PARAMS.dat' and saves in the 
   'eig_paramslist' structure array   */
void Read_EIG_non_sym_params( const char *file_name )
{
   FILE *fp;
   int i;
   char line[20];
   char *str;
   
   str = &line[0];
   fp = fopen( file_name, "r");
   if (fp == NULL){
    printf("Error: EIG non symmetric API params config file missing. Exiting.. \n");
    exit(-1);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].howmny = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].initv = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].job_seqr = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].eigsrc = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].initv = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].job = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].howmny_trsna = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].job_trsen = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].compq = *str;
   }

   /* Reading config params for 'trsyl' API  */
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].trana_real = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].trana_complex = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].tranb_real= *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].tranb_complex = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &eig_non_sym_paramslist[i].isgn);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%f", &eig_non_sym_paramslist[i].gghrd_threshold);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%f", &eig_non_sym_paramslist[i].ggbal_threshold);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%f", &eig_non_sym_paramslist[i].GenNonSymEigProblem_threshold);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].compq_hgeqz = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].compz_hgeqz = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].side_tgevc = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].jobvsl = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].jobvsr = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].sort = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].sense_ggesx = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].balance_ggevx = *str;
   }
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].sense_ggevx = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].sort_gees = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &eig_non_sym_paramslist[i].wantz);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &eig_non_sym_paramslist[i].wantq);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &eig_non_sym_paramslist[i].tgsen_ijob);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
	  fscanf(fp, "%s", str);
      eig_non_sym_paramslist[i].unmhr_trans = *str;
   }


   fclose(fp);
}

/* This function reads parameters needed for SVD APIs 
   from the config settings file 'SVD.dat' and saves in the 
   'svd_paramslist' structure array   */
void Read_SVD_param ( const char *file_name )
{
   FILE *fp;
   int i;
   char line[25];
   char *str;
   
   str = &line[0];
   fp = fopen( file_name, "r");
   if (fp == NULL){
    printf("Error: SVD config file missing. Exiting.. \n");
    exit(-1);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &svd_paramslist[i].matrix_layout);
   }

   str = &line[0];
   fscanf(fp, "%s", str);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      svd_paramslist[i].jobu = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      svd_paramslist[i].jobv = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      svd_paramslist[i].jobq = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &svd_paramslist[i].m);
   }
  
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &svd_paramslist[i].p);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &svd_paramslist[i].n);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%f", &svd_paramslist[i].tola);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%f", &svd_paramslist[i].tolb);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%f", &svd_paramslist[i].svd_threshold);
   }

   str = &line[0];
   fscanf(fp, "%s", str);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      svd_paramslist[i].jobu_gesvd = *str;
   }

   str = &line[0];
   fscanf(fp, "%s", str);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      svd_paramslist[i].jobvt_gesvd = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      svd_paramslist[i].joba_gejsv = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      svd_paramslist[i].jobu_gejsv = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      svd_paramslist[i].jobv_gejsv = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      svd_paramslist[i].jobr_gejsv = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      svd_paramslist[i].jobt_gejsv = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      svd_paramslist[i].jobp_gejsv = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &svd_paramslist[i].m_gejsv);
   }
  
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &svd_paramslist[i].n_gejsv);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      svd_paramslist[i].joba_gesvj = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      svd_paramslist[i].jobu_gesvj = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      svd_paramslist[i].jobv_gesvj = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &svd_paramslist[i].m_gesvj);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &svd_paramslist[i].n_gesvj);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &svd_paramslist[i].mv_gesvj);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%f", &svd_paramslist[i].ctol_gesvj);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", &svd_paramslist[i].jobu_gesvdx);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", &svd_paramslist[i].jobvt_gesvdx);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", &svd_paramslist[i].range_gesvdx);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &svd_paramslist[i].il);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%d", &svd_paramslist[i].iu);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%f", &svd_paramslist[i].vl);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%f", &svd_paramslist[i].vu);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      svd_paramslist[i].joba_gesvdq = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      svd_paramslist[i].jobu_gesvdq = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<NUM_SUB_TESTS; i++){
      fscanf(fp, "%s", str);
      svd_paramslist[i].jobv_gesvdq= *str;
   }

   fclose(fp);
}


int main(int argc, char **argv) 
{
  /*Read Linear API parameters from config file */
  Read_Lin_solver_param("Config/LIN_SLVR.dat");
  Read_Lin_driver_param("Config/LIN_DRVR.dat");
  
  /*Read eigen parameters from config file */
  Read_EIG_params("Config/EIG_PARAMS.dat");
  Read_EIG_non_sym_params("Config/EIG_NSYM_PARAMS.dat");
  Read_SVD_param ( "Config/SVD.dat" );
  
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

/* Constructor definition  double_common_parameters */
double_common_parameters:: double_common_parameters ( int nrow, int ncol)
   {
      int i, j;
      m = nrow;
      n = ncol; // set test matrix size
      nrhs = 1;
      lda = m;
      ldb = n;
      
      hModule = NULL;
      dModule = NULL;

      A = (double *)malloc(m*n*sizeof(double)) ;
      b = (double *)malloc(m*nrhs*sizeof(double)) ;
      ipiv = (lapack_int *)malloc(m*sizeof(lapack_int));        
      if ((ipiv==NULL) || (A==NULL) || (b==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        free(A); free(b); free(ipiv);
        exit(0); 
      }
      
      /* Allocation of memory for capturing reference o/ps */
      Aref = (double *)malloc(m*n*sizeof(double));
      bref = (double *)malloc(m*nrhs*sizeof(double)) ;
      ipivref = (lapack_int *)malloc(m*sizeof(lapack_int)) ;
      if ((ipivref==NULL) || (Aref==NULL) || (bref==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        free(Aref); free(bref); free(ipivref);
        exit(0); 
      }

      /* Initialization of input matrices */
      for( i = 0; i < m; i++ ) {
         for( j = 0; j < n; j++ ) {
           A[i+j*lda] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
           Aref[i+j*lda] = A[i+j*lda];
         }
      }
      for(i=0;i<m*nrhs;i++) {
         b[i] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
         bref[i] = b[i];
      }
          
   } /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
double_common_parameters:: ~double_common_parameters ()
{
   /* De-Allocate memory for the input matrices */  
   if (A!=NULL){
      free(A); 
   }
   if (Aref!=NULL){
      free(Aref); 
   }
   if (b!=NULL){
      free(b); 
   }
   if (bref!=NULL){
      free(bref); 
   }
   if (ipiv!=NULL){
      free(ipiv); 
   }
   if (ipivref!=NULL){
      free(ipivref); 
   }

}

/* Constructor definition  float_common_parameters */
float_common_parameters:: float_common_parameters ( int nrow, int ncol)
   {
      int i, j;
      m = nrow;
      n = ncol; // set test matrix size
      nrhs = 1;
      lda = m;
      ldb = n;
      
      hModule = NULL;
      dModule = NULL;

      A = (float *)malloc(m*n*sizeof(float)) ;
      b = (float *)malloc(m*nrhs*sizeof(float)) ;
      ipiv = (lapack_int *)malloc(m*sizeof(lapack_int));        
      if ((ipiv==NULL) || (A==NULL) || (b==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        free(A); free(b); free(ipiv);
        exit(0); 
      }
      
      /* Allocation of memory for capturing reference o/ps */
      Aref = (float *)malloc(m*n*sizeof(float)) ;
      bref = (float *)malloc(m*nrhs*sizeof(float)) ;
      ipivref = (lapack_int *)malloc(m*sizeof(lapack_int)) ;
      if ((ipivref==NULL) || (Aref==NULL) || (bref==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        free(Aref); free(bref); free(ipivref);
        exit(0); 
      }

      /* Initialization of input matrices */
      for( i = 0; i < m; i++ ) {
         for( j = 0; j < n; j++ ) {
           A[i+j*lda] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
           Aref[i+j*lda] = A[i+j*lda];
         }
      }
      for(i=0;i<m*nrhs;i++) {
         b[i] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
         bref[i] = b[i];
      }
          
   } /* end of Constructor  */


/* Destructor definition  'float_common_parameters' class  */
float_common_parameters:: ~float_common_parameters ()
{
   /* De-Allocate memory for the input matrices */  
   if (A!=NULL){
      free(A); 
   }
   if (Aref!=NULL){
      free(Aref); 
   }
   if (b!=NULL){
      free(b); 
   }
   if (bref!=NULL){
      free(bref); 
   }
   if (ipiv!=NULL){
      free(ipiv); 
   }
   if (ipivref!=NULL){
      free(ipivref); 
   }

}

/* Constructor definition  scomplex_common_parameters */
scomplex_common_parameters:: scomplex_common_parameters ( int nrow, int ncol)
{
   m = nrow;
   n = ncol; // set test matrix size
   nrhs = 1;
   lda = m;
   ldb = n;
 
   hModule = NULL;
   dModule = NULL;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &A,  &Aref,  (m*n));
   lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b,  &bref,  (m*nrhs));
   lapacke_gtest_alloc_int_buffer_pair( &ipiv,  &ipivref,  (m*sizeof(lapack_int)));

   if ((ipiv==NULL) || (A==NULL) || (b==NULL) || \
       (ipivref==NULL) || (Aref==NULL) || (bref==NULL)) {
      printf("error of memory allocation. Exiting ...\n");
      free(A); free(b); free(ipiv);
      free(Aref); free(bref); free(ipivref);
      exit(-1); 
   }
    
    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, (m*n));
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, (m*nrhs));

} /* end of Constructor  */

scomplex_common_parameters:: ~scomplex_common_parameters ()
{
   #if LAPACKE_TEST_VERBOSE
   printf(" scomplex_common_parameters object: destructor invoked. \n");
   #endif

   common_parameters_free();
}

/* Constructor definition  dcomplex_common_parameters */
dcomplex_common_parameters:: dcomplex_common_parameters ( int nrow, int ncol)
{
   m = nrow;
   n = ncol; // set test matrix size
   nrhs = 1;
   lda = m;
   ldb = n;
 
   hModule = NULL;
   dModule = NULL;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &A,  &Aref,  (m*n));
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b,  &bref,  (m*nrhs));
   lapacke_gtest_alloc_int_buffer_pair( &ipiv,  &ipivref,  (m*sizeof(lapack_int)));

   if ((ipiv==NULL) || (A==NULL) || (b==NULL) || \
       (ipivref==NULL) || (Aref==NULL) || (bref==NULL)) {
      printf("error of memory allocation. Exiting ...\n");
      free(A); free(b); free(ipiv);
      free(Aref); free(bref); free(ipivref);
      exit(-1); 
   }
    
    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, (m*n));
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, (m*nrhs));

} /* end of Constructor  */

dcomplex_common_parameters:: ~dcomplex_common_parameters ()
{
   #if LAPACKE_TEST_VERBOSE
   printf(" dcomplex_common_parameters object: destructor invoked. \n");
   #endif

   common_parameters_free();
}

