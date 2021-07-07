#include <cstdio>
#include <iostream>
#include <string>
#include <cstring>
#include "benchmark/benchmark.h"

#include "gtest/gtest.h"
//#include "gtest/gtest.h"
#include "Flame_gbench_main.h"
//#include "Flame_gbench_aux.h"
#include "lapack.h"

#include <stdio.h>
using namespace std;

/* Global variables */
Lin_solver_paramlist lin_solver_paramslist[30];
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
   int num_tests;
   int mode;
   int num_ranges;
   
   fp = fopen( file_name, "r");
   if (fp == NULL){
    printf("Error: Lin solver config file missing. Exiting.. \n");
    exit(-1);
   }
   
   /* Read the mode */
   fscanf(fp, "%s", &line[0]);
   fscanf(fp, "%d", &mode);
   fscanf(fp, "%*[^\n]\n");
   
   /* Read the number of Ranges */
   fscanf(fp, "%s", &line[0]);
   fscanf(fp, "%d", &num_ranges);

   fscanf(fp, "%s", &line[0]); // Range_start

   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].m_range_start);
	  lin_solver_paramslist[i].num_ranges = num_ranges;
	  lin_solver_paramslist[i].mode = mode;
   }

   fscanf(fp, "%s", &line[0]); // Range_end
   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].m_range_end);
   }

   fscanf(fp, "%s", &line[0]); // Range_step_size
   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].m_range_step_size);
   }

   fscanf(fp, "%s", &line[0]); // Range_start

   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].n_range_start);
   }

   fscanf(fp, "%s", &line[0]); // Range_end
   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].n_range_end);
   }

   fscanf(fp, "%s", &line[0]); // Range_step_size
   for (i=0; i<num_ranges; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].n_range_step_size);
   }

   fscanf(fp, "%s", &line[0]);
   fscanf(fp, "%d", &num_tests);
   for (i=0; i<num_tests; i++){
      lin_solver_paramslist[i].num_tests = num_tests;   
   }
   
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].matrix_layout);
   }

   str = &line[0];
   fscanf(fp, "%s", str);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      lin_solver_paramslist[i].transr = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      lin_solver_paramslist[i].Uplo = *str;
   }
   
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].m);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].n);
   }
  
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].nrhs);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].lda);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].ldb);
   }
  
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].ldab);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].kl);
   }
  
   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].ku);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].kd);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      lin_solver_paramslist[i].diag = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      lin_solver_paramslist[i].fact = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      lin_solver_paramslist[i].equed = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      lin_solver_paramslist[i].symm = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%f", &lin_solver_paramslist[i].solver_threhold);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      lin_solver_paramslist[i].equed_porfsx = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].n_err_bnds_porfsx);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].nparams_porfsx);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%s", str);
      lin_solver_paramslist[i].norm_gbcon = *str;
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].kl_gbcon);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
      fscanf(fp, "%d", &lin_solver_paramslist[i].ku_gbcon);
   }

   fscanf(fp, "%s", &line[0]);
   for (i=0; i<num_tests; i++){
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

void Read_config_parameters() 
{
  /*Read Linear API parameters from config file */
  Read_Lin_solver_param("Config/LIN_SLVR.dat");
  Read_Lin_driver_param("Config/LIN_DRVR.dat");
  
  /*Read eigen parameters from config file */
  Read_EIG_params("Config/EIG_PARAMS.dat");
  Read_EIG_non_sym_params("Config/EIG_NSYM_PARAMS.dat");
  Read_SVD_param ( "Config/SVD.dat" );
}

class MyInit
{
public:
    MyInit()
    {
		Read_config_parameters ();
	}
};

MyInit my_init;
// Invoking Benchmark main
BENCHMARK_MAIN();
