/******************************************************************************
* Copyright (C) 2021-2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file main.cc
 *  @brief Defines main function and other common functions to use in APIs
           of CPP template interface.
 *  */

#include <stdio.h>
#include <gtest/gtest.h> 
#include "main.h"

char NETLIB_LAPACK_LIB[60] = NETLIB_LAPACK_LIB_PATH;
char NETLIB_BLAS_LIB[60] = NETLIB_BLAS_LIB_PATH;

void *blasModule = NULL, *lapackModule = NULL;

// Global arrays to hold test parameters data.
EIG_paramlist eig_paramslist[NUM_SUB_TESTS];
Lin_solver_paramlist lin_solver_paramslist[NUM_SUB_TESTS];
Lin_driver_paramlist lin_driver_paramslist[NUM_SUB_TESTS];

/*! @brief  Read_EIG_params is function used to read and initialize 
      Symmetric Eigen values/vectors routines input data into global array.
 * @details
 * \b Purpose:
    \verbatim
    This function reads parameters needed for Eigen APIs 
    from the config settings file 'EIG_PARAMS.dat' and saves in the 
    'eig_paramslist' structure array
    \endverbatim
	
 * @param[in] file_name
    file_name is charater array.
    Used to specify the name of the file to read the data.

 * @return void
    Nothing.
 * */
void Read_EIG_params(const char *file_name)
{
  FILE *fp;
  int index;
  char line[30];
  char *str = &line[0];

  fp = fopen(file_name, "r");
  if (fp == NULL) {
  printf("Error: EIG params config file missing. Exiting.. \n");
  exit(-1);
  }

  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%s", str);
    eig_paramslist[index].jobz = *str;
  }

  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%s", str);
    eig_paramslist[index].uplo = *str;
  }

  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].n);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].subda);
  }

  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].sda);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].ldab);
  }

  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].ldz);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%s", str);
    eig_paramslist[index].jobz_2stage = *str;
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lwork);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lwork_hbevd);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lrwork_hbevd);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].liwork_hbevd);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%s", str);
    eig_paramslist[index].range = *str;
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].ldq);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%lf", &eig_paramslist[index].vl);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%lf", &eig_paramslist[index].vu);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].il);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].iu);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%lf", &eig_paramslist[index].abstol);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].ldz_hbgv);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].sdb);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].subdb);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].ldab_hbgv);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].ldbb);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lwork_hbgvd);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].ldq_hbgvx);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%s", str);
    eig_paramslist[index].vect = *str;
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%s", str);
    eig_paramslist[index].vect_hbgst = *str;
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].ldx);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lda);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lwork_heev_2stage);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lwork_heevd_2stage);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lwork_heevr);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lrwork_heevr);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].liwork_heevr);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].itype);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].ldb);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lwork_hegv);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].i1);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].i2);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lwork_hetrd);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lhous2);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lwork_hetrd_2stage);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%s", str);
    eig_paramslist[index].vect_hetrd_2stage = *str;
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%s", str);
    eig_paramslist[index].stage1 = *str;
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].kd);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lwork_hetrd_hb2st);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lhous);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lwork_hbtrd_he2hb);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%s", str);
    eig_paramslist[index].norm = *str;
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lwork_lansy);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].m);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lwork_lange);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lda_lange);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lwork_geqp3);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lwork_geqrf);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].k);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lwork_syev);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].lwork_syevd);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &eig_paramslist[index].liwork_syevd);
  }
  
  fclose(fp);
}

/*! @brief  Read_Lin_solver_params is function used to read and initialize 
      Linear solver routines input data into global array.
 * @details
 * \b Purpose:
    \verbatim
      This function reads parameters needed for Linear solver APIs 
     from the config settings file 'LIN_SLVR.dat' and saves in the 
     'lin_solver_paramslist' structure array.
    \endverbatim
	
 * @param[in] file_name
    file_name is charater array.
    Used to specify the name of the file to read the data.

 * @return void
    Nothing.
 * */
void Read_Lin_solver_params (const char *file_name)
{
  FILE *fp;
  int index;
  char line[20];
  char *str;

  fp = fopen(file_name, "r");
  if (fp == NULL){
    printf("Error: Lin solver config file missing. Exiting.. \n");
    exit(-1);
  }

  str = &line[0];
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%s", str);
    lin_solver_paramslist[index].uplo = *str;
  }

  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_solver_paramslist[index].n);
  }

  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_solver_paramslist[index].lda);
  }

  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%lf", &lin_solver_paramslist[index].anorm);
  }

  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_solver_paramslist[index].itype);
  }

  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_solver_paramslist[index].ldb);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_solver_paramslist[index].nrhs);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_solver_paramslist[index].ldaf);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_solver_paramslist[index].ldx);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_solver_paramslist[index].m);
  }

  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_solver_paramslist[index].lda_getrf);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%s", str);
    lin_solver_paramslist[index].trans = *str;
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%s", str);
    lin_solver_paramslist[index].norm = *str;
  }

  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_solver_paramslist[index].lwork);
  }
  
  fclose(fp);
}

/*! @brief  Read_Lin_driver_params is function used to read and initialize 
      Linear driver routines input data into global array.
 * @details
 * \b Purpose:
    \verbatim
     This function reads parameters needed for Linear driver APIs 
     from the config settings file 'LIN_DRVR.dat' and saves in the 
     'lin_driver_paramslist' structure array.
    \endverbatim
	
 * @param[in] file_name
    file_name is charater array.
    Used to specify the name of the file to read the data.

 * @return void
    Nothing.
 * */
void Read_Lin_driver_params (const char *file_name)
{
  FILE *fp;
  int index;
  char line[20];
  char *str;

  fp = fopen(file_name, "r");
  if (fp == NULL){
    printf("Error: Lin driver config file missing. Exiting.. \n");
    exit(-1);
  }

  str = &line[0];
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%s", str);
    lin_driver_paramslist[index].uplo = *str;
  }

  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_driver_paramslist[index].n);
  }

  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_driver_paramslist[index].nrhs);
  }

  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_driver_paramslist[index].lda);
  }

  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_driver_paramslist[index].ldb);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_driver_paramslist[index].lwork);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_driver_paramslist[index].lwork_hesv_aa);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_driver_paramslist[index].ltb);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_driver_paramslist[index].lwork_hesv_aa_2stage);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%s", str);
    lin_driver_paramslist[index].fact = *str;
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_driver_paramslist[index].ldaf);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_driver_paramslist[index].ldx);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%s", str);
    lin_driver_paramslist[index].fact_hesvxx = *str;
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%s", str);
    lin_driver_paramslist[index].equed = *str;
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_driver_paramslist[index].n_err_bnds);
  }
  
  fscanf(fp, "%s", &line[0]);
  for (index = 0; index < NUM_SUB_TESTS; index++) {
    fscanf(fp, "%d", &lin_driver_paramslist[index].nparams);
  }
  
  fclose(fp);
}

/*! @brief  closelibs is function used to close the blasModule and lapackModule
            dynamic library.
 * @details
 * \b Purpose:
    \verbatim
    This function is used to close the blasModule and lapackModule
    dynamic library.
    \endverbatim
	
 * @param void
    Nothing.

 * @return void
    Nothing.
 * */
void closelibs(void) {
  dlclose(blasModule);
  dlclose(lapackModule);
}

/*! @brief  main is function to initialize all input data into global arrays.
            And calls google test APIs to perform the testing and returns
            the result.
 * @details
 * \b Purpose:
    \verbatim
    This function to initialize all input data into global arrays.
    And calls google test APIs to perform the testing and returns
    the result.
    \endverbatim
	
 * @param[in] argc
    Count of the arguments passed along with main while execution.
 * @param[in] argv
    The arguments passed along with main while execution.

 * @return int
    Returns the result as integer.
 * */
int main(int argc, char **argv) 
{
  // Read eigen parameters from config file.
  Read_EIG_params("config/EIG_PARAMS.dat");
  
  // Read linear solver parameters from config file.
  Read_Lin_solver_params("config/LIN_SLVR.dat");
  
  // Read linear driver parameters from config file.
  Read_Lin_driver_params("config/LIN_DRVR.dat");
  
  // Check for opening BLAS & LAPACK library files.
  blasModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
  lapackModule = dlopen(NETLIB_LAPACK_LIB, RTLD_NOW);
  if ((blasModule == NULL) || (lapackModule == NULL)) {
    printf("Open BLAS & LAPACK Library failed. Exiting ....\n");
    exit (-1);
  } else {
	  printf("Open BLAS & LAPACK Library successful.\n");
  }

  // Runs all TEST() in test API files.
  testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();
  
  // Close the opened library files and return result of tests.
  closelibs();
  return result;
}
