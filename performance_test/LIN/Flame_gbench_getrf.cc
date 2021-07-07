// Standard headers
#include <cstdio>
#include <iostream>
#include <string>
#include <cstring>
#include <cassert>
#include <memory>
#include <math.h>
#include <stdio.h>
#include "benchmark/benchmark.h"
#include "gtest/gtest.h"

// Library headers
#include "../Flame_gbench_main.h"
#include "lapack.h"

using namespace std;

/* Macros */
#define getrf_parameters_free() \
   if (A!=NULL){ \
      free(A); \
   }\
   if (b!=NULL){ \
      free(b);  \
   } \
   if (ipiv!=NULL){ \
      free(ipiv);  \
   } 

class getrf_double_parameters{
   public:  
      int m, n, nrhs, lda, ldb, info, inforef;
      void *hModule, *dModule;

      /* Local arrays */
      double *A, *b;
      int *ipiv, *ipivref;
      double norm, normref; 
   public: 
      getrf_double_parameters ( int num_row, int num_col, int lda_, int ldb_);
      ~getrf_double_parameters ();   
    
};  /* end of getrf_double_parameters  class definition */



/* Constructor definition  getrf_double_parameters */
getrf_double_parameters:: getrf_double_parameters ( int nrow, int ncol, int lda_, int ldb_)
   {
      int i, j;
      m = nrow;
      n = ncol; // set test matrix size
      nrhs = 1;
      lda = lda_;
      ldb = ldb_;
      
      A = (double *)malloc(m*n*sizeof(double)) ;
      b = (double *)malloc(m*nrhs*sizeof(double)) ;
      ipiv = (int *)malloc(m*sizeof(int));        
      if ((ipiv==NULL) || (A==NULL) || (b==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        free(A); free(b); free(ipiv);
        exit(0); 
      }
      
      /* Initialization of input matrices */
      for( i = 0; i < m; i++ ) {
         for( j = 0; j < n; j++ ) {
           A[i+j*lda] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
         }
      }
      for(i=0;i<m*nrhs;i++) {
         b[i] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
      }
          
   } /* end of Constructor  */

/* Destructor definition  'getrf_double_parameters' class  */
getrf_double_parameters:: ~getrf_double_parameters ()
{
   //printf("\n buffer free happening \n");
   getrf_parameters_free ();
}


class dgetrf_test : public ::benchmark::Fixture {
 public:
   int m,n, lda, ldb;
   std::unique_ptr<getrf_double_parameters> data;

 void SetUp(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
    m = static_cast<size_t>(state.range(0));
    n = static_cast<size_t>(state.range(1));
	
    assert(data.get() == nullptr);
    data.reset(new getrf_double_parameters(m,n,m,n));

  }

  void TearDown(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      assert(data.get() != nullptr);
      data.reset();
  }

  ~dgetrf_test() { assert(data == nullptr); }

};

BENCHMARK_DEFINE_F(dgetrf_test, dgetrf1)(benchmark::State &state) {
  assert(data.get() != nullptr);
  int forced_loop_count = 1;
  
  if ( (data->m > FLAME_BIG_MATRIX_SIZE) || (data->n > FLAME_BIG_MATRIX_SIZE )){
	  forced_loop_count = FLAME_GBENCH_FORCED_ITERATION_COUNT;
  }
  
  //printf(" data->m: %d ,data->n : %d lda: %d, ldb: %d\n", data->m,data->n, data->lda, data->ldb);

  for (auto _ : state) {
	for (int i=0; i<forced_loop_count; i++) {
		/* Invoke the getrf API for benchmarking */  
		dgetrf_(&data->m, &data->n,  data->A, &data->lda, data->ipiv, &data->info);
		//state.PauseTiming();

		if( data->info < 0 ) {
			printf( "\n warning: The %d th argument dgetrf is wrong\n", data->info );
		}
		//state.ResumeTiming();
	}
  }
}

static void CustomArguments(benchmark::internal::Benchmark* b) {
	if (lin_solver_paramslist[0].mode == MODE_RANGE) {
	   for (int i = 0; i < lin_solver_paramslist[0].num_tests ; i++){
	       b->ArgsProduct({ 
              {benchmark::CreateDenseRange( lin_solver_paramslist[i].m_range_start, 
											lin_solver_paramslist[i].m_range_end, 
											lin_solver_paramslist[i].m_range_step_size)} , 

			  {benchmark::CreateDenseRange( lin_solver_paramslist[i].n_range_start, 
											lin_solver_paramslist[i].n_range_end, 
											lin_solver_paramslist[i].n_range_step_size)} 
			  });
       }
	}
	if (lin_solver_paramslist[0].mode == MODE_DISCRETE) {
	   for (int i = 0; i < lin_solver_paramslist[0].num_tests ; i++)
	   {
		   b->Args({lin_solver_paramslist[i].m, lin_solver_paramslist[i].n});
	   }
	}

	if (lin_solver_paramslist[0].mode == MODE_COMBINATIONAL) {
	   for (int i = 0; i < lin_solver_paramslist[0].num_tests ; i++)
	   {
		  for (int j = 0; j < lin_solver_paramslist[0].num_tests ; j++){
			b->Args({lin_solver_paramslist[i].m, lin_solver_paramslist[j].n});
		  }
	   }
	}

}

BENCHMARK_REGISTER_F(dgetrf_test, dgetrf1)->Apply(CustomArguments); 

