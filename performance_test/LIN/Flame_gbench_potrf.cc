// Standard headers
#include <cstdio>
#include <iostream>
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
#define potrf_free() \
       free (a   )

/* Begin potrf_double_parameters  class definition */
class potrf_double_parameters{

   public:
      /* Input parameters */
      char uplo; // upper or lower triangular part of A is stored
      int n; // No of rows,Columns
      int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      double *a; //The array ab contains the matrix A

      /* Return Values */
      int info;

   public: 
      potrf_double_parameters ( char uplo_i, int n_i, int lda_i);
      ~potrf_double_parameters (); 
};  /* end of potrf_double_parameters  class definition */


/* Constructor potrf_double_parameters definition */
potrf_double_parameters:: potrf_double_parameters ( char uplo_i,int n_i,
                                       int lda_i) {
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;

    /* Memory allocation of the buffers */
    a = (double *)malloc((lda*n)*sizeof(double)) ;

    if(a==NULL) {
       potrf_free();
       printf(" potrf_double_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
	for( int i = 0; i < n; i++ ) {
		for( int j = 0; j < lda; j++ ) {
			a[i+j*lda] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
		}
	}

   } /* end of Constructor  */

potrf_double_parameters:: ~potrf_double_parameters ()
{
#if FLAME_TEST_VERBOSE
   printf(" potrf_double_parameters object: destructor invoked. \n");
#endif
   potrf_free();
}


class dpotrf_test : public ::benchmark::Fixture {
 public:
   int n, lda;
   char uplo;
   std::unique_ptr<potrf_double_parameters> data;

 void SetUp(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
    n = static_cast<size_t>(state.range(0));
    uplo = (char)static_cast<size_t>(state.range(1));
    assert(data.get() == nullptr);
    data.reset(new potrf_double_parameters(uplo, n, n));

  }

  void TearDown(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      assert(data.get() != nullptr);
      data.reset();
  }

  ~dpotrf_test() { assert(data == nullptr); }

};

BENCHMARK_DEFINE_F(dpotrf_test, dpotrf1)(benchmark::State &state) {
  assert(data.get() != nullptr);
  int forced_loop_count = 1;
  
  if ( data->n > FLAME_BIG_MATRIX_SIZE ){
	  forced_loop_count = FLAME_GBENCH_FORCED_ITERATION_COUNT;
  }
  //printf(" data->m: %d ,data->lda : %d \n", data->n,data->lda);

  for (auto _ : state) {
	for (int i=0; i<forced_loop_count; i++) {
    /*  Invoke the POTRF api for benchmarking  */
	dpotrf_(&data->uplo, &data->n,  data->a, &data->lda, &data->info);

	if( data->info < 0 ) {
		printf( "\n warning: The i:%d th argument with dpotrf is wrong\n",
					data->info );
	}
	}
  }
}

static void CustomArguments(benchmark::internal::Benchmark* b) {
	if (lin_solver_paramslist[0].mode == MODE_RANGE) {
		benchmark::internal::Benchmark* b1;
	   for (int i = 0; i < lin_solver_paramslist[0].num_tests ; i++){
	       b->ArgsProduct({ 
              {benchmark::CreateDenseRange( lin_solver_paramslist[i].m_range_start, 
											lin_solver_paramslist[i].m_range_end, 
											lin_solver_paramslist[i].m_range_step_size)} , 
											{'U', 'L'}
		   });
       }
	}
	if (lin_solver_paramslist[0].mode == MODE_DISCRETE) {
	   for (int i = 0; i < lin_solver_paramslist[0].num_tests ; i++)
	   {
		   b->Args({lin_solver_paramslist[i].m, lin_solver_paramslist[i].Uplo});
	   }
	}

	if (lin_solver_paramslist[0].mode == MODE_COMBINATIONAL) {
	   for (int i = 0; i < lin_solver_paramslist[0].num_tests ; i++)
	   {
		  for (int j = 0; j < lin_solver_paramslist[0].num_tests ; j++){
			b->Args({lin_solver_paramslist[i].m, lin_solver_paramslist[j].Uplo});
		  }
	   }
	}

}
BENCHMARK_REGISTER_F(dpotrf_test, dpotrf1)->Apply(CustomArguments);
