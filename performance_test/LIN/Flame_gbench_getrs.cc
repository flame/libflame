// Gbench headers
#include "benchmark/benchmark.h"
#include "gtest/gtest.h"

// Test suite headers
#include "../Flame_gbench_main.h"
#include "../Flame_gbench_aux.h"

using namespace std;

/* Macros */
#define getrs_parameters_free() \
   if (A!=NULL){ \
      free(A); \
   }\
   if (B!=NULL){ \
      free(B); \
   }\
   if (ipiv!=NULL){ \
      free(ipiv);  \
   }

/* Begin getrs_parameters  class definition */
template <class T>

class getrs_parameters{
   public:
      int n, lda, ldb, info;
      int forced_loop_count;
	  
      char trans; // 'N' or 'T' or 'C'. Indicates the form of the equations:
      int nrhs; // The number of right-hand sides

      /* Local arrays */
      T *A, *B;
      int *ipiv;
   public:
    ~getrs_parameters ();
    getrs_parameters ( int n_i, int nrhs_i, char trans_i )
    {
      n = n_i; // set test matrix size
	  nrhs = nrhs_i;
      lda = n;
	  ldb = n;
	  trans = trans_i;
      forced_loop_count = 1;

      if (n > FLAME_BIG_MATRIX_SIZE ){
         forced_loop_count = FLAME_GBENCH_FORCED_ITERATION_COUNT;
      }

    /* Memory allocation of the buffers */
      A = (T *)malloc(n*n*sizeof(T)) ;
      B = (T *)malloc(n*nrhs*sizeof(T)) ;
      ipiv = (int *)malloc(n*sizeof(int));
      if ((ipiv==NULL) || (A==NULL) || (B==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        getrs_parameters_free();
        exit(0);
      }
    /* Initialization of input Buffers */
      Flame_gbench_init_buffer_rand( A, n*n );
      Flame_gbench_init_buffer_rand( B, n*nrhs );

    } /* end of Constructor  */
};  /* end of getrs_parameters  class definition */



/* Destructor definition  'getrs_parameters' class  */
template <class T>
getrs_parameters<T>:: ~getrs_parameters ()
{
#if FLAME_TEST_VERBOSE
   printf(" getrs_parameters object: destructor invoked. \n");
#endif
   getrs_parameters_free ();
}

/* Fixture definition  */
template <class T>
class getrs_test : public ::benchmark::Fixture {
 public:
   int n, nrhs;
   char trans;
   std::unique_ptr<getrs_parameters<T>> data;

 void SetUp(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
    n = static_cast<size_t>(state.range(0));
    nrhs = static_cast<size_t>(state.range(1));
    trans = static_cast<char>(state.range(2));

    assert(data.get() == nullptr);
    data.reset(new getrs_parameters<T>(n, nrhs, trans));

  }

  void TearDown(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      assert(data.get() != nullptr);
      data.reset();
  }

  ~getrs_test() { assert(data == nullptr); }

};

// Customized argument generator for getrs API
static void CustomArguments(benchmark::internal::Benchmark* b) {
   if (lin_solver_paramslist[0].mode == MODE_RANGE) {
       for (int i = 0; i < lin_solver_paramslist[0].num_tests ; i++){
           b->ArgsProduct({
              {benchmark::CreateDenseRange( lin_solver_paramslist[i].m_range_start,
                                            lin_solver_paramslist[i].m_range_end,
                                            lin_solver_paramslist[i].m_range_step_size)} ,

              {lin_solver_paramslist[i].nrhs},
              {'N', 'T', 'C'},
              });
       }
    }

    if (lin_solver_paramslist[0].mode == MODE_DISCRETE) {
       for (int i = 0; i < lin_solver_paramslist[0].num_tests ; i++)
       {
           b->Args({lin_solver_paramslist[i].m, lin_solver_paramslist[i].nrhs,
		            lin_solver_paramslist[i].transr});
       }
    }

    if (lin_solver_paramslist[0].mode == MODE_COMBINATIONAL) {
       for (int i = 0; i < lin_solver_paramslist[0].num_tests ; i++)
       {
          for (int j = 0; j < lin_solver_paramslist[0].num_tests ; j++){
			  for (int k = 0; k < lin_solver_paramslist[0].num_tests ; k++){
                 b->Args({lin_solver_paramslist[i].m, lin_solver_paramslist[j].nrhs,
                          lin_solver_paramslist[j].transr});
			  }
          }
       }
    }
}

/* Template Fixture for dgetrs API  */
BENCHMARK_TEMPLATE_DEFINE_F(getrs_test, dgetrs, double)(benchmark::State &state) {
  assert(data.get() != nullptr);
  
  /* Prepare for getrs by prior invoking getrf */
  dgetrf_(&data->n, &data->n, data->A, &data->lda, 
		         data->ipiv, &data->info);

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the getrs API for benchmarking */
        dgetrs_(&data->trans, &data->n, &data->nrhs, data->A, &data->lda, 
		         data->ipiv, data->B, &data->ldb, &data->info);
        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument dgetrs is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(getrs_test, dgetrs)->Apply(CustomArguments);

/* Template Fixture for sgettrf API  */
BENCHMARK_TEMPLATE_DEFINE_F(getrs_test, sgetrs, float)(benchmark::State &state) {
  assert(data.get() != nullptr);

  /* Prepare for getrs by prior invoking getrf */
  sgetrf_(&data->n, &data->n, data->A, &data->lda, 
		         data->ipiv, &data->info);

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the getrs API for benchmarking */
        sgetrs_(&data->trans, &data->n, &data->nrhs, data->A, &data->lda, 
		         data->ipiv, data->B, &data->ldb, &data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument sgetrs is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(getrs_test, sgetrs)->Apply(CustomArguments);

/* Template Fixture for cgetrs API  */
BENCHMARK_TEMPLATE_DEFINE_F(getrs_test, cgetrs, float _Complex)(benchmark::State &state) {
   assert(data.get() != nullptr);

  /* Prepare for getrs by prior invoking getrf */
  cgetrf_(&data->n, &data->n, data->A, &data->lda, 
		         data->ipiv, &data->info);

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the getrs API for benchmarking */
        cgetrs_(&data->trans, &data->n, &data->nrhs, data->A, &data->lda, 
		         data->ipiv, data->B, &data->ldb, &data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument cgetrs is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(getrs_test, cgetrs)->Apply(CustomArguments);

/* Template Fixture for zgetrs API  */
BENCHMARK_TEMPLATE_DEFINE_F(getrs_test, zgetrs, double _Complex)(benchmark::State &state) {
  assert(data.get() != nullptr);

  /* Prepare for getrs by prior invoking getrf */
  zgetrf_(&data->n, &data->n, data->A, &data->lda, 
		         data->ipiv, &data->info);

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the getrs API for benchmarking */
        zgetrs_(&data->trans, &data->n, &data->nrhs, data->A, &data->lda, 
		         data->ipiv, data->B, &data->ldb, &data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument zgetrs is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(getrs_test, zgetrs)->Apply(CustomArguments);
