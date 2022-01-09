// Gbench headers
#include "benchmark/benchmark.h"
#include "gtest/gtest.h"

// Test suite headers
#include "../Flame_gbench_main.h"
#include "../Flame_gbench_aux.h"

using namespace std;

/* Macros */
#define getrf_parameters_free() \
   if (A!=NULL){ \
      free(A); \
   }\
   if (ipiv!=NULL){ \
      free(ipiv);  \
   }

/* Begin getrf_parameters  class definition */
template <class T>

class getrf_parameters{
   public:
      int m, n, lda, info;
      int forced_loop_count;

      /* Local arrays */
      T *A;
      int *ipiv, *ipivref;
      T norm;
   public:
    ~getrf_parameters ();
    getrf_parameters ( int nrow, int ncol, int lda_ )
    {
      m = nrow;
      n = ncol; // set test matrix size
      lda = lda_;
      forced_loop_count = 1;

      if ( (m > FLAME_BIG_MATRIX_SIZE) || (n > FLAME_BIG_MATRIX_SIZE )){
         forced_loop_count = FLAME_GBENCH_FORCED_ITERATION_COUNT;
      }

    /* Memory allocation of the buffers */
      A = (T *)malloc(m*n*sizeof(T)) ;
      ipiv = (int *)malloc(m*sizeof(int));
      if ((ipiv==NULL) || (A==NULL) ){
        printf("error of memory allocation. Exiting ...\n");
        getrf_parameters_free();
        exit(0);
      }
    /* Initialization of input Buffers */
      Flame_gbench_init_buffer_rand( A, m*n );

    } /* end of Constructor  */
};  /* end of getrf_parameters  class definition */



/* Destructor definition  'getrf_parameters' class  */
template <class T>
getrf_parameters<T>:: ~getrf_parameters ()
{
#if FLAME_TEST_VERBOSE
   printf(" getrf_parameters object: destructor invoked. \n");
#endif
   getrf_parameters_free ();
}

/* Fixture definition  */
template <class T>
class getrf_test : public ::benchmark::Fixture {
 public:
   int m,n, lda, ldb;
   std::unique_ptr<getrf_parameters<T>> data;

 void SetUp(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
    m = static_cast<size_t>(state.range(0));
    n = static_cast<size_t>(state.range(1));

    assert(data.get() == nullptr);
    data.reset(new getrf_parameters<T>(m,n,m));

  }

  void TearDown(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      assert(data.get() != nullptr);
      data.reset();
  }

  ~getrf_test() { assert(data == nullptr); }

};

// Customized argument generator for getrf API
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

/* Template Fixture for dgetrf API  */
BENCHMARK_TEMPLATE_DEFINE_F(getrf_test, dgetrf, double)(benchmark::State &state) {
  assert(data.get() != nullptr);

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the getrf API for benchmarking */
        dgetrf_(&data->m, &data->n,  data->A, &data->lda, data->ipiv, &data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument dgetrf is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(getrf_test, dgetrf)->Apply(CustomArguments);

/* Template Fixture for sgettrf API  */
BENCHMARK_TEMPLATE_DEFINE_F(getrf_test, sgetrf, float)(benchmark::State &state) {
  assert(data.get() != nullptr);

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the getrf API for benchmarking */
        sgetrf_(&data->m, &data->n,  data->A, &data->lda, data->ipiv, &data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument sgetrf is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(getrf_test, sgetrf)->Apply(CustomArguments);

/* Template Fixture for cgetrf API  */
BENCHMARK_TEMPLATE_DEFINE_F(getrf_test, cgetrf, float _Complex)(benchmark::State &state) {
   assert(data.get() != nullptr);

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the getrf API for benchmarking */
        cgetrf_(&data->m, &data->n,  data->A, &data->lda, data->ipiv, &data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument cgetrf is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(getrf_test, cgetrf)->Apply(CustomArguments);

/* Template Fixture for zgetrf API  */
BENCHMARK_TEMPLATE_DEFINE_F(getrf_test, zgetrf, double _Complex)(benchmark::State &state) {
  assert(data.get() != nullptr);

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the getrf API for benchmarking */
        zgetrf_(&data->m, &data->n,  data->A, &data->lda, data->ipiv, &data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument zgetrf is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(getrf_test, zgetrf)->Apply(CustomArguments);
