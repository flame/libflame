// Gbench headers
#include "benchmark/benchmark.h"
#include "gtest/gtest.h"

// Test suite headers
#include "../Flame_gbench_main.h"
#include "../Flame_gbench_aux.h"

using namespace std;

/* Macros */
#define potrf_free() \
       free (a   )

/* Begin potrf_parameters  class definition */
template <class T>

class potrf_parameters{

   public:
      int forced_loop_count;
      /* Input parameters */
      char uplo; // upper or lower triangular part of A is stored
      int n; // No of rows,Columns
      int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      T *a; //The array ab contains the matrix A

      /* Return Values */
      int info;

   public:
      potrf_parameters ( char uplo_i, int n_i, int lda_i);
      ~potrf_parameters ();
};  /* end of potrf_parameters  class definition */

template <class T>
/* Constructor potrf_parameters definition */
potrf_parameters<T>:: potrf_parameters ( char uplo_i,int n_i,
                                       int lda_i) {
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
      forced_loop_count = 1;

      if ( n > FLAME_BIG_MATRIX_SIZE ){
         forced_loop_count = FLAME_GBENCH_FORCED_ITERATION_COUNT;
      }

    /* Memory allocation of the buffers */
    a = (T *)malloc((lda*n)*sizeof(T)) ;

    if(a==NULL) {
       potrf_free();
       printf(" potrf_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }
    /* Initialization of input Buffers */
      Flame_gbench_init_buffer_rand( a, lda*n );

   } /* end of Constructor  */

template <class T>
potrf_parameters<T>:: ~potrf_parameters ()
{
#if FLAME_TEST_VERBOSE
   printf(" potrf_parameters object: destructor invoked. \n");
#endif
   potrf_free();
}

/* Fixture definition  */
template <class T>

class potrf_test : public ::benchmark::Fixture {
 public:
   int n, lda;
   char uplo;
   std::unique_ptr<potrf_parameters<T>> data;

 void SetUp(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
    n = static_cast<size_t>(state.range(0));
    uplo = (char)static_cast<size_t>(state.range(1));
    assert(data.get() == nullptr);
    data.reset(new potrf_parameters<T>(uplo, n, n));

  }

  void TearDown(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      assert(data.get() != nullptr);
      data.reset();
  }

  ~potrf_test() { assert(data == nullptr); }

};

// Customized argument generator for potrf API
static void CustomArguments(benchmark::internal::Benchmark* b) {
    if (lin_solver_paramslist[0].mode == MODE_RANGE) {
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

/* Template Fixture for dpotrf API  */
BENCHMARK_TEMPLATE_DEFINE_F(potrf_test, dpotrf, double)(benchmark::State &state) {
  assert(data.get() != nullptr);

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
    /*  Invoke the POTRF api for benchmarking  */
    dpotrf_(&data->uplo, &data->n,  data->a, &data->lda, &data->info);

    if( data->info < 0 ) {
        printf( "\n warning: The i:%d th argument with dpotrf is wrong\n",
                    data->info );
    }
    }
  }
}

BENCHMARK_REGISTER_F(potrf_test, dpotrf)->Apply(CustomArguments);

/* Template Fixture for spotrf API  */
BENCHMARK_TEMPLATE_DEFINE_F(potrf_test, spotrf, float)(benchmark::State &state) {
  assert(data.get() != nullptr);

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
    /*  Invoke the POTRF api for benchmarking  */
    spotrf_(&data->uplo, &data->n,  data->a, &data->lda, &data->info);

    if( data->info < 0 ) {
        printf( "\n warning: The i:%d th argument with spotrf is wrong\n",
                    data->info );
    }
    }
  }
}

BENCHMARK_REGISTER_F(potrf_test, spotrf)->Apply(CustomArguments);

/* Template Fixture for cpotrf API  */
BENCHMARK_TEMPLATE_DEFINE_F(potrf_test, cpotrf, float _Complex)(benchmark::State &state) {
  assert(data.get() != nullptr);

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
    /*  Invoke the POTRF api for benchmarking  */
    cpotrf_(&data->uplo, &data->n,  data->a, &data->lda, &data->info);

    if( data->info < 0 ) {
        printf( "\n warning: The i:%d th argument with cpotrf is wrong\n",
                    data->info );
    }
    }
  }
}

BENCHMARK_REGISTER_F(potrf_test, cpotrf)->Apply(CustomArguments);

/* Template Fixture for zpotrf API  */
BENCHMARK_TEMPLATE_DEFINE_F(potrf_test, zpotrf, double _Complex)(benchmark::State &state) {
  assert(data.get() != nullptr);

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
    /*  Invoke the POTRF api for benchmarking  */
    zpotrf_(&data->uplo, &data->n,  data->a, &data->lda, &data->info);

    if( data->info < 0 ) {
        printf( "\n warning: The i:%d th argument with zpotrf is wrong\n",
                    data->info );
    }
    }
  }
}

BENCHMARK_REGISTER_F(potrf_test, zpotrf)->Apply(CustomArguments);
