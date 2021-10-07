// Gbench headers
#include "benchmark/benchmark.h"
#include "gtest/gtest.h"

// Test suite headers
#include "../Flame_gbench_main.h"
#include "../Flame_gbench_aux.h"

using namespace std;

/* Macros */
#define pocon_parameters_free() \
   if (A!=NULL){ \
      free(A); \
   }\
   if (IWork!=NULL){ \
      free(IWork); \
   }\
   if (RWork!=NULL){ \
      free(RWork); \
   }\
   if (Work!=NULL){ \
      free(Work);  \
   }

/* Begin pocon_parameters  class definition */
template <class T, class T2>

class pocon_parameters{
   public:
   
   char uplo;
   T2 anorm;
   T2 rcond;
   int n, lda, info;
   int forced_loop_count;
	  
   /* Local arrays */
   T  *A;
   T *Work;
   T2 *RWork;
   int *IWork;

   public:
    ~pocon_parameters ();
    pocon_parameters ( int n_, char uplo_ )
    {
      n = n_; // set test matrix size
      lda = n;
	  uplo = uplo_;
      forced_loop_count = 1;

      if ( n > FLAME_BIG_MATRIX_SIZE ){
         forced_loop_count = FLAME_GBENCH_FORCED_ITERATION_COUNT;
      }

    /* Memory allocation of the buffers */
      A = (T *)malloc(n*n*sizeof(T)) ;
      Work = (T *)malloc(4*n*sizeof(T)) ;
	  
      IWork = (int *)malloc(n*sizeof(int));
      RWork = (T2 *)malloc(2*n*sizeof(T2)) ;
      if ((IWork==NULL) || (A==NULL) || (Work==NULL)|| (RWork==NULL) ){
        printf("error of memory allocation. Exiting ...\n");
        pocon_parameters_free();
        exit(0);
      }
    /* Initialization of input Buffers */
      Flame_gbench_init_buffer_rand( A, n*n );
	  
	  if((uplo == '1') || (uplo == 'O')) {
		  anorm = n;
	  }
	  else{
		  anorm = 1;
	  }

    } /* end of Constructor  */
};  /* end of pocon_parameters  class definition */



/* Destructor definition  'pocon_parameters' class  */
template <class T, class T2>
pocon_parameters<T, T2>:: ~pocon_parameters ()
{
#if FLAME_TEST_VERBOSE
   printf(" pocon_parameters object: destructor invoked. \n");
#endif
   pocon_parameters_free ();
}

/* Fixture definition  */
template <class T, class T2>
class pocon_test : public ::benchmark::Fixture {
 public:
   int n;
   char uplo;
   std::unique_ptr<pocon_parameters<T, T2>> data;

 void SetUp(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
    n = static_cast<size_t>(state.range(0));
    uplo = (char)static_cast<size_t>(state.range(1));
    assert(data.get() == nullptr);
    data.reset(new pocon_parameters<T, T2>(n,uplo));

  }

  void TearDown(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      assert(data.get() != nullptr);
      data.reset();
  }

  ~pocon_test() { assert(data == nullptr); }

};

// Customized argument generator for pocon API
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

/* Template Fixture for dpocon API  */
BENCHMARK_TEMPLATE_DEFINE_F(pocon_test, dpocon, double, double)(benchmark::State &state) {
  assert(data.get() != nullptr);

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the pocon API for benchmarking */
        dpocon_(&data->uplo, &data->n,  data->A, &data->lda, &data->anorm, 
		        &data->rcond, data->Work, data->IWork, &data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument dpocon is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(pocon_test, dpocon)->Apply(CustomArguments);

/* Template Fixture for sgettrf API  */
BENCHMARK_TEMPLATE_DEFINE_F(pocon_test, spocon, float, float)(benchmark::State &state) {
  assert(data.get() != nullptr);

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the pocon API for benchmarking */
        spocon_(&data->uplo, &data->n,  data->A, &data->lda, &data->anorm, 
		        &data->rcond, data->Work, data->IWork, &data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument spocon is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(pocon_test, spocon)->Apply(CustomArguments);


/* Template Fixture for cpocon API  */
BENCHMARK_TEMPLATE_DEFINE_F(pocon_test, cpocon, float _Complex, float)(benchmark::State &state) {
   assert(data.get() != nullptr);

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the pocon API for benchmarking */
        cpocon_(&data->uplo, &data->n,  data->A, &data->lda, &data->anorm, 
		        &data->rcond, data->Work, data->RWork, &data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument cpocon is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(pocon_test, cpocon)->Apply(CustomArguments);

/* Template Fixture for zpocon API  */
BENCHMARK_TEMPLATE_DEFINE_F(pocon_test, zpocon, double _Complex, double)(benchmark::State &state) {
  assert(data.get() != nullptr);

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the pocon API for benchmarking */
        zpocon_(&data->uplo, &data->n,  data->A, &data->lda, &data->anorm, 
		        &data->rcond, data->Work, data->RWork, &data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument zpocon is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(pocon_test, zpocon)->Apply(CustomArguments);
