// Gbench headers
#include "benchmark/benchmark.h"
#include "gtest/gtest.h"

// Test suite headers
#include "../Flame_gbench_main.h"
#include "../Flame_gbench_aux.h"

using namespace std;

/* Macros */
#define gecon_parameters_free() \
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

/* Begin gecon_parameters  class definition */
template <class T, class T2>

class gecon_parameters{
   public:
   
   char norm;
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
    ~gecon_parameters ();
    gecon_parameters ( int n_, char norm_ )
    {
      n = n_; // set test matrix size
      lda = n;
	  norm = norm_;
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
        gecon_parameters_free();
        exit(0);
      }
    /* Initialization of input Buffers */
      Flame_gbench_init_buffer_rand( A, n*n );
	  
	  if((norm == '1') || (norm == 'O')) {
		  anorm = n;
	  }
	  else{
		  anorm = 1;
	  }

    } /* end of Constructor  */
};  /* end of gecon_parameters  class definition */



/* Destructor definition  'gecon_parameters' class  */
template <class T, class T2>
gecon_parameters<T, T2>:: ~gecon_parameters ()
{
#if FLAME_TEST_VERBOSE
   printf(" gecon_parameters object: destructor invoked. \n");
#endif
   gecon_parameters_free ();
}

/* Fixture definition  */
template <class T, class T2>
class gecon_test : public ::benchmark::Fixture {
 public:
   int n;
   char norm;
   std::unique_ptr<gecon_parameters<T, T2>> data;

 void SetUp(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
    n = static_cast<size_t>(state.range(0));
    norm = (char)static_cast<size_t>(state.range(1));
    assert(data.get() == nullptr);
    data.reset(new gecon_parameters<T, T2>(n,norm));

  }

  void TearDown(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      assert(data.get() != nullptr);
      data.reset();
  }

  ~gecon_test() { assert(data == nullptr); }

};

// Customized argument generator for gecon API
static void CustomArguments(benchmark::internal::Benchmark* b) {
   if (lin_solver_paramslist[0].mode == MODE_RANGE) {
       for (int i = 0; i < lin_solver_paramslist[0].num_tests ; i++){
           b->ArgsProduct({
              {benchmark::CreateDenseRange( lin_solver_paramslist[i].m_range_start,
                                            lin_solver_paramslist[i].m_range_end,
                                            lin_solver_paramslist[i].m_range_step_size)} ,
                                            {'1', 'I'}

              });
       }
    }

    if (lin_solver_paramslist[0].mode == MODE_DISCRETE) {
       for (int i = 0; i < lin_solver_paramslist[0].num_tests ; i++)
       {
           b->Args({lin_solver_paramslist[i].m, lin_solver_paramslist[i].norm_gbcon});
       }
    }

    if (lin_solver_paramslist[0].mode == MODE_COMBINATIONAL) {
       for (int i = 0; i < lin_solver_paramslist[0].num_tests ; i++)
       {
          for (int j = 0; j < lin_solver_paramslist[0].num_tests ; j++){
            b->Args({lin_solver_paramslist[i].m, lin_solver_paramslist[j].norm_gbcon});
          }
       }
    }
}

/* Template Fixture for dgecon API  */
BENCHMARK_TEMPLATE_DEFINE_F(gecon_test, dgecon, double, double)(benchmark::State &state) {
  assert(data.get() != nullptr);

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the gecon API for benchmarking */
        dgecon_(&data->norm, &data->n,  data->A, &data->lda, &data->anorm, 
		        &data->rcond, data->Work, data->IWork, &data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument dgecon is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(gecon_test, dgecon)->Apply(CustomArguments);

/* Template Fixture for sgettrf API  */
BENCHMARK_TEMPLATE_DEFINE_F(gecon_test, sgecon, float, float)(benchmark::State &state) {
  assert(data.get() != nullptr);

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the gecon API for benchmarking */
        sgecon_(&data->norm, &data->n,  data->A, &data->lda, &data->anorm, 
		        &data->rcond, data->Work, data->IWork, &data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument sgecon is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(gecon_test, sgecon)->Apply(CustomArguments);


/* Template Fixture for cgecon API  */
BENCHMARK_TEMPLATE_DEFINE_F(gecon_test, cgecon, float _Complex, float)(benchmark::State &state) {
   assert(data.get() != nullptr);

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the gecon API for benchmarking */
        cgecon_(&data->norm, &data->n,  data->A, &data->lda, &data->anorm, 
		        &data->rcond, data->Work, data->RWork, &data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument cgecon is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(gecon_test, cgecon)->Apply(CustomArguments);

/* Template Fixture for zgecon API  */
BENCHMARK_TEMPLATE_DEFINE_F(gecon_test, zgecon, double _Complex, double)(benchmark::State &state) {
  assert(data.get() != nullptr);

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the gecon API for benchmarking */
        zgecon_(&data->norm, &data->n,  data->A, &data->lda, &data->anorm, 
		        &data->rcond, data->Work, data->RWork, &data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument zgecon is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(gecon_test, zgecon)->Apply(CustomArguments);
