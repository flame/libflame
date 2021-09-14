// Gbench headers
#include "benchmark/benchmark.h"
#include "gtest/gtest.h"

// Test suite headers
#include "../Flame_gbench_main.h"
#include "../Flame_gbench_aux.h"

using namespace std;

/* Macros */
#define geqrf_parameters_free() \
   if (A!=NULL){ \
      free(A); \
   }\
   if (Tau!=NULL){ \
      free(Tau);  \
   }\
   if (Work!=NULL){ \
      free(Work);  \
   }\

#define geqrf_parameters_free_2() \
   if (data->A!=NULL){ \
      free(data->A); \
   }\
   if (data->Tau!=NULL){ \
      free(data->Tau);  \
   }\
   if (data->Work!=NULL){ \
      free(data->Work);  \
   }\

/* Begin geqrf_parameters  class definition */
template <class T>

class geqrf_parameters{
   public:
      int m, n, lda,lwork,info;
      int forced_loop_count;

      /* Local arrays */
      T *A,*Tau,*Work;
   public:
    ~geqrf_parameters ();
    geqrf_parameters ( int nrow, int ncol)
    {
      m = nrow;
      n = ncol; // set test matrix size
      lda = m;
      forced_loop_count = 1;
      lwork=-1;
      if ( (m > FLAME_BIG_MATRIX_SIZE) || (n > FLAME_BIG_MATRIX_SIZE )){
         forced_loop_count = FLAME_GBENCH_FORCED_ITERATION_COUNT;
      }

    /* Memory allocation of the buffers */
      A = (T *)malloc(m*n*sizeof(T)) ;
      Tau = (T *)malloc(m*sizeof(T));
      Work=NULL;
      if ((A==NULL) || (Tau==NULL) ){
        printf("error of memory allocation. Exiting ...\n");
        geqrf_parameters_free();
        exit(0);
      }
    /* Initialization of input Buffers */
      Flame_gbench_init_buffer_rand( A, m*n );


    } /* end of Constructor  */
};  /* end of geqrf_parameters  class definition */



/* Destructor definition  'geqrf_parameters' class  */
template <class T>
geqrf_parameters<T>:: ~geqrf_parameters ()
{
#if FLAME_TEST_VERBOSE
   printf(" geqrf_parameters object: destructor invoked. \n");
#endif
   geqrf_parameters_free ();
}

/* Fixture definition  */
template <class T>
class geqrf_test : public ::benchmark::Fixture {
 public:
   int m,n;
   std::unique_ptr<geqrf_parameters<T>> data;

 void SetUp(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
    m = static_cast<size_t>(state.range(0));
    n = static_cast<size_t>(state.range(1));

    assert(data.get() == nullptr);
    data.reset(new geqrf_parameters<T>(m,n));

  }

  void TearDown(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      assert(data.get() != nullptr);
      data.reset();
  }

  ~geqrf_test() { assert(data == nullptr); }

};

// Customized argument generator for geqrf API
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

/* Template Fixture for dgeqrf API  */
BENCHMARK_TEMPLATE_DEFINE_F(geqrf_test, dgeqrf, double)(benchmark::State &state) {
  assert(data.get() != nullptr);
  double work_query = -1;
  /* Query optimal working array size */
  dgeqrf_(&data->m, &data->n,data->A,&data->lda,data->Tau,&work_query,&data->lwork,&data->info);

  data->lwork = (int)work_query;
  data->lwork =max(data->n,data->lwork);
  data->Work = (double *)malloc( sizeof(double) * (data->lwork));
  if (data->Work == NULL){
        printf("error of memory allocation. in Work Exiting ...\n");
        geqrf_parameters_free_2();
        exit(0);
    }
  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the geqrf API for benchmarking */
        dgeqrf_(&data->m, &data->n,data->A,&data->lda,data->Tau,data->Work,&data->lwork,&data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument dgeqrf is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(geqrf_test, dgeqrf)->Apply(CustomArguments);

/* Template Fixture for sgeqrf API  */
BENCHMARK_TEMPLATE_DEFINE_F(geqrf_test, sgeqrf, float)(benchmark::State &state) {
  assert(data.get() != nullptr);
  float work_query = -1;
  /* Query optimal working array size */
  sgeqrf_(&data->m, &data->n,data->A,&data->lda,data->Tau,&work_query,&data->lwork,&data->info);

  data->lwork = (int)work_query;
  data->lwork =max(data->n,data->lwork);
  data->Work = (float *)malloc( sizeof(float) * (data->lwork));
  if (data->Work == NULL){
        printf("error of memory allocation. in Work Exiting ...\n");
        geqrf_parameters_free_2();
        exit(0);
    }
  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the geqrf API for benchmarking */
        sgeqrf_(&data->m, &data->n,data->A,&data->lda,data->Tau,data->Work,&data->lwork,&data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument sgeqrf is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(geqrf_test, sgeqrf)->Apply(CustomArguments);

/* Template Fixture for cgeqrf API  */
BENCHMARK_TEMPLATE_DEFINE_F(geqrf_test, cgeqrf, lapack_complex_float)(benchmark::State &state) {
   assert(data.get() != nullptr);
 float work_query = -1;
 /* Query optimal working array size */
 cgeqrf_(&data->m, &data->n,data->A,&data->lda,data->Tau,(lapack_complex_float *)&work_query,&data->lwork,&data->info);

  data->lwork = (int)work_query;
  data->lwork =max(data->n,data->lwork);
  data->Work = (lapack_complex_float *)malloc( sizeof(lapack_complex_float) * (data->lwork));
  if (data->Work == NULL){
        printf("error of memory allocation. in Work Exiting ...\n");
        geqrf_parameters_free_2();
        exit(0);
    }
  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the geqrf API for benchmarking */
        cgeqrf_(&data->m, &data->n,data->A,&data->lda,data->Tau,data->Work,&data->lwork,&data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument cgeqrf is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(geqrf_test, cgeqrf)->Apply(CustomArguments);

/* Template Fixture for zgeqrf API  */
BENCHMARK_TEMPLATE_DEFINE_F(geqrf_test, zgeqrf, lapack_complex_double)(benchmark::State &state) {
  assert(data.get() != nullptr);
 double work_query = -1;
 /* Query optimal working array size */
 zgeqrf_(&data->m, &data->n,data->A,&data->lda,data->Tau,(lapack_complex_double *)&work_query,&data->lwork,&data->info);

  data->lwork = (int)work_query;
  data->lwork =max(data->n,data->lwork);
  data->Work = (lapack_complex_double *)malloc( sizeof(lapack_complex_double) * (data->lwork));
  if (data->Work == NULL){
        printf("error of memory allocation. in Work Exiting ...\n");
        geqrf_parameters_free_2();
        exit(0);
    }
  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the geqrf API for benchmarking */
        zgeqrf_(&data->m, &data->n,data->A,&data->lda,data->Tau,data->Work,&data->lwork,&data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument zgeqrf is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(geqrf_test, zgeqrf)->Apply(CustomArguments);