// Gbench headers
#include "benchmark/benchmark.h"
#include "gtest/gtest.h"

// Test suite headers
#include "../Flame_gbench_main.h"
#include "../Flame_gbench_aux.h"

using namespace std;

/* Macros */
#define geqr_parameters_free() \
   if (A!=NULL){ \
      free(A); \
   }\
   if (Tau!=NULL){ \
      free(Tau);  \
   }\
   if (Work!=NULL){ \
      free(Work);  \
   }\

#define geqr_parameters_free_2() \
   if (data->A!=NULL){ \
      free(data->A); \
   }\
   if (data->Tau!=NULL){ \
      free(data->Tau);  \
   }\
   if (data->Work!=NULL){ \
      free(data->Work);  \
   }\

/* Begin geqr_parameters  class definition */
template <class T>

class geqr_parameters{
   public:
      int m, n, lda,lwork,tsize,info;
      int forced_loop_count;

      /* Local arrays */
      T *A,*Tau,*Work;
   public:
    ~geqr_parameters ();
    geqr_parameters ( int nrow, int ncol)
    {
      m = nrow;
      n = ncol; // set test matrix size
      lda = m;
      forced_loop_count = 1;
      lwork=-1;
      tsize=-1;
      info=0;
      if ( (m > FLAME_BIG_MATRIX_SIZE) || (n > FLAME_BIG_MATRIX_SIZE )){
         forced_loop_count = FLAME_GBENCH_FORCED_ITERATION_COUNT;
      }

    /* Memory allocation of the buffers */
      A = (T *)malloc(m*n*sizeof(T)) ;
      Work= NULL;
      Tau = NULL;
      if (A==NULL){
        printf("error of memory allocation. Exiting ...\n");
        geqr_parameters_free();
        exit(0);
      }
    /* Initialization of input Buffers */
      Flame_gbench_init_buffer_rand( A, m*n );


    } /* end of Constructor  */
};  /* end of geqr_parameters  class definition */



/* Destructor definition  'geqr_parameters' class  */
template <class T>
geqr_parameters<T>:: ~geqr_parameters ()
{
#if FLAME_TEST_VERBOSE
   printf(" geqr_parameters object: destructor invoked. \n");
#endif
   geqr_parameters_free ();
}

/* Fixture definition  */
template <class T>
class geqr_test : public ::benchmark::Fixture {
 public:
   int m,n;
   std::unique_ptr<geqr_parameters<T>> data;

 void SetUp(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
    m = static_cast<size_t>(state.range(0));
    n = static_cast<size_t>(state.range(1));

    assert(data.get() == nullptr);
    data.reset(new geqr_parameters<T>(m,n));

  }

  void TearDown(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      assert(data.get() != nullptr);
      data.reset();
  }

  ~geqr_test() { assert(data == nullptr); }

};

// Customized argument generator for geqr API
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

/* Template Fixture for dgeqr API  */
BENCHMARK_TEMPLATE_DEFINE_F(geqr_test, dgeqr, double)(benchmark::State &state) {
  assert(data.get() != nullptr);
  double work_query = -1;
  double t_work_query = -1;
  /* Query optimal working array size */
  dgeqr_(&data->m, &data->n,data->A,&data->lda,&t_work_query,&data->tsize,&work_query,&data->lwork,&data->info);

  data->lwork = (int)work_query;
  data->tsize = (int)t_work_query;
  data->lwork = fla_max(data->n,data->lwork);
  data->tsize = fla_max((data->n) + 5,data->tsize);
  data->Work = (double *)malloc( sizeof(double) * (data->lwork));
  data->Tau = (double *)malloc( sizeof(double) * (data->tsize));

  if ((data->Work == NULL) || (data->Tau == NULL)){
        printf("error of memory allocation. in Work Exiting ...\n");
        geqr_parameters_free_2();
        exit(0);
    }
  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the geqr API for benchmarking */
        dgeqr_(&data->m, &data->n,data->A,&data->lda,data->Tau,&data->tsize,data->Work,&data->lwork,&data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument dgeqr is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(geqr_test, dgeqr)->Apply(CustomArguments);

/* Template Fixture for sgeqr API  */
BENCHMARK_TEMPLATE_DEFINE_F(geqr_test, sgeqr, float)(benchmark::State &state) {
  assert(data.get() != nullptr);
  float work_query = -1;
  float t_work_query = -1;
  /* Query optimal working array size */
  sgeqr_(&data->m, &data->n,data->A,&data->lda,&t_work_query,&data->tsize,&work_query,&data->lwork,&data->info);
  
  data->lwork = (int)work_query;
  data->tsize = (int)t_work_query;
  data->lwork = fla_max(data->n,data->lwork);
  data->tsize = fla_max((data->n) + 5,data->tsize);
  data->Work = (float *)malloc( sizeof(float) * (data->lwork));
  data->Tau = (float *)malloc( sizeof(float) * (data->tsize));
  if ((data->Work == NULL) || (data->Tau == NULL) ){
        printf("error of memory allocation. in Work Exiting ...\n");
        geqr_parameters_free_2();
        exit(0);
    }
  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the geqr API for benchmarking */
        sgeqr_(&data->m, &data->n,data->A,&data->lda,data->Tau,&data->tsize,data->Work,&data->lwork,&data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument sgeqr is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(geqr_test, sgeqr)->Apply(CustomArguments);

/* Template Fixture for cgeqr API  */
BENCHMARK_TEMPLATE_DEFINE_F(geqr_test, cgeqr, lapack_complex_float)(benchmark::State &state) {
   assert(data.get() != nullptr);
 float work_query = -1;
 float t_work_query = -1;
  /* Query optimal working array size */
  cgeqr_(&data->m, &data->n,data->A,&data->lda,(lapack_complex_float *)&t_work_query,&data->tsize,
         (lapack_complex_float *)&work_query,&data->lwork,&data->info);
  
  data->lwork = (int)work_query;
  data->tsize = (int)t_work_query;
  data->lwork = fla_max(data->n,data->lwork);
  data->tsize = fla_max((data->n) + 5,data->tsize);
  data->Work = (lapack_complex_float *)malloc( sizeof(lapack_complex_float) * (data->lwork));
  data->Tau = (lapack_complex_float *)malloc( sizeof(lapack_complex_float) * (data->tsize));
  if ((data->Work == NULL) || (data->Tau == NULL) ){
        printf("error of memory allocation. in Work Exiting ...\n");
        geqr_parameters_free_2();
        exit(0);
    }
  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the geqr API for benchmarking */
        cgeqr_(&data->m, &data->n,data->A,&data->lda,data->Tau,&data->tsize,data->Work,&data->lwork,&data->info);
        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument cgeqr is wrong\n", data->info);
        }
    }
    
    
  }
}

BENCHMARK_REGISTER_F(geqr_test, cgeqr)->Apply(CustomArguments);

/* Template Fixture for zgeqr API  */
BENCHMARK_TEMPLATE_DEFINE_F(geqr_test, zgeqr, lapack_complex_double)(benchmark::State &state) {
  assert(data.get() != nullptr);
 double work_query = -1;
 double t_work_query = -1;
  /* Query optimal working array size */
  zgeqr_(&data->m, &data->n,data->A,&data->lda,(lapack_complex_double *)&t_work_query,&data->tsize,(lapack_complex_double *)&work_query,&data->lwork,&data->info);
  
  data->lwork = (int)work_query;
  data->tsize = (int)t_work_query;
  data->lwork = fla_max(data->n,data->lwork);
  data->tsize = fla_max((data->n) + 5,data->tsize);
  data->Work = (lapack_complex_double *)malloc( sizeof(lapack_complex_double) * (data->lwork));
  data->Tau = (lapack_complex_double *)malloc( sizeof(lapack_complex_double) * (data->tsize));
  if ((data->Work == NULL) || (data->Tau == NULL) ){
        printf("error of memory allocation. in Work Exiting ...\n");
        geqr_parameters_free_2();
        exit(0);
    }
  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the geqr API for benchmarking */
        zgeqr_(&data->m, &data->n,data->A,&data->lda,data->Tau,&data->tsize,data->Work,&data->lwork,&data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument zgeqr is wrong\n", data->info);

        }
    }
  }
}

BENCHMARK_REGISTER_F(geqr_test, zgeqr)->Apply(CustomArguments);