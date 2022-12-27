// Gbench headers
#include "benchmark/benchmark.h"
#include "gtest/gtest.h"

// Test suite headers
#include "../Flame_gbench_main.h"
#include "../Flame_gbench_aux.h"

using namespace std;

/* Macros */
#define syevd_parameters_free() \
   if (A!=NULL){ \
      free(A); \
   }\
   if (W!=NULL){ \
      free(W);  \
   } \
   if (Work!=NULL){ \
      free(Work); \
   }\
   if (iwork!=NULL){ \
      free(iwork);  \
   } \

#define syevd_parameters_free_2() \
   if (data->A!=NULL){ \
      free(data->A); \
   }\
   if (data->W!=NULL){ \
      free(data->W);  \
   } \
   if (data->Work!=NULL){ \
      free(data->Work); \
   }\
   if (data->iwork!=NULL){ \
      free(data->iwork);  \
   } \


/* Begin syevd_parameters  class definition */
template <class T>

class syevd_parameters{
   public:
      char jobz;
      char uplo; 
      int  n,lda,liwork,lwork,info;
      int forced_loop_count;

      /* Local arrays */
      T *A;
      T *W,*Work;
      int *iwork;
   public:
    ~syevd_parameters ();
    syevd_parameters ( char jobz_i,char uplo_i, int n_i)
    {
      jobz=jobz_i;
      uplo=uplo_i;
      n = n_i; // set test matrix size
      lda = n;
      liwork= -1;
      lwork= -1;
      forced_loop_count = 1;

      if (n > FLAME_BIG_MATRIX_SIZE ){
         forced_loop_count = FLAME_GBENCH_FORCED_ITERATION_COUNT;
      }

    /* Memory allocation of the buffers */
         A = (T *)malloc((n*lda)*sizeof(T)) ;
         W = (T *)malloc((n)*sizeof(T)) ;
         Work=NULL;
         iwork=NULL;
      if ((A==NULL) || (W==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        syevd_parameters_free();
        exit(0);
      }
    /* Initialization of input Buffers */
      Flame_gbench_init_buffer_rand( A, n*lda);

    } /* end of Constructor  */
};  /* end of syevd_parameters  class definition */



/* Destructor definition  'syevd_parameters' class  */
template <class T>
syevd_parameters<T>:: ~syevd_parameters ()
{
#if FLAME_TEST_VERBOSE
   printf(" syevd_parameters object: destructor invoked. \n");
#endif
   syevd_parameters_free ();
}

/* Fixture definition  */
template <class T>
class syevd_test : public ::benchmark::Fixture {
 public:
   int n,lda,lwork,liwork;
   char jobz,uplo;
   std::unique_ptr<syevd_parameters<T>> data;

 void SetUp(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      n   = static_cast<size_t>(state.range(0));
      jobz = (char)static_cast<size_t>(state.range(1));  
      uplo = (char)static_cast<size_t>(state.range(2));
          
    assert(data.get() == nullptr);
    data.reset(new syevd_parameters<T>(jobz,uplo,n));

  }

  void TearDown(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      assert(data.get() != nullptr);
      data.reset();
  }

  ~syevd_test() { assert(data == nullptr); }

};

// Customized argument generator for syevd API
static void CustomArguments(benchmark::internal::Benchmark* b) {
    
    if (eig_paramslist[0].mode == MODE_RANGE) {
       for (int i = 0; i < eig_paramslist[0].num_tests ; i++){
           b->ArgsProduct({
              {benchmark::CreateDenseRange( eig_paramslist[i].n_range_start,
                                            eig_paramslist[i].n_range_end,
                                            eig_paramslist[i].n_range_step_size)},
                                            {'N', 'V'},
                                            {'U', 'L'}
              });
       }
    }
   
    
    if (eig_paramslist[0].mode == MODE_DISCRETE) {
       for (int i = 0; i < eig_paramslist[0].num_tests ; i++)
       {
           b->Args({eig_paramslist[i].n, eig_paramslist[i].jobz, eig_paramslist[i].uplo});
       }
    }
    if (eig_paramslist[0].mode == MODE_COMBINATIONAL) {
       for (int i = 0; i < eig_paramslist[0].num_tests ; i++)
       {
          for (int j = 0; j < eig_paramslist[0].num_tests ; j++){
              for(int k=0; k< eig_paramslist[0].num_tests ; k++){
            b->Args({eig_paramslist[i].n, eig_paramslist[j].jobz,eig_paramslist[k].uplo});
              }
          }
       }
    }
}

/* Template Fixture for dsyevd API  */
BENCHMARK_TEMPLATE_DEFINE_F(syevd_test, dsyevd, double)(benchmark::State &state) {
  assert(data.get() != nullptr);
      int  iwork_query=-1;
      double work_query=-1;
    /* Query optimal working array(s) size */
   dsyevd_(&data->jobz, &data->uplo,&data->n,data->A, &data->lda,data->W, &work_query, &data->lwork,
   &iwork_query, &data->liwork, &data->info);

    data->liwork = (int)iwork_query;
    data->lwork = (int)work_query;
    if(data->jobz == 'N')
    {
        data->liwork=fla_max(1,data->liwork);
        data->lwork=fla_max((2*(data->n))+1,data->lwork);
    }
    else{
        data->liwork=fla_max(3+(5*(data->n)),data->liwork);
        data->lwork=fla_max((1+(6*(data->n)))+(2*((data->n)^2)),data->lwork);
    }
    /* Allocate memory for work arrays */
    data->iwork = (int*)malloc( sizeof(int) * (data->liwork) );
    data->Work = (double*)malloc( sizeof(double) * (data->lwork));
    if ((data->iwork == NULL) || (data->Work==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        syevd_parameters_free_2();
        exit(0);
    }
    

  for (auto _ : state) {

    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the syevd API for benchmarking */
        dsyevd_(&data->jobz, &data->uplo,&data->n,data->A, &data->lda, data->W,data->Work,&data->lwork,data->iwork,&data->liwork,&data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument dsyevd is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(syevd_test, dsyevd)->Apply(CustomArguments);

/* Template Fixture for ssyevd API  */
BENCHMARK_TEMPLATE_DEFINE_F(syevd_test, ssyevd, float)(benchmark::State &state) {
    assert(data.get() != nullptr);
      int  iwork_query=-1;
      float work_query=-1;
    /* Query optimal working array(s) size */
    ssyevd_(&data->jobz, &data->uplo,&data->n,data->A, &data->lda,data->W, &work_query, &data->lwork,
    &iwork_query, &data->liwork, &data->info);
                                    
    data->liwork = (int)iwork_query;
    data->lwork = (int)work_query;
    if(data->jobz == 'N')
    {
        data->liwork=fla_max(1,data->liwork);
        data->lwork=fla_max((2*(data->n))+1,data->lwork);
    }
    else{
        data->liwork=fla_max(3+(5*(data->n)),data->liwork);
        data->lwork=fla_max((1+(6*(data->n)))+(2*((data->n)^2)),data->lwork);
    }
    /* Allocate memory for work arrays */
    data->iwork = (int*)malloc( sizeof(int) * (data->liwork) );
    data->Work = (float*)malloc( sizeof(float) * (data->lwork));
    if ((data->iwork == NULL) || (data->Work==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        syevd_parameters_free_2();
        exit(0);
    }
    for (auto _ : state) {
        /* Invoke the syevd API for benchmarking */
        for (int i=0; i<data->forced_loop_count; i++) {
            ssyevd_(&data->jobz, &data->uplo,&data->n,data->A, &data->lda, data->W,data->Work,&data->lwork,data->iwork,&data->liwork,&data->info);
        
            if( data->info < 0 ) {
                printf( "\n warning: The %d th argument ssyevd is wrong\n", data->info );
            }
        }
    }
}

BENCHMARK_REGISTER_F(syevd_test, ssyevd)->Apply(CustomArguments);