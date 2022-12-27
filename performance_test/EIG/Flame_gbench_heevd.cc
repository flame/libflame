// Gbench headers
#include "benchmark/benchmark.h"
#include "gtest/gtest.h"

// Test suite headers
#include "../Flame_gbench_main.h"
#include "../Flame_gbench_aux.h"

using namespace std;

/* Macros */
#define heevd_parameters_free() \
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
   if (Rwork!=NULL){ \
      free(Rwork);  \
   } \

#define heevd_parameters_free_2() \
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
   if (data->Rwork!=NULL){ \
      free(data->Rwork);  \
   } \
   


/* Begin heevd_parameters  class definition */
template <class T,class S>

class heevd_parameters{
   public:
      char jobz;
      char uplo; 
      int  n,lda,liwork,lwork,lrwork,info;
      int forced_loop_count;

      /* Local arrays */
      T *A,*Work;
      int *iwork;
      S *Rwork,*W;
   public:
    ~heevd_parameters ();
    heevd_parameters ( char jobz_i,char uplo_i, int n_i)
    {
      jobz=jobz_i;
      uplo=uplo_i;
      n = n_i; // set test matrix size
      lda = n;
      liwork= -1;
      lwork= -1;
      lrwork = -1;
      forced_loop_count = 1;

      if (n > FLAME_BIG_MATRIX_SIZE ){
         forced_loop_count = FLAME_GBENCH_FORCED_ITERATION_COUNT;
      }

    /* Memory allocation of the buffers */
         A = (T *)malloc((n*lda)*sizeof(T)) ;
         W = (S *)malloc((n)*sizeof(S)) ;
         Work=NULL;
         iwork=NULL;
         Rwork =NULL;
      if ((A==NULL) || (W==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        heevd_parameters_free();
        exit(0);
      }
    /* Initialization of input Buffers */
      Flame_gbench_init_buffer_rand( A, n*lda);

    } /* end of Constructor  */
};  /* end of heevd_parameters  class definition */



/* Destructor definition  'heevd_parameters' class  */
template <class T,class S>
heevd_parameters<T,S>:: ~heevd_parameters ()
{
#if FLAME_TEST_VERBOSE
   printf(" heevd_parameters object: destructor invoked. \n");
#endif
   heevd_parameters_free ();
}

/* Fixture definition  */
template <class T,class S>
class heevd_test : public ::benchmark::Fixture {
 public:
   int n;
   char jobz,uplo;
   std::unique_ptr<heevd_parameters<T,S>> data;

 void SetUp(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      n   = static_cast<size_t>(state.range(0));
      jobz = (char)static_cast<size_t>(state.range(1));  
      uplo = (char)static_cast<size_t>(state.range(2));
          
    assert(data.get() == nullptr);
    data.reset(new heevd_parameters<T,S>(jobz,uplo,n));

  }

  void TearDown(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      assert(data.get() != nullptr);
      data.reset();
  }

  ~heevd_test() { assert(data == nullptr); }

};

// Customized argument generator for heevd API
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

/* Template Fixture for cheevd API  */
BENCHMARK_TEMPLATE_DEFINE_F(heevd_test, cheevd,lapack_complex_float,float )(benchmark::State &state) {
  assert(data.get() != nullptr);
      int  iwork_query=-1;
      float work_query=-1;
      float rwork_query=-1;
    /* Query optimal working array(s) size */
   cheevd_(&data->jobz, &data->uplo,&data->n,data->A, &data->lda,data->W,(lapack_complex_float *) &work_query, &data->lwork,
   &rwork_query,&data->lrwork,&iwork_query, &data->liwork, &data->info);

    data->liwork = (int)iwork_query;
    data->lwork = (int)work_query;
    data->lrwork = (int)rwork_query;
    if(data->jobz == 'N')
    {
        data->liwork=fla_max(1,data->liwork);
        data->lwork=fla_max((data->n)+1,data->lwork);
        data->lrwork = fla_max(data->n,data->lrwork);
    }
    else{
        data->liwork=fla_max(3+(5*(data->n)),data->liwork);
        data->lwork=fla_max((2*(data->n))+((data->n)^2),data->lwork);
        data->lrwork=fla_max(1+(5*(data->n))+(2*((data->n)^2)),data->lrwork);
    }
    /* Allocate memory for work arrays */
    data->iwork = (int*)malloc( sizeof(int) * (data->liwork) );
    data->Work = (lapack_complex_float *)malloc( sizeof(lapack_complex_float) * (data->lwork));
    data->Rwork = (float *)malloc( sizeof(float) * (data->lrwork) );

    if ((data->iwork == NULL) || (data->Work==NULL)|| (data->Rwork==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        heevd_parameters_free_2();
        exit(0);
    }
    

  for (auto _ : state) {
    
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the heevd API for benchmarking */
        cheevd_(&data->jobz, &data->uplo,&data->n,data->A, &data->lda, data->W,data->Work,&data->lwork,data->Rwork,&data->lrwork,
        data->iwork,&data->liwork,&data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument cheevd is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(heevd_test,cheevd)->Apply(CustomArguments);

/* Template Fixture for zheevd API  */
BENCHMARK_TEMPLATE_DEFINE_F(heevd_test, zheevd,lapack_complex_double,double )(benchmark::State &state) {
  assert(data.get() != nullptr);
      int  iwork_query=-1;
      double work_query=-1;
      double rwork_query=-1;
    /* Query optimal working array(s) size */
   zheevd_(&data->jobz, &data->uplo,&data->n,data->A, &data->lda,data->W,(lapack_complex_double *) &work_query, &data->lwork,
   &rwork_query,&data->lrwork,&iwork_query, &data->liwork, &data->info);

    data->liwork = (int)iwork_query;
    data->lwork = (int)work_query;
    data->lrwork = (int)rwork_query;
    if(data->jobz == 'N')
    {
        data->liwork=fla_max(1,data->liwork);
        data->lwork=fla_max((data->n)+1,data->lwork);
        data->lrwork = fla_max(data->n,data->lrwork);
    }
    else{
        data->liwork=fla_max(3+(5*(data->n)),data->liwork);
        data->lwork=fla_max((2*(data->n))+((data->n)^2),data->lwork);
        data->lrwork=fla_max(1+(5*(data->n))+(2*((data->n)^2)),data->lrwork);
    }
    /* Allocate memory for work arrays */
    data->iwork = (int*)malloc( sizeof(int) * (data->liwork) );
    data->Work = (lapack_complex_double *)malloc( sizeof(lapack_complex_double) * (data->lwork));
    data->Rwork = (double *)malloc( sizeof(double) * (data->lrwork) );

    if ((data->iwork == NULL) || (data->Work==NULL)|| (data->Rwork==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        heevd_parameters_free_2();
        exit(0);
    }
    

  for (auto _ : state) {

    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the heevd API for benchmarking */
        zheevd_(&data->jobz, &data->uplo,&data->n,data->A, &data->lda, data->W,data->Work,&data->lwork,data->Rwork,&data->lrwork,
        data->iwork,&data->liwork,&data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument zheevd is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(heevd_test,zheevd)->Apply(CustomArguments);