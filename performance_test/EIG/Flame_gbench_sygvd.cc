// Gbench headers
#include "benchmark/benchmark.h"
#include "gtest/gtest.h"

// Test suite headers
#include "../Flame_gbench_main.h"
#include "../Flame_gbench_aux.h"

using namespace std;

/* Macros */
#define sygvd_parameters_free() \
   if (A!=NULL){ \
      free(A); \
   }\
   if (B!=NULL){ \
      free(B);  \
   } \
   if (W!=NULL){ \
      free(W);  \
   } \
   if (Work!=NULL){ \
      free(Work); \
   }\
   if (iwork!=NULL){ \
      free(iwork);  \
   } \

#define sygvd_parameters_free_2() \
   if (data->A!=NULL){ \
      free(data->A); \
   }\
   if (data->B!=NULL){ \
      free(data->B);  \
   } \
   if (data->W!=NULL){ \
      free(data->W);  \
   } \
   if (data->Work!=NULL){ \
      free(data->Work); \
   }\
   if (data->iwork!=NULL){ \
      free(data->iwork);  \
   } \

/* Begin sygvd_parameters  class definition */
template <class T>

class sygvd_parameters{
   public:
      char jobz;
      char uplo; 
      int  n,lda,ldb,lwork,liwork,info,itype;
      int forced_loop_count;

      /* Local arrays */
      T *A,*B,*W;
      T *Work;
      int *iwork;
   public:
    ~sygvd_parameters ();
    sygvd_parameters ( int itype_i,char jobz_i,char uplo_i, int n_i)
    {
      itype = itype_i;
      jobz=jobz_i;
      uplo=uplo_i;
      n = n_i; // set test matrix size
      lda = n;
      ldb = n;
      lwork= -1;
      liwork= -1;
      forced_loop_count = 1;

      if (n > FLAME_BIG_MATRIX_SIZE ){
         forced_loop_count = FLAME_GBENCH_FORCED_ITERATION_COUNT;
      }

    /* Memory allocation of the buffers */
         A = (T *)malloc((n*lda)*sizeof(T)) ;
         B = (T *)malloc((n*ldb)*sizeof(T)) ;
         W = (T *)malloc((n)*sizeof(T)) ;
         Work=NULL;
         iwork=NULL;
         
      if ((A==NULL) || (B==NULL)|| (W==NULL)){ 
        printf("error of memory allocation. Exiting ...\n");
        sygvd_parameters_free();
        exit(0);
      }
    /* Initialization of input Buffers */
      Flame_gbench_init_buffer_rand( A, n*lda);
      Flame_gbench_init_buffer_rand( B, n*ldb);

    } /* end of Constructor  */
};  /* end of sygvd_parameters  class definition */



/* Destructor definition  'sygvd_parameters' class  */
template <class T>
sygvd_parameters<T>:: ~sygvd_parameters ()
{
#if FLAME_TEST_VERBOSE
   printf(" sygvd_parameters object: destructor invoked. \n");
#endif
   sygvd_parameters_free ();
}

/* Fixture definition  */
template <class T>
class sygvd_test : public ::benchmark::Fixture {
 public:
   int n,itype;
   char jobz,uplo;
   std::unique_ptr<sygvd_parameters<T>> data;

 void SetUp(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      n   = static_cast<size_t>(state.range(0));
      itype = static_cast<size_t>(state.range(1));
      jobz = (char)static_cast<size_t>(state.range(2));  
      uplo = (char)static_cast<size_t>(state.range(3));
          
    assert(data.get() == nullptr);
    data.reset(new sygvd_parameters<T>(itype,jobz,uplo,n));

  }

  void TearDown(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      assert(data.get() != nullptr);
      data.reset();
  }

  ~sygvd_test() { assert(data == nullptr); }

};

// Customized argument generator for sygvd API
static void CustomArguments(benchmark::internal::Benchmark* b) {
    
    if (eig_paramslist[0].mode == MODE_RANGE) {
       for (int i = 0; i < eig_paramslist[0].num_tests ; i++){
           b->ArgsProduct({
              {benchmark::CreateDenseRange( eig_paramslist[i].n_range_start,
                                            eig_paramslist[i].n_range_end,
                                            eig_paramslist[i].n_range_step_size)},
                                            {1,2,3},
                                            {'N', 'V'},
                                            {'U', 'L'}
              });
       }
    }
     
    if (eig_paramslist[0].mode == MODE_DISCRETE) {
       for (int i = 0; i < eig_paramslist[0].num_tests ; i++)
       {
           b->Args({eig_paramslist[i].n,eig_paramslist[i].itype, eig_paramslist[i].jobz, eig_paramslist[i].uplo});
       }
    }

    if (eig_paramslist[0].mode == MODE_COMBINATIONAL) {
       for (int i = 0; i < eig_paramslist[0].num_tests ; i++)
       {
          for (int j = 0; j < eig_paramslist[0].num_tests ; j++){
              for(int k=0; k< eig_paramslist[0].num_tests ; k++){
                  for(int l=0; l< eig_paramslist[0].num_tests ; l++){

                    b->Args({eig_paramslist[i].n, eig_paramslist[j].itype,eig_paramslist[k].jobz,eig_paramslist[l].uplo});
                  
                  }
               }
           }
        }
    }
}

/* Template Fixture for ssygvd API  */
BENCHMARK_TEMPLATE_DEFINE_F(sygvd_test, ssygvd,float)(benchmark::State &state) {
  assert(data.get() != nullptr);
      float work_query=-1;
      int iwork_query = -1;
    /* Query optimal working array(s) size */
   ssygvd_(&data->itype,&data->jobz,&data->uplo,&data->n,data->A, &data->lda, data->B,
           &data->ldb,data->W,&work_query,&data->lwork,&iwork_query,&data->liwork,&data->info); 

    data->lwork = (int)work_query;
    data->liwork = (int)iwork_query;
    if(data->jobz == 'N')
    {
        data->liwork=max(1,data->liwork);
        data->lwork=max((2*(data->n)+1),data->lwork);
    }
    else{
        data->liwork=max(3+(5*(data->n)),data->liwork);
        data->lwork=max((1+(6*(data->n))+(2*(data->n)^2)),data->lwork);
    }
    /* Allocate memory for work arrays */
    data->iwork = (int*)malloc( sizeof(int) * (data->liwork) );
    data->Work = (float *)malloc( sizeof(float) * (data->lwork));
    if ((data->Work == NULL) || (data->iwork == NULL) ){
        printf("error of memory allocation. in Work Exiting ...\n");
        sygvd_parameters_free_2();
        exit(0);
    }
    
  for (auto _ : state) {

    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the sygvd API for benchmarking */
        ssygvd_(&data->itype,&data->jobz,&data->uplo,&data->n,data->A, &data->lda, data->B,
           &data->ldb,data->W,data->Work,&data->lwork,data->iwork,&data->liwork,&data->info); 
        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument ssygvd is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(sygvd_test, ssygvd)->Apply(CustomArguments);

/* Template Fixture for dsygvd API  */
BENCHMARK_TEMPLATE_DEFINE_F(sygvd_test, dsygvd,double)(benchmark::State &state) {
  assert(data.get() != nullptr);
      double work_query=-1;
      int iwork_query = -1;
    /* Query optimal working array(s) size */
   dsygvd_(&data->itype,&data->jobz,&data->uplo,&data->n,data->A, &data->lda, data->B,
           &data->ldb,data->W,&work_query,&data->lwork,&iwork_query,&data->liwork,&data->info); 

    data->lwork = (int)work_query;
    data->liwork = (int)iwork_query;
    if(data->jobz == 'N')
    {
        data->liwork=max(1,data->liwork);
        data->lwork=max((2*(data->n)+1),data->lwork);
    }
    else{
        data->liwork=max(3+(5*(data->n)),data->liwork);
        data->lwork=max((1+(6*(data->n))+(2*(data->n)^2)),data->lwork);
    }
    /* Allocate memory for work arrays */
    data->iwork = (int*)malloc( sizeof(int) * (data->liwork) );
    data->Work = (double *)malloc( sizeof(double) * (data->lwork));
    if ((data->Work == NULL) || (data->iwork == NULL) ){
        printf("error of memory allocation. in Work Exiting ...\n");
        sygvd_parameters_free_2();
        exit(0);
    }
    
  for (auto _ : state) {

    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the sygvd API for benchmarking */
        dsygvd_(&data->itype,&data->jobz,&data->uplo,&data->n,data->A, &data->lda, data->B,
           &data->ldb,data->W,data->Work,&data->lwork,data->iwork,&data->liwork,&data->info); 
        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument dsygvd is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(sygvd_test, dsygvd)->Apply(CustomArguments);