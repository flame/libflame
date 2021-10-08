// Gbench headers
#include "benchmark/benchmark.h"
#include "gtest/gtest.h"

// Test suite headers
#include "../Flame_gbench_main.h"
#include "../Flame_gbench_aux.h"

using namespace std;

/* Macros */
#define gesdd_parameters_free() \
   if (A!=NULL){ \
      free(A); \
   }\
   if (Se!=NULL){ \
      free(Se);  \
   } \
   if (U!=NULL){ \
      free(U);  \
   } \
   if (Vt!=NULL){ \
      free(Vt);  \
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

#define gesdd_parameters_free_2() \
   if (data->A!=NULL){ \
      free(data->A); \
   }\
   if (data->Se!=NULL){ \
      free(data->Se);  \
   } \
   if (data->U!=NULL){ \
      free(data->U);  \
   } \
   if (data->Vt!=NULL){ \
      free(data->Vt);  \
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
   


/* Begin gesdd_parameters  class definition */
template <class T,class S>

class gesdd_parameters{
   public:
      char jobz;
      int  m,n,lda,lwork,ldvt,ldu,info;
      int forced_loop_count;

      /* Local arrays */
      T *A,*Work,*U,*Vt;
      int *iwork;
      S *Rwork,*Se;
   public:
    ~gesdd_parameters ();
    gesdd_parameters ( char jobz_i,int m_i, int n_i)
    {
      jobz=jobz_i;
      m = m_i;
      n = n_i; // set test matrix size
      int mx=max(m,n);
      lda = mx;
      ldu = mx;
      ldvt = mx; 
      lwork= -1;
      forced_loop_count = 1;

      if (n > FLAME_BIG_MATRIX_SIZE ){
         forced_loop_count = FLAME_GBENCH_FORCED_ITERATION_COUNT;
      }

    /* Memory allocation of the buffers */
         A = (T *)malloc((mx*lda)*sizeof(T)) ;
         U = (T *)malloc((mx*ldu)*sizeof(T)) ;
         Vt = (T *)malloc((mx*ldvt)*sizeof(T)) ;
         Se = (S *)malloc((mx)*sizeof(S)) ;
         Work=NULL;
         Rwork=NULL;
         iwork= (int *)malloc((8*mx)*sizeof(int)) ;
      if ((A==NULL) || (U==NULL) || (Vt==NULL)|| (Se==NULL) || (iwork==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        gesdd_parameters_free();
        exit(0);
      }
    /* Initialization of input Buffers */
      Flame_gbench_init_buffer_rand( A, n*lda);

    } /* end of Constructor  */
};  /* end of gesdd_parameters  class definition */



/* Destructor definition  'gesdd_parameters' class  */
template <class T,class S>
gesdd_parameters<T,S>:: ~gesdd_parameters ()
{
#if FLAME_TEST_VERBOSE
   printf(" gesdd_parameters object: destructor invoked. \n");
#endif
     gesdd_parameters_free ();
}

/* Fixture definition  */
template <class T,class S>
class gesdd_test : public ::benchmark::Fixture {
 public:
   int m,n;
   char jobz;
   std::unique_ptr<gesdd_parameters<T,S>> data;

 void SetUp(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      n   = static_cast<size_t>(state.range(0));
      m   = static_cast<size_t>(state.range(1));
      jobz = (char)static_cast<size_t>(state.range(2));  
          
    assert(data.get() == nullptr);
    data.reset(new gesdd_parameters<T,S>(jobz,m,n));

  }

  void TearDown(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      assert(data.get() != nullptr);
      data.reset();
  }

  ~gesdd_test() { assert(data == nullptr); }

};

// Customized argument generator for gesdd API
static void CustomArguments(benchmark::internal::Benchmark* b) {
    
    if (svd_paramslist[0].mode == MODE_RANGE) {
       for (int i = 0; i < svd_paramslist[0].num_tests ; i++){
           b->ArgsProduct({
              {benchmark::CreateDenseRange( svd_paramslist[i].n_range_start,
                                            svd_paramslist[i].n_range_end,
                                            svd_paramslist[i].n_range_step_size)},
              {benchmark::CreateDenseRange( svd_paramslist[i].m_range_start,
                                            svd_paramslist[i].m_range_end,
                                            svd_paramslist[i].m_range_step_size)},
                                            {'N','O','S','A'}
                                            
                                          });
       }
    }
   
    
    if (svd_paramslist[0].mode == MODE_DISCRETE) {
       for (int i = 0; i < svd_paramslist[0].num_tests ; i++)
       {
           b->Args({svd_paramslist[i].n, svd_paramslist[i].m,svd_paramslist[i].jobu_gesvd});
       }
    }
    if (svd_paramslist[0].mode == MODE_COMBINATIONAL) {
       for (int i = 0; i < svd_paramslist[0].num_tests ; i++)
       {
          for (int j = 0; j < svd_paramslist[0].num_tests ; j++){
              for(int k=0; k< svd_paramslist[0].num_tests ; k++){
            b->Args({svd_paramslist[i].n, svd_paramslist[j].m,svd_paramslist[k].jobu_gesvd});
              }
          }
       }
    }
}

/* Template Fixture for sgesdd API  */
BENCHMARK_TEMPLATE_DEFINE_F(gesdd_test, sgesdd,float,float )(benchmark::State &state) {
  assert(data.get() != nullptr);
      float work_query=-1;
      int mx=max(data->m,data->n);
      int mn=min(data->m,data->n);
      
    /* Query optimal working array(s) size */
   sgesdd_(&data->jobz,&data->m,&data->n,data->A, &data->lda,data->Se,data->U,&data->ldu,data->Vt,&data->ldvt,
   (float *) &work_query, &data->lwork,data->iwork, &data->info);

    data->lwork = (int)work_query;
    data->lwork=max((2*mn*mn+2*mn+mx),data->lwork);
        
    /* Allocate memory for work arrays */
    data->Work = (float *)malloc( sizeof(float) * (data->lwork));
    if ((data->Work==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        gesdd_parameters_free_2();
        exit(0);
    }
    

  for (auto _ : state) {
    
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the gesdd API for benchmarking */
        sgesdd_(&data->jobz,&data->m,&data->n,data->A,&data->lda,data->Se,data->U,&data->ldu,data->Vt,&data->ldvt,
                data->Work, &data->lwork,data->iwork, &data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument sgesdd is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(gesdd_test,sgesdd)->Apply(CustomArguments);

/* Template Fixture for dgesdd API  */
BENCHMARK_TEMPLATE_DEFINE_F(gesdd_test, dgesdd,double,double )(benchmark::State &state) {
  assert(data.get() != nullptr);
      double work_query=-1;
      int mx=max(data->m,data->n);
      int mn=min(data->m,data->n);
      
    /* Query optimal working array(s) size */
   dgesdd_(&data->jobz,&data->m,&data->n,data->A, &data->lda,data->Se,data->U,&data->ldu,data->Vt,&data->ldvt,
   (double *) &work_query, &data->lwork,data->iwork, &data->info);

    data->lwork = (int)work_query;
    data->lwork=max((2*mn*mn+2*mn+mx),data->lwork);
        
    /* Allocate memory for work arrays */
    data->Work = (double *)malloc( sizeof(double) * (data->lwork));
    if ((data->Work==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        gesdd_parameters_free_2();
        exit(0);
    }
    

  for (auto _ : state) {
    
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the gesdd API for benchmarking */
        dgesdd_(&data->jobz,&data->m,&data->n,data->A,&data->lda,data->Se,data->U,&data->ldu,data->Vt,&data->ldvt,
                data->Work, &data->lwork,data->iwork, &data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument dgesdd is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(gesdd_test,dgesdd)->Apply(CustomArguments);

/* Template Fixture for cgesdd API  */
BENCHMARK_TEMPLATE_DEFINE_F(gesdd_test, cgesdd,lapack_complex_float,float )(benchmark::State &state) {
  assert(data.get() != nullptr);
      float work_query=-1;
      int mx=max(data->m,data->n);
      int mn=min(data->m,data->n);
      int lrwork= max((5*(mn*mn)+5*mn),((2*mx*mn)+(2*mn*mn)+mn));
      data->Rwork = (float *)malloc((lrwork)*sizeof(float)) ;

      if ((data->Rwork==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        gesdd_parameters_free_2();
        exit(0);
    }

    /* Query optimal working array(s) size */
   cgesdd_(&data->jobz,&data->m,&data->n,data->A, &data->lda,data->Se,data->U,&data->ldu,data->Vt,&data->ldvt,
   (lapack_complex_float *) &work_query, &data->lwork,data->Rwork,data->iwork, &data->info);

    data->lwork = (int)work_query;
    data->lwork=max((2*mn*mn+2*mn+mx),data->lwork);
        
    /* Allocate memory for work arrays */
    data->Work = (lapack_complex_float *)malloc( sizeof(lapack_complex_float) * (data->lwork));
    if ((data->Work==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        gesdd_parameters_free_2();
        exit(0);
    }
    

  for (auto _ : state) {
    
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the gesdd API for benchmarking */
        cgesdd_(&data->jobz,&data->m,&data->n,data->A,&data->lda,data->Se,data->U,&data->ldu,data->Vt,&data->ldvt,
                data->Work, &data->lwork,data->Rwork,data->iwork, &data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument cgesdd is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(gesdd_test,cgesdd)->Apply(CustomArguments);

/* Template Fixture for zgesdd API  */
BENCHMARK_TEMPLATE_DEFINE_F(gesdd_test, zgesdd,lapack_complex_double,double )(benchmark::State &state) {
  assert(data.get() != nullptr);
      double work_query=-1;
      int mx=max(data->m,data->n);
      int mn=min(data->m,data->n);
      int lrwork= max((5*(mn*mn)+5*mn),((2*mx*mn)+(2*mn*mn)+mn));
      data->Rwork = (double *)malloc((lrwork)*sizeof(double)) ;

      if ((data->Rwork==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        gesdd_parameters_free_2();
        exit(0);
    }
    /* Query optimal working array(s) size */
   zgesdd_(&data->jobz,&data->m,&data->n,data->A, &data->lda,data->Se,data->U,&data->ldu,data->Vt,&data->ldvt,
   (lapack_complex_double *) &work_query, &data->lwork,data->Rwork,data->iwork, &data->info);

    data->lwork = (int)work_query;

    
    data->lwork=max((2*mn*mn+2*mn+mx),data->lwork);
        
    /* Allocate memory for work arrays */
    data->Work = (lapack_complex_double *)malloc( sizeof(lapack_complex_double) * (data->lwork));
    if ((data->Work==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        gesdd_parameters_free_2();
        exit(0);
    }
    

  for (auto _ : state) {
    
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the gesdd API for benchmarking */
        zgesdd_(&data->jobz,&data->m,&data->n,data->A,&data->lda,data->Se,data->U,&data->ldu,data->Vt,&data->ldvt,
                data->Work, &data->lwork,data->Rwork,data->iwork, &data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument zgesdd is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(gesdd_test,zgesdd)->Apply(CustomArguments);