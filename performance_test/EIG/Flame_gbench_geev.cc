// Gbench headers
#include "benchmark/benchmark.h"
#include "gtest/gtest.h"

// Test suite headers
#include "../Flame_gbench_main.h"
#include "../Flame_gbench_aux.h"

using namespace std;

/* Macros */
#define geev_parameters_free() \
   if (A!=NULL){ \
      free(A); \
   }\
   if (Vl!=NULL){ \
      free(Vl);  \
   } \
   if (Vr!=NULL){ \
      free(Vr);  \
   } \
   if (W!=NULL){ \
      free(W);  \
   } \
   if (Work!=NULL){ \
      free(Work); \
   }\
   if (Rwork!=NULL){ \
      free(Rwork);  \
   } \

#define geev_parameters_free_2() \
   if (data->A!=NULL){ \
      free(data->A); \
   }\
   if (data->Vl!=NULL){ \
      free(data->Vl);  \
   } \
   if (data->Vr!=NULL){ \
      free(data->Vr);  \
   } \
   if (data->W!=NULL){ \
      free(data->W);  \
   } \
   if (data->Work!=NULL){ \
      free(data->Work); \
   }\
   if (data->Rwork!=NULL){ \
      free(data->Rwork);  \
   } \

/* Begin geev_parameters  class definition */
template <class T,class S>

class geev_parameters{
   public:
      char jobvl;
      char jobvr; 
      int  n,lda,ldvl,ldvr,lwork,info;
      int forced_loop_count;

      /* Local arrays */
      T *A;
      T *W,*Vl,*Vr,*Work;
      S *Rwork; // NOTE: utilized as 'WI' for float, double variants
     // int *iwork;
   public:
    ~geev_parameters ();
    geev_parameters ( char jobvl_i,char jobvr_i, int n_i)
    {
      jobvl=jobvl_i;
      jobvr=jobvr_i;
      n = n_i; // set test matrix size
      lda = n;
      ldvl = n;
      ldvr = n;
      lwork= -1;
      forced_loop_count = 1;

      if (n > FLAME_BIG_MATRIX_SIZE ){
         forced_loop_count = FLAME_GBENCH_FORCED_ITERATION_COUNT;
      }

    /* Memory allocation of the buffers */
         A = (T *)malloc((n*lda)*sizeof(T)) ;
         Vl = (T *)malloc((n*ldvl)*sizeof(T)) ;
         Vr = (T *)malloc((n*ldvr)*sizeof(T)) ;
         W = (T *)malloc((n)*sizeof(T)) ;
        // Work = (T *)malloc((n)*sizeof(T)) ;
         Work=NULL;
         Rwork = (S *)malloc((2*n)*sizeof(S)) ;
      if ((A==NULL) || (Vl==NULL)|| (Vr==NULL)|| (W==NULL)|| (Rwork==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        geev_parameters_free();
        exit(0);
      }
    /* Initialization of input Buffers */
      Flame_gbench_init_buffer_rand( A, n*lda);
    } /* end of Constructor  */
};  /* end of geev_parameters  class definition */



/* Destructor definition  'geev_parameters' class  */
template <class T,class S>
geev_parameters<T,S>:: ~geev_parameters ()
{
#if FLAME_TEST_VERBOSE
   printf(" geev_parameters object: destructor invoked. \n");
#endif
   geev_parameters_free ();
}

/* Fixture definition  */
template <class T,class S>
class geev_test : public ::benchmark::Fixture {
 public:
   int n;
   char jobvl,jobvr;
   std::unique_ptr<geev_parameters<T,S>> data;

 void SetUp(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      n   = static_cast<size_t>(state.range(0));
      jobvl = (char)static_cast<size_t>(state.range(1));  
      jobvr = (char)static_cast<size_t>(state.range(2));
          
    assert(data.get() == nullptr);
    data.reset(new geev_parameters<T,S>(jobvl,jobvr,n));

  }

  void TearDown(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      assert(data.get() != nullptr);
      data.reset();
  }

  ~geev_test() { assert(data == nullptr); }

};

// Customized argument generator for geev API
static void CustomArguments(benchmark::internal::Benchmark* b) {
    
    if (eig_non_sym_paramslist[0].mode == MODE_RANGE) {
       for (int i = 0; i < eig_non_sym_paramslist[0].num_tests ; i++){
           b->ArgsProduct({
              {benchmark::CreateDenseRange( eig_non_sym_paramslist[i].n_range_start,
                                            eig_non_sym_paramslist[i].n_range_end,
                                            eig_non_sym_paramslist[i].n_range_step_size)},
                                            {'N', 'V'},
                                            {'N', 'V'}
              });
       }
    }
    
    if (eig_non_sym_paramslist[0].mode == MODE_DISCRETE) {
       for (int i = 0; i < eig_non_sym_paramslist[0].num_tests ; i++)
       {
           b->Args({eig_non_sym_paramslist[i].n, eig_non_sym_paramslist[i].jobvsl, eig_non_sym_paramslist[i].jobvsr});
       }
    }
    if (eig_non_sym_paramslist[0].mode == MODE_COMBINATIONAL) {
       for (int i = 0; i < eig_non_sym_paramslist[0].num_tests ; i++)
       {
          for (int j = 0; j < eig_non_sym_paramslist[0].num_tests ; j++){
              for(int k=0; k< eig_non_sym_paramslist[0].num_tests ; k++){
            b->Args({eig_non_sym_paramslist[i].n, eig_non_sym_paramslist[j].jobvsl,eig_non_sym_paramslist[k].jobvsr});
              }
          }
       }
    }
}

/* Template Fixture for sgeev API  */
BENCHMARK_TEMPLATE_DEFINE_F(geev_test, sgeev, float,float )(benchmark::State &state) {
  assert(data.get() != nullptr);
      //int  iwork_query=-1;
      float work_query=-1;
    /* Query for optimal working array(s) size */
   sgeev_( &data->jobvl, &data->jobvr,&data->n,data->A, &data->lda,
           data->W,data->Rwork, data->Vl, &data->ldvl, data->Vr,
		   &data->ldvr, (float *)&work_query, &data->lwork, &data->info); 

    data->lwork = (int)work_query;
    data->lwork = max(1, data->lwork);
    /* Allocate memory for work arrays */
    //data->lwork=2*(data->n);
    data->Work = (float *)malloc( sizeof(float) * (data->lwork));
    if (data->Work == NULL){
        printf("error of memory allocation. in Work Exiting ...\n");
        geev_parameters_free_2();
        exit(0);
    }

  for (auto _ : state) {

    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the geev API for benchmarking */
        sgeev_( &data->jobvl, &data->jobvr, &data->n, data->A,
		        &data->lda, data->W, data->Rwork, data->Vl,
                &data->ldvl,data->Vr,&data->ldvr,data->Work,
				&data->lwork,&data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument sgeev is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(geev_test, sgeev)->Apply(CustomArguments);

/* Template Fixture for dgeev API  */
BENCHMARK_TEMPLATE_DEFINE_F(geev_test, dgeev, double,double )(benchmark::State &state) {
   assert(data.get() != nullptr);

   double work_query=-1;
   /* Query for optimal working array(s) size */
   dgeev_( &data->jobvl, &data->jobvr,&data->n,data->A, &data->lda,
           data->W,data->Rwork, data->Vl, &data->ldvl, data->Vr,
		   &data->ldvr, (double *)&work_query, &data->lwork, &data->info); 

    data->lwork = (int)work_query;
    data->lwork = max(1, data->lwork);

    /* Allocate memory for work arrays */
    data->Work = (double *)malloc( sizeof(double) * (data->lwork));
    if (data->Work == NULL){
        printf("error of memory allocation. in Work Exiting ...\n");
        geev_parameters_free_2();
        exit(0);
    }

  for (auto _ : state) {

    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the geev API for benchmarking */
        dgeev_( &data->jobvl, &data->jobvr, &data->n, data->A,
		        &data->lda, data->W, data->Rwork, data->Vl,
                &data->ldvl,data->Vr,&data->ldvr,data->Work,
				&data->lwork,&data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument dgeev is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(geev_test, dgeev)->Apply(CustomArguments);
/* Template Fixture for cgeev API  */
BENCHMARK_TEMPLATE_DEFINE_F(geev_test, cgeev, float _Complex, float )(benchmark::State &state) {
  assert(data.get() != nullptr);
      //int  iwork_query=-1;
      float work_query=-1;
    /* Query for optimal working array(s) size */
   cgeev_(&data->jobvl, &data->jobvr,&data->n,data->A, &data->lda, data->W,data->Vl,
        &data->ldvl,data->Vr,&data->ldvr,(float _Complex *)&work_query,&data->lwork,data->Rwork,&data->info); 

    data->lwork = (int)work_query;
    data->lwork = max(1, data->lwork);
    /* Allocate memory for work arrays */
    //data->lwork=2*(data->n);
    data->Work = (float _Complex *)malloc( sizeof(float _Complex) * (data->lwork));
    if (data->Work == NULL){
        printf("error of memory allocation. in Work Exiting ...\n");
        geev_parameters_free_2();
        exit(0);
    }
    

  for (auto _ : state) {

    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the geev API for benchmarking */
        cgeev_(&data->jobvl, &data->jobvr,&data->n,data->A, &data->lda, data->W,data->Vl,
        &data->ldvl,data->Vr,&data->ldvr,data->Work,&data->lwork,data->Rwork,&data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument cgeev is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(geev_test, cgeev)->Apply(CustomArguments);

/* Template Fixture for zgeev API  */
BENCHMARK_TEMPLATE_DEFINE_F(geev_test, zgeev, double _Complex, double)(benchmark::State &state) {
    assert(data.get() != nullptr);

   double work_query=-1;
    /* Query for optimal working array(s) size */
    zgeev_(&data->jobvl, &data->jobvr,&data->n,data->A, &data->lda, data->W,data->Vl,
        &data->ldvl,data->Vr,&data->ldvr,(double _Complex *)&work_query,&data->lwork,data->Rwork,&data->info); 
                                    
    data->lwork = (int)work_query;
    data->lwork = max(1, data->lwork);
    //data->lwork=2*(data->n);
    /* Allocate memory for work arrays */
    data->Work = (double _Complex *)malloc( sizeof(double _Complex) * (data->lwork));
    if ((data->Work==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        geev_parameters_free_2();
        exit(0);
    }
    for (auto _ : state) {
        
        for (int i=0; i<data->forced_loop_count; i++) {
           /* Invoke the geev API for benchmarking */
            zgeev_(&data->jobvl, &data->jobvr,&data->n,data->A, &data->lda, data->W,data->Vl,
            &data->ldvl,data->Vr,&data->ldvr,data->Work,&data->lwork,data->Rwork,&data->info);        
            if( data->info < 0 ) {
                printf( "\n warning: The %d th argument zgeev is wrong\n", data->info );
            }
        }
    }
}

BENCHMARK_REGISTER_F(geev_test, zgeev)->Apply(CustomArguments);
