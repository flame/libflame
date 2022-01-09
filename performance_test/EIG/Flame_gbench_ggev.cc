// Gbench headers
#include "benchmark/benchmark.h"
#include "gtest/gtest.h"

// Test suite headers
#include "../Flame_gbench_main.h"
#include "../Flame_gbench_aux.h"

using namespace std;

/* Macros */
#define ggev_parameters_free() \
   if (A!=NULL){ \
      free(A); \
   }\
   if (B!=NULL){ \
      free(B);  \
   } \
   if (Vl!=NULL){ \
      free(Vl);  \
   } \
   if (Vr!=NULL){ \
      free(Vr);  \
   } \
   if (Alpha!=NULL){ \
      free(Alpha);  \
   } \
   if (Beta!=NULL){ \
      free(Beta);  \
   } \
   if (Work!=NULL){ \
      free(Work); \
   }\
   if (Rwork!=NULL){ \
      free(Rwork);  \
   } \

#define ggev_parameters_free_2() \
   if (data->A!=NULL){ \
      free(data->A); \
   }\
   if (data->B!=NULL){ \
      free(data->B);  \
   } \
   if (data->Vl!=NULL){ \
      free(data->Vl);  \
   } \
   if (data->Vr!=NULL){ \
      free(data->Vr);  \
   } \
   if (data->Alpha!=NULL){ \
      free(data->Alpha);  \
   } \
   if (data->Beta!=NULL){ \
      free(data->Beta);  \
   } \
   if (data->Work!=NULL){ \
      free(data->Work); \
   }\
   if (data->Rwork!=NULL){ \
      free(data->Rwork);  \
   } \

/* Begin ggev_parameters  class definition */
template <class T,class S>

class ggev_parameters{
   public:
      char jobvl;
      char jobvr; 
      int  n,lda,ldb,ldvl,ldvr,lwork,info;
      int forced_loop_count;

      /* Local arrays */
      T *A,*B;
      T *Alpha,*Beta,*Vl,*Vr,*Work;
      S *Rwork;
   public:
    ~ggev_parameters ();
    ggev_parameters ( char jobvl_i,char jobvr_i, int n_i)
    {
      jobvl=jobvl_i;
      jobvr=jobvr_i;
      n = n_i; // set test matrix size
      lda = n;
      ldb = n;
      ldvl = n;
      ldvr = n;
      lwork= -1;
      forced_loop_count = 1;

      if (n > FLAME_BIG_MATRIX_SIZE ){
         forced_loop_count = FLAME_GBENCH_FORCED_ITERATION_COUNT;
      }

    /* Memory allocation of the buffers */
         A = (T *)malloc((n*lda)*sizeof(T)) ;
         B = (T *)malloc((n*ldb)*sizeof(T)) ;
         Vl = (T *)malloc((n*ldvl)*sizeof(T)) ;
         Vr = (T *)malloc((n*ldvr)*sizeof(T)) ;
         Alpha = (T *)malloc((n)*sizeof(T)) ;
         Beta = (T *)malloc((n)*sizeof(T)) ;
         Work=NULL;
         Rwork = (S *)malloc((8*n)*sizeof(S)) ;
      if ((A==NULL) || (B==NULL)|| (Vl==NULL)|| (Vr==NULL)|| (Alpha==NULL)|| (Beta==NULL)|| (Rwork==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        ggev_parameters_free();
        exit(0);
      }
    /* Initialization of input Buffers */
      Flame_gbench_init_buffer_rand( A, n*lda);
      Flame_gbench_init_buffer_rand( B, n*ldb);

    } /* end of Constructor  */
};  /* end of ggev_parameters  class definition */



/* Destructor definition  'ggev_parameters' class  */
template <class T,class S>
ggev_parameters<T,S>:: ~ggev_parameters ()
{
#if FLAME_TEST_VERBOSE
   printf(" ggev_parameters object: destructor invoked. \n");
#endif
   ggev_parameters_free ();
}

/* Fixture definition  */
template <class T,class S>
class ggev_test : public ::benchmark::Fixture {
 public:
   int n;
   char jobvl,jobvr;
   std::unique_ptr<ggev_parameters<T,S>> data;

 void SetUp(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      n   = static_cast<size_t>(state.range(0));
      jobvl = (char)static_cast<size_t>(state.range(1));  
      jobvr = (char)static_cast<size_t>(state.range(2));
          
    assert(data.get() == nullptr);
    data.reset(new ggev_parameters<T,S>(jobvl,jobvr,n));

  }

  void TearDown(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      assert(data.get() != nullptr);
      data.reset();
  }

  ~ggev_test() { assert(data == nullptr); }

};

// Customized argument generator for ggev API
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

/* Template Fixture for cggev API  */
BENCHMARK_TEMPLATE_DEFINE_F(ggev_test, cggev,lapack_complex_float,float)(benchmark::State &state) {
  assert(data.get() != nullptr);
      float work_query=-1;
    /* Query optimal working array(s) size */
   cggev_(&data->jobvl, &data->jobvr,&data->n,data->A, &data->lda, data->B,&data->ldb,data->Alpha,data->Beta,data->Vl,
        &data->ldvl,data->Vr,&data->ldvr,(lapack_complex_float *)&work_query,&data->lwork,data->Rwork,&data->info); 

    data->lwork = (int)work_query;
    data->lwork =max(8*(data->n),data->lwork);
    /* Allocate memory for work arrays */
    data->Work = (lapack_complex_float *)malloc( sizeof(lapack_complex_float) * (data->lwork));
    if (data->Work == NULL){
        printf("error of memory allocation. in Work Exiting ...\n");
        ggev_parameters_free_2();
        exit(0);
    }
    
  for (auto _ : state) {

    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the ggev API for benchmarking */
        cggev_(&data->jobvl, &data->jobvr,&data->n,data->A, &data->lda, data->B,&data->ldb,data->Alpha,data->Beta,data->Vl,
        &data->ldvl,data->Vr,&data->ldvr,data->Work,&data->lwork,data->Rwork,&data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument cggev is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(ggev_test, cggev)->Apply(CustomArguments);

/* Template Fixture for zggev API  */
BENCHMARK_TEMPLATE_DEFINE_F(ggev_test, zggev, lapack_complex_double,double)(benchmark::State &state) {
    assert(data.get() != nullptr);

   double work_query=-1;
    /* Query optimal working array(s) size */
    zggev_(&data->jobvl, &data->jobvr,&data->n,data->A, &data->lda,data->B,&data->ldb,data->Alpha,data->Beta,data->Vl,
        &data->ldvl,data->Vr,&data->ldvr,(lapack_complex_double *)&work_query,&data->lwork,data->Rwork,&data->info); 
                                    
    data->lwork = (int)work_query;
    data->lwork =max(8*(data->n),data->lwork);
    /* Allocate memory for work arrays */
    data->Work = (lapack_complex_double *)malloc( sizeof(lapack_complex_double) * (data->lwork));
    if (data->Work==NULL){
        printf("error of memory allocation. Exiting ...\n");
        ggev_parameters_free_2();
        exit(0);
    }
    for (auto _ : state) {
        
        for (int i=0; i<data->forced_loop_count; i++) {
           /* Invoke the ggev API for benchmarking */
            zggev_(&data->jobvl, &data->jobvr,&data->n,data->A, &data->lda, data->B,&data->ldb,data->Alpha,data->Beta,data->Vl,
            &data->ldvl,data->Vr,&data->ldvr,data->Work,&data->lwork,data->Rwork,&data->info);        
            if( data->info < 0 ) {
                printf( "\n warning: The %d th argument zggev is wrong\n", data->info );
            }
        }
    }
}

BENCHMARK_REGISTER_F(ggev_test, zggev)->Apply(CustomArguments);