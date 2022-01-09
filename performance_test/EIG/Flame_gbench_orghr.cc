// Gbench headers
#include "benchmark/benchmark.h"
#include "gtest/gtest.h"

// Test suite headers
#include "../Flame_gbench_main.h"
#include "../Flame_gbench_aux.h"

using namespace std;

/* Macros */
#define orghr_parameters_free() \
   if (A!=NULL){ \
      free(A); \
   }\
   if (Tau!=NULL){ \
      free(Tau);  \
   } \
   if (Work!=NULL){ \
      free(Work); \
   }\
   
#define orghr_parameters_free_2() \
   if (data->A!=NULL){ \
      free(data->A); \
   }\
   if (data->Tau!=NULL){ \
      free(data->Tau);  \
   } \
   if (data->Work!=NULL){ \
      free(data->Work); \
   }\
   
/* Begin orghr_parameters  class definition */
template <class T>

class orghr_parameters{
   public: 
      int  n,lda,ilo,ihi,lwork,info;
      int forced_loop_count;

      /* Local arrays */
      T *A;
      T *Tau,*Work;
   public:
    ~orghr_parameters ();
    orghr_parameters ( int n_i)
    {
      n = n_i; // set test matrix size
      lda = n;
      ilo = max(1,n);
      ihi = max(1,n);
      lwork= -1;
      forced_loop_count = 1;

      if (n > FLAME_BIG_MATRIX_SIZE ){
         forced_loop_count = FLAME_GBENCH_FORCED_ITERATION_COUNT;
      }

    /* Memory allocation of the buffers */
         A = (T *)malloc((n*lda)*sizeof(T)) ;
         Tau = (T *)malloc((n)*sizeof(T)) ;
         Work=NULL;
      if ((A==NULL) || (Tau==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        orghr_parameters_free();
        exit(0);
      }
    /* Initialization of input Buffers */
      Flame_gbench_init_buffer_rand( A, n*lda);
    } /* end of Constructor  */
};  /* end of orghr_parameters  class definition */



/* Destructor definition  'orghr_parameters' class  */
template <class T>
orghr_parameters<T>:: ~orghr_parameters ()
{
#if FLAME_TEST_VERBOSE
   printf(" orghr_parameters object: destructor invoked. \n");
#endif
   orghr_parameters_free ();
}

/* Fixture definition  */
template <class T>
class orghr_test : public ::benchmark::Fixture {
 public:
   int n;
   std::unique_ptr<orghr_parameters<T>> data;

 void SetUp(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      n   = static_cast<size_t>(state.range(0));          
    assert(data.get() == nullptr);
    data.reset(new orghr_parameters<T>(n));

  }

  void TearDown(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      assert(data.get() != nullptr);
      data.reset();
  }

  ~orghr_test() { assert(data == nullptr); }

};

// Customized argument generator for orghr API
static void CustomArguments(benchmark::internal::Benchmark* b) {
    
    if (eig_non_sym_paramslist[0].mode == MODE_RANGE) {
       for (int i = 0; i < eig_non_sym_paramslist[0].num_tests ; i++){
           b->ArgsProduct(
              {benchmark::CreateDenseRange( eig_non_sym_paramslist[i].n_range_start,
                                            eig_non_sym_paramslist[i].n_range_end,
                                            eig_non_sym_paramslist[i].n_range_step_size)}
              );
       }
    }
     
    if (eig_non_sym_paramslist[0].mode == MODE_DISCRETE) {
       for (int i = 0; i < eig_non_sym_paramslist[0].num_tests ; i++)
       {
           b->Args({eig_non_sym_paramslist[i].n});
       }
    }

    if (eig_non_sym_paramslist[0].mode == MODE_COMBINATIONAL) {
       for (int i = 0; i < eig_non_sym_paramslist[0].num_tests ; i++)
       {
            b->Args({eig_non_sym_paramslist[i].n});
              
       }
    }
}
/* Template Fixture for dorghr API  */
BENCHMARK_TEMPLATE_DEFINE_F(orghr_test, dorghr,double)(benchmark::State &state) {
  assert(data.get() != nullptr);
      double work_query=-1;
    /* Query optimal working array(s) size */
   dorghr_(&data->n,&data->ilo,&data->ihi,data->A, &data->lda,data->Tau,&work_query,&data->lwork,&data->info); 

    data->lwork = (int)work_query;
    data->lwork =max(n,data->lwork);
    /* Allocate memory for work arrays */
    data->Work = (double *)malloc( sizeof(double) * (data->lwork));
    if (data->Work == NULL){
        printf("error of memory allocation. in Work Exiting ...\n");
        orghr_parameters_free_2();
        exit(0);
    }    
  for (auto _ : state) {
     
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the orghr API for benchmarking */
        dorghr_(&data->n,&data->ilo,&data->ihi,data->A, &data->lda,data->Tau,data->Work,&data->lwork,&data->info); 
        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument dorghr is wrong\n", data->info );
        }
    }
  }
  
}

BENCHMARK_REGISTER_F(orghr_test, dorghr)->Apply(CustomArguments);

/* Template Fixture for sorghr API  */
BENCHMARK_TEMPLATE_DEFINE_F(orghr_test, sorghr,float)(benchmark::State &state) {
  assert(data.get() != nullptr);
      float work_query=-1;
    /* Query optimal working array(s) size */
   sorghr_(&data->n,&data->ilo,&data->ihi,data->A, &data->lda,data->Tau,&work_query,&data->lwork,&data->info); 

    data->lwork = (int)work_query;
    data->lwork =max(n,data->lwork);
    /* Allocate memory for work arrays */
    data->Work = (float *)malloc( sizeof(float) * (data->lwork));
    if (data->Work == NULL){
        printf("error of memory allocation. in Work Exiting ...\n");
        orghr_parameters_free_2();
        exit(0);
    }    
  for (auto _ : state) {
     
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the orghr API for benchmarking */
        sorghr_(&data->n,&data->ilo,&data->ihi,data->A, &data->lda,data->Tau,data->Work,&data->lwork,&data->info); 
        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument sorghr is wrong\n", data->info );
        }
    }
  }
  
}

BENCHMARK_REGISTER_F(orghr_test, sorghr)->Apply(CustomArguments);