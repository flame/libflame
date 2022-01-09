// Gbench headers
#include "benchmark/benchmark.h"
#include "gtest/gtest.h"

// Test suite headers
#include "../Flame_gbench_main.h"
#include "../Flame_gbench_aux.h"

using namespace std;

/* Macros */
#define hseqr_parameters_free() \
   if (H!=NULL){ \
      free(H); \
   }\
   if (Z!=NULL){ \
      free(Z); \
   }\
   if (W!=NULL){ \
      free(W);  \
   } \
   if (Wi!=NULL){ \
      free(Wi);  \
   } \
   if (Work!=NULL){ \
      free(Work); \
   }\
   
#define hseqr_parameters_free_2() \
   if (data->H!=NULL){ \
      free(data->H); \
   }\
   if (data->Z!=NULL){ \
      free(data->Z); \
   }\
   if (data->W!=NULL){ \
      free(data->W);  \
   } \
   if (data->Wi!=NULL){ \
      free(data->Wi);  \
   } \
   if (data->Work!=NULL){ \
      free(data->Work); \
   }\
   
/* Begin hseqr_parameters  class definition */
template <class T>

class hseqr_parameters{
   public: 
      int  n,ldh,ldz,ilo,ihi,lwork,info;
      int forced_loop_count;
      char job,compz;

      /* Local arrays */
      T *H,*Z;
      T *W,*Work,*Wi;
   public:
    ~hseqr_parameters ();
    hseqr_parameters ( int n_i,char job_i,char compz_i)
    {
      n = n_i; // set test matrix size
      job = job_i;
      compz = compz_i;
      ldh = n;
      ldz = n;
      ilo = max(1,n);
      ihi = max(1,n);
      lwork= -1;
      forced_loop_count = 1;

      if (n > FLAME_BIG_MATRIX_SIZE ){
         forced_loop_count = FLAME_GBENCH_FORCED_ITERATION_COUNT;
      }

    /* Memory allocation of the buffers */
         H = (T *)malloc((n*ldh)*sizeof(T)) ;
         Z = (T *)malloc((n*ldz)*sizeof(T)) ;
         W = (T *)malloc((n)*sizeof(T)) ;
         Work=NULL;
         Wi = NULL;
      if ((H==NULL) || (Z==NULL) || (W==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        hseqr_parameters_free();
        exit(0);
      }
    /* Initialization of input Buffers */
      Flame_gbench_init_buffer_rand( H, n*ldh);
      Flame_gbench_init_buffer_rand( Z, n*ldz);
    } /* end of Constructor  */
};  /* end of hseqr_parameters  class definition */



/* Destructor definition  'hseqr_parameters' class  */
template <class T>
hseqr_parameters<T>:: ~hseqr_parameters ()
{
#if FLAME_TEST_VERBOSE
   printf(" hseqr_parameters object: destructor invoked. \n");
#endif
   hseqr_parameters_free ();
}

/* Fixture definition  */
template <class T>
class hseqr_test : public ::benchmark::Fixture {
 public:
   int n;
   char job,compz;
   std::unique_ptr<hseqr_parameters<T>> data;

 void SetUp(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      n   = static_cast<size_t>(state.range(0));
      job = (char) static_cast<size_t>(state.range(1));
      compz = (char) static_cast<size_t>(state.range(2));          
    assert(data.get() == nullptr);
    data.reset(new hseqr_parameters<T>(n,job,compz));

  }

  void TearDown(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      assert(data.get() != nullptr);
      data.reset();
  }

  ~hseqr_test() { assert(data == nullptr); }

};

// Customized argument generator for hseqr API
static void CustomArguments(benchmark::internal::Benchmark* b) {
    
    if (eig_paramslist[0].mode == MODE_RANGE) {
       for (int i = 0; i < eig_paramslist[0].num_tests ; i++){
           b->ArgsProduct({
              {benchmark::CreateDenseRange( eig_paramslist[i].n_range_start,
                                            eig_paramslist[i].n_range_end,
                                            eig_paramslist[i].n_range_step_size)},
                                            {'E','S'},
                                            {'N','I','V'}
           });
       }
    }
     
    if (eig_paramslist[0].mode == MODE_DISCRETE) {
       for (int i = 0; i < eig_paramslist[0].num_tests ; i++)
       {
           b->Args({{eig_paramslist[i].n},{eig_paramslist[i].job_seqr},{eig_paramslist[i].compz}});
       }
    }

    if (eig_paramslist[0].mode == MODE_COMBINATIONAL) {
       for (int i = 0; i < eig_paramslist[0].num_tests ; i++)
       {
           for (int j = 0; j < eig_paramslist[0].num_tests ; j++)
           {
               for (int k = 0; k < eig_paramslist[0].num_tests ; k++)
               {
                    b->Args({{eig_paramslist[i].n},{eig_paramslist[j].job_seqr},
                    {eig_paramslist[k].compz}}  );
               }
           }
              
       }
    }
}
/* Template Fixture for shseqr API  */
BENCHMARK_TEMPLATE_DEFINE_F(hseqr_test, shseqr,float)(benchmark::State &state) {
  assert(data.get() != nullptr);
      float work_query=-1;
      data->Wi = (float *)malloc( sizeof(float) * (data->n));
      if (data->Wi == NULL){
        printf("error of memory allocation. in Work Exiting ...\n");
        hseqr_parameters_free_2();
        exit(0);
    }    
    /* Query optimal working array(s) size */
   shseqr_(&data->job,&data->compz,&data->n,&data->ilo,&data->ihi,data->H, &data->ldh,data->W,data->Wi,
   data->Z,&data->ldz,(float *)&work_query,&data->lwork,&data->info); 

    data->lwork = (int)work_query;
    data->lwork =max(n,data->lwork);
    /* Allocate memory for work arrays */
    data->Work = (float *)malloc( sizeof(float) * (data->lwork));
    if (data->Work == NULL){
        printf("error of memory allocation. in Work Exiting ...\n");
        hseqr_parameters_free_2();
        exit(0);
    }    
  for (auto _ : state) {
     
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the hseqr API for benchmarking */
        shseqr_(&data->job,&data->compz,&data->n,&data->ilo,&data->ihi,data->H, &data->ldh,data->W,data->Wi,
                data->Z,&data->ldz,data->Work,&data->lwork,&data->info); 
        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument shseqr is wrong\n", data->info );
        }
    }
  }
  
}

BENCHMARK_REGISTER_F(hseqr_test, shseqr)->Apply(CustomArguments);

/* Template Fixture for dhseqr API  */
BENCHMARK_TEMPLATE_DEFINE_F(hseqr_test, dhseqr,double)(benchmark::State &state) {
  assert(data.get() != nullptr);
      double work_query=-1;
      data->Wi = (double *)malloc( sizeof(double) * (data->n));
      if (data->Wi == NULL){
        printf("error of memory allocation. in Work Exiting ...\n");
        hseqr_parameters_free_2();
        exit(0);
    }    

    /* Query optimal working array(s) size */
   dhseqr_(&data->job,&data->compz,&data->n,&data->ilo,&data->ihi,data->H, &data->ldh,data->W,data->Wi,
   data->Z,&data->ldz,(double *)&work_query,&data->lwork,&data->info); 

    data->lwork = (int)work_query;
    data->lwork =max(n,data->lwork);
    /* Allocate memory for work arrays */
    data->Work = (double *)malloc( sizeof(double) * (data->lwork));
    if (data->Work == NULL){
        printf("error of memory allocation. in Work Exiting ...\n");
        hseqr_parameters_free_2();
        exit(0);
    }    
  for (auto _ : state) {
     
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the hseqr API for benchmarking */
        dhseqr_(&data->job,&data->compz,&data->n,&data->ilo,&data->ihi,data->H, &data->ldh,data->W,data->Wi,
                data->Z,&data->ldz,data->Work,&data->lwork,&data->info); 
        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument dhseqr is wrong\n", data->info );
        }
    }
  }
  
}

BENCHMARK_REGISTER_F(hseqr_test, dhseqr)->Apply(CustomArguments);

/* Template Fixture for chseqr API  */
BENCHMARK_TEMPLATE_DEFINE_F(hseqr_test, chseqr,lapack_complex_float)(benchmark::State &state) {
  assert(data.get() != nullptr);
      float work_query=-1;
    /* Query optimal working array(s) size */
   chseqr_(&data->job,&data->compz,&data->n,&data->ilo,&data->ihi,data->H, &data->ldh,data->W,
   data->Z,&data->ldz,(lapack_complex_float *)&work_query,&data->lwork,&data->info); 

    data->lwork = (int)work_query;
    data->lwork =max(n,data->lwork);
    /* Allocate memory for work arrays */
    data->Work = (lapack_complex_float *)malloc( sizeof(lapack_complex_float) * (data->lwork));
    if (data->Work == NULL){
        printf("error of memory allocation. in Work Exiting ...\n");
        hseqr_parameters_free_2();
        exit(0);
    }    
  for (auto _ : state) {
     
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the hseqr API for benchmarking */
        chseqr_(&data->job,&data->compz,&data->n,&data->ilo,&data->ihi,data->H, &data->ldh,data->W,
                data->Z,&data->ldz,data->Work,&data->lwork,&data->info); 
        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument chseqr is wrong\n", data->info );
        }
    }
  }
  
}

BENCHMARK_REGISTER_F(hseqr_test, chseqr)->Apply(CustomArguments);

/* Template Fixture for zhseqr API  */
BENCHMARK_TEMPLATE_DEFINE_F(hseqr_test, zhseqr,lapack_complex_double)(benchmark::State &state) {
  assert(data.get() != nullptr);
      double work_query=-1;
    /* Query optimal working array(s) size */
   zhseqr_(&data->job,&data->compz,&data->n,&data->ilo,&data->ihi,data->H, &data->ldh,data->W,
   data->Z,&data->ldz,(lapack_complex_double *)&work_query,&data->lwork,&data->info); 

    data->lwork = (int)work_query;
    data->lwork =max(n,data->lwork);
    /* Allocate memory for work arrays */
    data->Work = (lapack_complex_double *)malloc( sizeof(lapack_complex_double) * (data->lwork));
    if (data->Work == NULL){
        printf("error of memory allocation. in Work Exiting ...\n");
        hseqr_parameters_free_2();
        exit(0);
    }    
  for (auto _ : state) {
     
    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the hseqr API for benchmarking */
        zhseqr_(&data->job,&data->compz,&data->n,&data->ilo,&data->ihi,data->H, &data->ldh,data->W,
                data->Z,&data->ldz,data->Work,&data->lwork,&data->info); 
        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument zhseqr is wrong\n", data->info );
        }
    }
  }
  
}

BENCHMARK_REGISTER_F(hseqr_test, zhseqr)->Apply(CustomArguments);