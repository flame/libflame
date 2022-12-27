// Gbench headers
#include "benchmark/benchmark.h"
#include "gtest/gtest.h"

// Test suite headers
#include "../Flame_gbench_main.h"
#include "../Flame_gbench_aux.h"

using namespace std;

/* Macros */
#define geevx_parameters_free() \
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
   if (Scale!=NULL){ \
      free(Scale);  \
   } \
   if (Rconde!=NULL){ \
      free(Rconde);  \
   } \
   if (Rcondv!=NULL){ \
      free(Rcondv);  \
   } \
   if (Rwork!=NULL){ \
      free(Rwork);  \
   } \
   if (iwork!=NULL){ \
      free(iwork);  \
   } \

#define geevx_parameters_free_2() \
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
   if (data->Scale!=NULL){ \
      free(data->Scale);  \
   } \
   if (data->Rconde!=NULL){ \
      free(data->Rconde);  \
   } \
   if (data->Rcondv!=NULL){ \
      free(data->Rcondv);  \
   } \
   if (data->Rwork!=NULL){ \
      free(data->Rwork);  \
   } \
   if (data->iwork!=NULL){ \
      free(data->iwork);  \
   } \

/* Begin geevx_parameters  class definition */
template <class T,class S>

class geevx_parameters{
   public:
      char balanc,jobvl,jobvr,sense; 
      int  n,lda,ldvl,ldvr,lwork,info,ilo,ihi;
      int forced_loop_count;
      S Abnrm;

      /* Local arrays */
      T *A;
      T *W,*Vl,*Vr,*Work;
      S *Scale,*Rwork,*Rconde,*Rcondv;
      int *iwork;
   public:
    ~geevx_parameters ();
    geevx_parameters ( char balanc_i,char jobvl_i,char jobvr_i,char sense_i, int n_i)
    {
      balanc=balanc_i;
      jobvl=jobvl_i;
      jobvr=jobvr_i;
      sense=sense_i;
      if(sense=='E' || sense == 'B')
      {
          jobvl='V';
          jobvr='V';
      }
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
         Work=NULL;
         iwork = NULL;
         Scale = (S *)malloc((n)*sizeof(S)) ;
         Rconde = (S *)malloc((n)*sizeof(S)) ;
         Rcondv = (S *)malloc((n)*sizeof(S)) ;
         Rwork = (S *)malloc((2*n)*sizeof(S)) ;
      if ((A==NULL) || (Vl==NULL)|| (Vr==NULL)|| (W==NULL)|| (Scale==NULL)|| (Rconde==NULL) || 
          (Rcondv==NULL) || (Rwork==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        geevx_parameters_free();
        exit(0);
      }
    /* Initialization of input Buffers */
      Flame_gbench_init_buffer_rand( A, n*lda);
    } /* end of Constructor  */
};  /* end of geevx_parameters  class definition */



/* Destructor definition  'geevx_parameters' class  */
template <class T,class S>
geevx_parameters<T,S>:: ~geevx_parameters ()
{
#if FLAME_TEST_VERBOSE
   printf(" geevx_parameters object: destructor invoked. \n");
#endif
   geevx_parameters_free ();
}

/* Fixture definition  */
template <class T,class S>
class geevx_test : public ::benchmark::Fixture {
 public:
   int n;
   char balanc,jobvl,jobvr,sense;
   std::unique_ptr<geevx_parameters<T,S>> data;

 void SetUp(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      n   = static_cast<size_t>(state.range(0));
      balanc = (char)static_cast<size_t>(state.range(1));
      jobvl = (char)static_cast<size_t>(state.range(2));  
      jobvr = (char)static_cast<size_t>(state.range(3));
      sense = (char)static_cast<size_t>(state.range(4));
          
    assert(data.get() == nullptr);
    data.reset(new geevx_parameters<T,S>(balanc,jobvl,jobvr,sense,n));

  }

  void TearDown(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      assert(data.get() != nullptr);
      data.reset();
  }

  ~geevx_test() { assert(data == nullptr); }

};

// Customized argument generator for geevx API
static void CustomArguments(benchmark::internal::Benchmark* b) {
    
    if (eig_non_sym_paramslist[0].mode == MODE_RANGE) {
       for (int i = 0; i < eig_non_sym_paramslist[0].num_tests ; i++){
           b->ArgsProduct({
              {benchmark::CreateDenseRange( eig_non_sym_paramslist[i].n_range_start,
                                            eig_non_sym_paramslist[i].n_range_end,
                                            eig_non_sym_paramslist[i].n_range_step_size)},
                                            {'N','P','S','B'},
                                            {'N', 'V'},
                                            {'N', 'V'},
                                            {'N','E','V','B'}
              });
       }
    }
    
    if (eig_non_sym_paramslist[0].mode == MODE_DISCRETE) {
       for (int i = 0; i < eig_non_sym_paramslist[0].num_tests ; i++)
       {
           b->Args({eig_non_sym_paramslist[i].n,eig_non_sym_paramslist[i].balance_ggevx,eig_non_sym_paramslist[i].jobvsl, 
                    eig_non_sym_paramslist[i].jobvsr,eig_non_sym_paramslist[i].sense_ggevx});
       }
    }
    if (eig_non_sym_paramslist[0].mode == MODE_COMBINATIONAL) {
       for (int i = 0; i < eig_non_sym_paramslist[0].num_tests ; i++)
       {
          for (int j = 0; j < eig_non_sym_paramslist[0].num_tests ; j++){
              for(int k=0; k< eig_non_sym_paramslist[0].num_tests ; k++){
                  for(int l=0; l< eig_non_sym_paramslist[0].num_tests ; l++){
                      for(int m=0; m< eig_non_sym_paramslist[0].num_tests ; m++){
                        b->Args({eig_non_sym_paramslist[i].n,eig_non_sym_paramslist[j].balance_ggevx, 
                                eig_non_sym_paramslist[k].jobvsl,eig_non_sym_paramslist[l].jobvsr,
                                eig_non_sym_paramslist[m].sense_ggevx});
                      }
                  }
              }  
                    
            }
       }
    }
}

/* Template Fixture for sgeevx API  */
BENCHMARK_TEMPLATE_DEFINE_F(geevx_test, sgeevx, float,float )(benchmark::State &state) {
  assert(data.get() != nullptr);
      float work_query=-1;

      data->iwork=(int *)malloc(sizeof(int)*(2*(data->n)));
      if (data->iwork == NULL){
        printf("error of memory allocation. in Work Exiting ...\n");
        geevx_parameters_free_2();
        exit(0);
    }
    /* Query for optimal working array(s) size */
   sgeevx_( &data->balanc,&data->jobvl, &data->jobvr,&data->sense,&data->n,data->A, &data->lda,
           data->W,data->Rwork,data->Vl, &data->ldvl, data->Vr,&data->ldvr, &data->ilo,&data->ihi,data->Scale,
           &data->Abnrm,data->Rconde,data->Rcondv,(float *)&work_query, &data->lwork,data->iwork,&data->info); 

    data->lwork = (int)work_query;
    data->lwork = fla_max(((2*data->n)+((data->n)^2)), data->lwork);
    /* Allocate memory for work arrays */
    data->Work = (float *)malloc( sizeof(float) * (data->lwork));
    if (data->Work == NULL){
        printf("error of memory allocation. in Work Exiting ...\n");
        geevx_parameters_free_2();
        exit(0);
    }

  for (auto _ : state) {

    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the geevx API for benchmarking */
        sgeevx_( &data->balanc,&data->jobvl, &data->jobvr,&data->sense,&data->n,data->A, &data->lda,
           data->W,data->Rwork,data->Vl, &data->ldvl, data->Vr,&data->ldvr, &data->ilo,&data->ihi,data->Scale,
           &data->Abnrm,data->Rconde,data->Rcondv,data->Work, &data->lwork,data->iwork,&data->info); 
        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument sgeevx is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(geevx_test, sgeevx)->Apply(CustomArguments);

/* Template Fixture for dgeevx API  */
BENCHMARK_TEMPLATE_DEFINE_F(geevx_test, dgeevx, double,double )(benchmark::State &state) {
   assert(data.get() != nullptr);

   double work_query=-1;
   data->iwork=(int *)malloc(sizeof(int)*(2*(data->n)));
   if (data->iwork == NULL){
        printf("error of memory allocation. in Work Exiting ...\n");
        geevx_parameters_free_2();
        exit(0);
    }
   /* Query for optimal working array(s) size */
   dgeevx_( &data->balanc,&data->jobvl, &data->jobvr,&data->sense,&data->n,data->A, &data->lda,
           data->W,data->Rwork,data->Vl, &data->ldvl, data->Vr,&data->ldvr, &data->ilo,&data->ihi,data->Scale,
           &data->Abnrm,data->Rconde,data->Rcondv,(double *)&work_query, &data->lwork,data->iwork,&data->info); 
 

    data->lwork = (int)work_query;
    data->lwork = fla_max(((2*data->n)+((data->n)^2)), data->lwork);

    /* Allocate memory for work arrays */
    data->Work = (double *)malloc( sizeof(double) * (data->lwork));
    if (data->Work == NULL){
        printf("error of memory allocation. in Work Exiting ...\n");
        geevx_parameters_free_2();
        exit(0);
    }

  for (auto _ : state) {

    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the geevx API for benchmarking */
        dgeevx_( &data->balanc,&data->jobvl, &data->jobvr,&data->sense,&data->n,data->A, &data->lda,
           data->W,data->Rwork,data->Vl, &data->ldvl, data->Vr,&data->ldvr, &data->ilo,&data->ihi,data->Scale,
           &data->Abnrm,data->Rconde,data->Rcondv,data->Work, &data->lwork,data->iwork,&data->info); 

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument dgeevx is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(geevx_test, dgeevx)->Apply(CustomArguments);

/* Template Fixture for cgeevx API  */
BENCHMARK_TEMPLATE_DEFINE_F(geevx_test, cgeevx,lapack_complex_float, float )(benchmark::State &state) {
  assert(data.get() != nullptr);
      float work_query=-1;
    /* Query for optimal working array(s) size */
   cgeevx_( &data->balanc,&data->jobvl, &data->jobvr,&data->sense,&data->n,data->A, &data->lda,
           data->W,data->Vl, &data->ldvl, data->Vr,&data->ldvr, &data->ilo,&data->ihi,data->Scale,
           &data->Abnrm,data->Rconde,data->Rcondv,(lapack_complex_float *)&work_query, &data->lwork,data->Rwork,&data->info); 
 

    data->lwork = (int)work_query;
    data->lwork = fla_max(((2*data->n)+((data->n)^2)), data->lwork);
    /* Allocate memory for work arrays */
    data->Work = (lapack_complex_float *)malloc( sizeof(lapack_complex_float) * (data->lwork));
    if (data->Work == NULL){
        printf("error of memory allocation. in Work Exiting ...\n");
        geevx_parameters_free_2();
        exit(0);
    }
    

  for (auto _ : state) {

    for (int i=0; i<data->forced_loop_count; i++) {
        /* Invoke the geevx API for benchmarking */
        cgeevx_( &data->balanc,&data->jobvl, &data->jobvr,&data->sense,&data->n,data->A, &data->lda,
           data->W,data->Vl, &data->ldvl, data->Vr,&data->ldvr, &data->ilo,&data->ihi,data->Scale,
           &data->Abnrm,data->Rconde,data->Rcondv,data->Work, &data->lwork,data->Rwork,&data->info);

        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument cgeevx is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(geevx_test, cgeevx)->Apply(CustomArguments);

/* Template Fixture for zgeevx API  */
BENCHMARK_TEMPLATE_DEFINE_F(geevx_test, zgeevx,lapack_complex_double, double)(benchmark::State &state) {
    assert(data.get() != nullptr);

   double work_query=-1;
    /* Query for optimal working array(s) size */
    zgeevx_( &data->balanc,&data->jobvl, &data->jobvr,&data->sense,&data->n,data->A, &data->lda,
           data->W,data->Vl, &data->ldvl, data->Vr,&data->ldvr, &data->ilo,&data->ihi,data->Scale,
           &data->Abnrm,data->Rconde,data->Rcondv,(lapack_complex_double *)&work_query, &data->lwork,data->Rwork,&data->info); 
 
                                    
    data->lwork = (int)work_query;
    data->lwork = fla_max(((2*data->n)+((data->n)^2)), data->lwork);
    /* Allocate memory for work arrays */
    data->Work = (lapack_complex_double *)malloc( sizeof(lapack_complex_double) * (data->lwork));
    if ((data->Work==NULL)){
        printf("error of memory allocation. Exiting ...\n");
        geevx_parameters_free_2();
        exit(0);
    }
    for (auto _ : state) {
        
        for (int i=0; i<data->forced_loop_count; i++) {
           /* Invoke the geevx API for benchmarking */
            zgeevx_( &data->balanc,&data->jobvl, &data->jobvr,&data->sense,&data->n,data->A, &data->lda,
           data->W,data->Vl, &data->ldvl, data->Vr,&data->ldvr, &data->ilo,&data->ihi,data->Scale,
           &data->Abnrm,data->Rconde,data->Rcondv,data->Work, &data->lwork,data->Rwork,&data->info);        
            if( data->info < 0 ) {
                printf( "\n warning: The %d th argument zgeevx is wrong\n", data->info );
            }
        }
    }
}

BENCHMARK_REGISTER_F(geevx_test, zgeevx)->Apply(CustomArguments);
