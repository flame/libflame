// Gbench headers
#include "benchmark/benchmark.h"
#include "gtest/gtest.h"

// Test suite headers
#include "../Flame_gbench_main.h"
#include "../Flame_gbench_aux.h"

using namespace std;

/* Macros */
#define gesvd_free() \
    if (a!=NULL)        free(a); \
    if (s!=NULL)        free(s); \
    if (u!=NULL)        free(u); \
    if (vt!=NULL)       free(vt); \
    if (rwork!=NULL)     free(rwork); \
    if (work!=NULL)     free(work); \
    if (superb!=NULL)   free(superb); \

#define gesvd_free_2() \
    if (data->a!=NULL)        free(data->a); \
    if (data->s!=NULL)        free(data->s); \
    if (data->u!=NULL)        free(data->u); \
    if (data->vt!=NULL)       free(data->vt); \
    if (data->rwork!=NULL)    free(data->rwork); \
    if (data->work!=NULL)     free(data->work); \
    if (data->superb!=NULL)   free(data->superb); \


/* Begin gesvd_parameters  class definition */
template <class T, class T2>
class gesvd_parameters{
   public:

    /* INPUT PARAMETERS */
    char jobu; // Must be 'A', 'S', 'O', or 'N'.
    char jobvt; // Must be 'A', 'S', 'O', or 'N'.
    int m; // The number of rows of the matrix A
    int n; // The number of columns of the matrix A
    int lda; //  The leading dimension of a
    int ldu; // leading dimension of the output array u . ldu≥ max(1, m)
    int ldvt; // The leading dimension of the output array vt . ldvt≥ max(1, p)
    T2 wkopt;
    int lwork;

    /* Helper variables */
    int forced_loop_count;

    /* Input / Output Buffers */
    T* a; // contains m-by-n matrix A.

    /* Output parameters */
    T2 *superb; // the dimension of subblocks.
    // Below buffers contain the o/p orthogonal/unitary matrices
    T *u;
    T *vt;
    T2 *s;
    T* work;
    T2 *rwork;

    /*Return Values */
    int info;

    public:
    //gesvd_parameters<T>::gesvd_parameters ( int nrow, int ncol, int lda_ );
    ~gesvd_parameters ()
    {
    #if FLAME_TEST_VERBOSE
       printf(" gesvd_parameters object: destructor invoked. \n");
    #endif
       //gesvd_free ();
    if (a!=NULL)        free(a);
    if (s!=NULL)        free(s);
    if (u!=NULL)        free(u);
    if (vt!=NULL)       free(vt);
    if (superb!=NULL)   free(superb);
       
    }
    gesvd_parameters ( char jobu_i, char jobvt_i, int m_i,
                    int n_i )
    {
       jobu = jobu_i;
       jobvt = jobvt_i;
       n = n_i;
       m = m_i;
       ldu = m;
       ldvt = n;
       lda = m;//n;//max(m,n);
       work = NULL;
       forced_loop_count = 1;
       if ( (m > FLAME_BIG_MATRIX_SIZE) || (n > FLAME_BIG_MATRIX_SIZE )){
          forced_loop_count = FLAME_GBENCH_FORCED_ITERATION_COUNT;
       }

    /* Memory allocation of the buffers */
       a = (T *)malloc(m*n*sizeof(T));
       s = (T2 *)malloc(n*sizeof(T2));
       rwork = (T2 *)malloc(5*(n*sizeof(T2))); // as per the API doc
       u = (T *)malloc(ldu*m*sizeof(T));
       vt = (T *)malloc(ldvt*n*sizeof(T));
       superb = (T2 *)malloc(fla_min(m,n)*sizeof(T2));
       lwork = -1;
       wkopt = 0.0;

       if( (a==NULL) || (s==NULL) ||  \
           (u==NULL) || (vt==NULL) || \
           (superb==NULL)  ){
          EXPECT_FALSE( true) << "gesvd_float_parameters object: malloc error. Exiting ";
          gesvd_free();
          exit(-1);
       }

       /* Initialization of input Buffers */
       Flame_gbench_init_buffer_rand( a, m*n );
/*     
       Flame_gbench_buffer_memset<T>( u, m*m, 0.0);
       Flame_gbench_buffer_memset<T>( vt, n*n, 0.0);
       Flame_gbench_buffer_memset<T>( s, n, 0.0);
       Flame_gbench_buffer_memset <T>( superb, fla_min(m,n), 0.0); */

    } /* end of Constructor  */
};  /* end of getrf_parameters  class definition */


/* Fixture definition  */
template <class T, class T2>
class gesvd_test : public ::benchmark::Fixture {
 public:
   int m,n;
   char u, vt;
   std::unique_ptr<gesvd_parameters<T, T2> >data;

 void SetUp(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
    m = static_cast<size_t>(state.range(0));
    n = static_cast<size_t>(state.range(1));
    u = (char)static_cast<size_t>(state.range(2));
    vt = (char)static_cast<size_t>(state.range(3));

    assert(data.get() == nullptr);
    data.reset(new gesvd_parameters<T, T2>(u, vt, m, n));

  }

  void TearDown(const ::benchmark::State& state) BENCHMARK_OVERRIDE {
      assert(data.get() != nullptr);
      data.reset();
  }

  ~gesvd_test() { assert(data == nullptr); }

};

// Customized argument generator for gesvd API
static void CustomArguments(benchmark::internal::Benchmark* b) {
    if (svd_paramslist[0].mode == MODE_RANGE) {
       for (int i = 0; i < svd_paramslist[0].num_tests ; i++){
           b->ArgsProduct({
              {benchmark::CreateDenseRange( svd_paramslist[i].m_range_start,
                                            svd_paramslist[i].m_range_end,
                                            svd_paramslist[i].m_range_step_size)} ,

              {benchmark::CreateDenseRange( svd_paramslist[i].n_range_start,
                                            svd_paramslist[i].n_range_end,
                                            svd_paramslist[i].n_range_step_size)},
              {svd_paramslist[i].jobu_gesvd},
                  {svd_paramslist[i].jobvt_gesvd}
              });
       }
    }
    if (svd_paramslist[0].mode == MODE_DISCRETE) {
       for (int i = 0; i < svd_paramslist[0].num_tests ; i++)
       {
           b->Args({svd_paramslist[i].m, svd_paramslist[i].n,
                    svd_paramslist[i].jobu_gesvd,
                    svd_paramslist[i].jobvt_gesvd});
       }
    }

    if (svd_paramslist[0].mode == MODE_COMBINATIONAL) {
       for (int i = 0; i < svd_paramslist[0].num_tests ; i++)
       {
          for (int j = 0; j < svd_paramslist[0].num_tests ; j++){
            b->Args({svd_paramslist[i].m, svd_paramslist[j].n,
            svd_paramslist[i].jobu_gesvd, svd_paramslist[i].jobvt_gesvd});
          }
       }
    }
}

/* Template Fixture for dgesvd API  */
BENCHMARK_TEMPLATE_DEFINE_F(gesvd_test, dgesvd, double, double)(benchmark::State &state) {
   assert(data.get() != nullptr);
   // This call is to compute the work buffer size.
   dgesvd_( (const char *)&data->jobu, (const char *)&data->jobvt, &data->m, &data->n, data->a, &data->lda, data->s,
                data->u, &data->ldu, data->vt, &data->ldvt, &data->wkopt, &data->lwork,
                &data->info );
   /* Allocate the memory for the work buffer*/
      data->lwork = (int)data->wkopt;
      data->work = (double*)malloc( (data->lwork)*sizeof(double) );
   
   if(data->work == NULL)
   {
      printf( "\n warning: dgesvd work buffer is not allocated\n", data->info );
	  gesvd_free_2();
	  exit(0);
   }

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
       /* Invoke the gesvd API for benchmarking */
       dgesvd_( (const char *)&data->jobu, (const char *)&data->jobvt, &data->m, &data->n, data->a, &data->lda, data->s,
                 data->u, &data->ldu, data->vt, &data->ldvt, data->work, &data->lwork,
                 &data->info );
        
        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument dgesvd is wrong\n", data->info );
        }
    }
  }
}

/* Template Fixture for sgesvd API  */
BENCHMARK_TEMPLATE_DEFINE_F(gesvd_test, sgesvd, float, float)(benchmark::State &state) {
  assert(data.get() != nullptr);
   // This call is to compute the work buffer size.
   sgesvd_( (const char *)&data->jobu, (const char *)&data->jobvt, &data->m, &data->n, data->a, &data->lda, data->s,
                data->u, &data->ldu, data->vt, &data->ldvt, &data->wkopt, &data->lwork,
                &data->info );
   /* Allocate the memory for the work buffer*/
      data->lwork = (int)(data->wkopt);
      data->work = (float*)malloc( (data->lwork)*sizeof(float) );
   
   if(data->work == NULL)
   {
      printf( "\n warning: sgesvd work buffer is not allocated\n", data->info );
	  gesvd_free_2();
	  exit(0);
   }

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
       /* Invoke the gesvd API for benchmarking */
       sgesvd_( (const char *)&data->jobu, (const char *)&data->jobvt, &data->m, &data->n, data->a, &data->lda, data->s,
                 data->u, &data->ldu, data->vt, &data->ldvt, data->work, &data->lwork,
                 &data->info );
        
        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument sgesvd is wrong\n", data->info );
        }
    }
  }
}

/* Template Fixture for cgetrf API  */
BENCHMARK_TEMPLATE_DEFINE_F(gesvd_test, cgesvd, lapack_complex_float, float)(benchmark::State &state) {
  assert(data.get() != nullptr);
   // This call is to compute the work buffer size.
   cgesvd_( (const char *)&data->jobu, (const char *)&data->jobvt, &data->m, &data->n, data->a, &data->lda, data->s,
                data->u, &data->ldu, data->vt, &data->ldvt, (lapack_complex_float*)&data->wkopt, &data->lwork, data->rwork,
                &data->info );
   /* Allocate the memory for the work buffer*/
      data->lwork = (int)data->wkopt;
      data->work = (lapack_complex_float*)malloc( (data->lwork)*sizeof(lapack_complex_float) );
   
   if(data->work == NULL)
   {
      printf( "\n warning: cgesvd work buffer is not allocated\n", data->info );
	  gesvd_free_2();
	  exit(0);
   }

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
       /* Invoke the gesvd API for benchmarking */
       cgesvd_( (const char *)&data->jobu, (const char *)&data->jobvt, &data->m, &data->n, data->a, &data->lda, data->s,
                 data->u, &data->ldu, data->vt, &data->ldvt, data->work, &data->lwork, data->rwork,
                 &data->info );
        
        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument cgesvd is wrong\n", data->info );
        }
    }
  }
}


/* Template Fixture for zgesvd API  */
BENCHMARK_TEMPLATE_DEFINE_F(gesvd_test, zgesvd, lapack_complex_double, double)(benchmark::State &state) {
  assert(data.get() != nullptr);
   // This call is to compute the work buffer size.
   zgesvd_( (const char *)&data->jobu, (const char *)&data->jobvt, &data->m, &data->n, data->a, &data->lda, data->s,
                data->u, &data->ldu, data->vt, &data->ldvt, (lapack_complex_double*)&data->wkopt, &data->lwork, data->rwork,
                &data->info );
   /* Allocate the memory for the work buffer*/
      data->lwork = (int)data->wkopt;
      data->work = (lapack_complex_double*)malloc( (data->lwork)*sizeof(lapack_complex_double) );
   
   if(data->work == NULL)
   {
      printf( "\n warning: zgesvd work buffer is not allocated\n", data->info );
	  gesvd_free_2();
	  exit(0);
   }

  for (auto _ : state) {
    for (int i=0; i<data->forced_loop_count; i++) {
       /* Invoke the gesvd API for benchmarking */
       zgesvd_( (const char *)&data->jobu, (const char *)&data->jobvt, &data->m, &data->n, data->a, &data->lda, data->s,
                 data->u, &data->ldu, data->vt, &data->ldvt, data->work, &data->lwork, data->rwork,
                 &data->info );
        
        if( data->info < 0 ) {
            printf( "\n warning: The %d th argument zgesvd is wrong\n", data->info );
        }
    }
  }
}

BENCHMARK_REGISTER_F(gesvd_test, dgesvd)->Apply(CustomArguments);
BENCHMARK_REGISTER_F(gesvd_test, zgesvd)->Apply(CustomArguments);
BENCHMARK_REGISTER_F(gesvd_test, cgesvd)->Apply(CustomArguments);
BENCHMARK_REGISTER_F(gesvd_test, sgesvd)->Apply(CustomArguments);
