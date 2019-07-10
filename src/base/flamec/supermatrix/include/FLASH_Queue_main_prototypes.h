/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#ifndef FLASH_QUEUE_MAIN_PROTOTYPES_H
#define FLASH_QUEUE_MAIN_PROTOTYPES_H


void           FLASH_Queue_begin( void );
void           FLASH_Queue_end( void );
unsigned int   FLASH_Queue_stack_depth( void );

FLA_Error      FLASH_Queue_enable( void );
FLA_Error      FLASH_Queue_disable( void );
FLA_Bool       FLASH_Queue_get_enabled( void );

void           FLASH_Queue_set_num_threads( unsigned int n_threads );
unsigned int   FLASH_Queue_get_num_threads( void );


#ifdef FLA_ENABLE_SUPERMATRIX


void           FLASH_Queue_init( void );
void           FLASH_Queue_finalize( void );

unsigned int   FLASH_Queue_get_num_tasks( void );

void           FLASH_Queue_set_verbose_output( FLASH_Verbose verbose );
FLASH_Verbose  FLASH_Queue_get_verbose_output( void );
void           FLASH_Queue_set_sorting( FLA_Bool sorting );
FLA_Bool       FLASH_Queue_get_sorting( void );
void           FLASH_Queue_set_caching( FLA_Bool caching );
FLA_Bool       FLASH_Queue_get_caching( void );
void           FLASH_Queue_set_work_stealing( FLA_Bool work_stealing );
FLA_Bool       FLASH_Queue_get_work_stealing( void );
void           FLASH_Queue_set_data_affinity( FLASH_Data_aff data_affinity );
FLASH_Data_aff FLASH_Queue_get_data_affinity( void );
double         FLASH_Queue_get_total_time( void );
double         FLASH_Queue_get_parallel_time( void );

void           FLASH_Queue_exec( void );


// --- helper functions -------------------------------------------------------

void           FLASH_Queue_set_parallel_time( double dtime );
void           FLASH_Queue_set_block_size( dim_t size );
dim_t          FLASH_Queue_get_block_size( void );
void           FLASH_Queue_set_cache_size( dim_t size );
dim_t          FLASH_Queue_get_cache_size( void );
void           FLASH_Queue_set_cache_line_size( dim_t size );
dim_t          FLASH_Queue_get_cache_line_size( void );
void           FLASH_Queue_set_cores_per_cache( int cores );
int            FLASH_Queue_get_cores_per_cache( void );
void           FLASH_Queue_set_cores_per_queue( int cores );
int            FLASH_Queue_get_cores_per_queue( void );
void           FLASH_Queue_reset( void );
FLASH_Task*    FLASH_Queue_get_head_task( void );
FLASH_Task*    FLASH_Queue_get_tail_task( void );
void           FLASH_Queue_push( void *func, void *cntl, char *name,
                                 FLA_Bool enabled_gpu,
                                 int n_int_args, int n_fla_args,
                                 int n_input_args, int n_output_args, ... );
void           FLASH_Queue_push_input( FLA_Obj obj, FLASH_Task* t );
void           FLASH_Queue_push_output( FLA_Obj obj, FLASH_Task* t );
FLASH_Task*    FLASH_Task_alloc( void *func, void *cntl, char *name,
                                 FLA_Bool enabled_gpu,
                                 int n_int_args, int n_fla_args,
                                 int n_input_args, int n_output_args );
void           FLASH_Task_free( FLASH_Task *t );
void           FLASH_Queue_exec_task( FLASH_Task *t );
void           FLASH_Queue_verbose_output( void );

void           FLASH_Queue_init_tasks( void *arg );
void           FLASH_Queue_wait_enqueue( FLASH_Task *t, void *arg );
FLASH_Task*    FLASH_Queue_wait_dequeue( int queue, int cache, void *arg );
FLASH_Task*    FLASH_Queue_wait_dequeue_block( int queue, int cache, void *arg );
void           FLASH_Queue_update_cache( FLASH_Task *t, void *arg );
void           FLASH_Queue_update_cache_block( FLA_Obj obj, int cache, FLA_Bool output, void *arg );
void           FLASH_Queue_prefetch( int cache, void *arg );
void           FLASH_Queue_prefetch_block( FLA_Obj obj );
FLASH_Task*    FLASH_Queue_work_stealing( int queue, void *arg );
#ifdef FLA_ENABLE_GPU
void           FLASH_Queue_create_gpu( int thread, void *arg );
void           FLASH_Queue_destroy_gpu( int thread, void *arg );
FLA_Bool       FLASH_Queue_exec_gpu( FLASH_Task *t, void *arg );
FLA_Bool       FLASH_Queue_check_gpu( FLASH_Task *t, void *arg );
FLA_Bool       FLASH_Queue_check_block_gpu( FLA_Obj obj, int thread, void *arg );
void           FLASH_Queue_update_gpu( FLASH_Task *t, void **input_arg, void **output_arg, void *arg );
void           FLASH_Queue_update_block_gpu( FLA_Obj obj, void **buffer_gpu, int thread, void *arg );
void           FLASH_Queue_mark_gpu( FLASH_Task *t, void *arg );
void           FLASH_Queue_invalidate_block_gpu( FLA_Obj obj, int thread, void *arg );
void           FLASH_Queue_flush_block_gpu( FLA_Obj obj, int thread, void *arg );
void           FLASH_Queue_flush_gpu( int thread, void *arg );
#endif
void           FLASH_Queue_exec_parallel( void *arg );
void*          FLASH_Queue_exec_parallel_function( void *arg );
FLASH_Task*    FLASH_Task_update_dependencies( FLASH_Task *t, void *arg );
FLASH_Task*    FLASH_Task_update_binding( FLASH_Task *t, FLASH_Task *r, void *arg );
void           FLASH_Task_free_parallel( FLASH_Task *t, void *arg );

void           FLASH_Queue_exec_simulation( void *arg );


#endif // FLA_ENABLE_SUPERMATRIX


#endif // FLASH_QUEUE_MAIN_PROTOTYPES_H
