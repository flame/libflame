/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"


#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP
#ifdef FLA_ENABLE_TIDSP
#include <ti/omp/omp.h>
#else
#include <omp.h>
#endif
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS
#include <pthread.h>
#endif


#ifdef FLA_ENABLE_SUPERMATRIX

FLASH_Queue           _tq;

static FLA_Bool       flash_queue_initialized     = FALSE;

static int            flash_queue_n_read_blocks   = 0;
static int            flash_queue_n_write_blocks  = 0;

static FLASH_Verbose  flash_queue_verbose         = FLASH_QUEUE_VERBOSE_NONE;
static FLA_Bool       flash_queue_sorting         = FALSE;
static FLA_Bool       flash_queue_caching         = FALSE;
static FLA_Bool       flash_queue_work_stealing   = FALSE;
static FLASH_Data_aff flash_queue_data_affinity   = FLASH_QUEUE_AFFINITY_NONE;

static double         flash_queue_total_time      = 0.0;
static double         flash_queue_parallel_time   = 0.0;

static dim_t          flash_queue_block_size      = 0;
static dim_t          flash_queue_cache_size      = 2 * 1024 * 1024;
static dim_t          flash_queue_cache_line_size = 64;

static int            flash_queue_cores_per_cache = 1;
static int            flash_queue_cores_per_queue = 0;

#endif


static unsigned int   flash_queue_stack           = 0;
static FLA_Bool       flash_queue_enabled         = TRUE;

static unsigned int   flash_queue_n_threads       = 1;


void FLASH_Queue_begin( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_begin

----------------------------------------------------------------------------*/
{
#ifdef FLA_ENABLE_SUPERMATRIX
   if ( flash_queue_stack == 0 )
   {
      // Save the starting time for the total execution time.
      flash_queue_total_time = FLA_Clock();
   }
#endif

   // Push onto the stack.
   flash_queue_stack++;

   return;
}


void FLASH_Queue_end( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_end

----------------------------------------------------------------------------*/
{
   // Pop off the stack.
   flash_queue_stack--;

#ifdef FLA_ENABLE_SUPERMATRIX
   if ( flash_queue_stack == 0 )
   {
      // Execute tasks if encounter the outermost parallel region.
      FLASH_Queue_exec();

      // Find the total execution time.
      flash_queue_total_time = FLA_Clock() - flash_queue_total_time;
   }
#endif

   return;
}


unsigned int FLASH_Queue_stack_depth( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_stack_depth

----------------------------------------------------------------------------*/
{
   return flash_queue_stack;
}


FLA_Error FLASH_Queue_enable( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_enable

----------------------------------------------------------------------------*/
{
#ifdef FLA_ENABLE_SUPERMATRIX
   if ( flash_queue_stack == 0 )
   {
      // Enable if not begin parallel region yet.
      flash_queue_enabled = TRUE;   
      return FLA_SUCCESS;
   }
   else
   {
      // Cannot change status during parallel region.
      return FLA_FAILURE;
   }
#else
   // Raise an exception when SuperMatrix is not configured.
   FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED );
   return FLA_FAILURE;
#endif
}


FLA_Error FLASH_Queue_disable( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_disable

----------------------------------------------------------------------------*/
{
#ifdef FLA_ENABLE_SUPERMATRIX
   if ( flash_queue_stack == 0 )
   {
      // Disable if not begin parallel region yet.
      flash_queue_enabled = FALSE;   
      return FLA_SUCCESS;      
   }
   else
   {
      // Cannot change status during parallel region.
      return FLA_FAILURE;
   }
#else
   // Allow disabling enqueuing even when SuperMatrix is not configured.
   flash_queue_enabled = FALSE;   
   return FLA_SUCCESS;
#endif
}


FLA_Bool FLASH_Queue_get_enabled( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_enabled

----------------------------------------------------------------------------*/
{
   // Return if enabled, but always false if SuperMatrix is not configured.
#ifdef FLA_ENABLE_SUPERMATRIX
   return flash_queue_enabled;
#else
   return FALSE;
#endif
}


void FLASH_Queue_set_num_threads( unsigned int n_threads )
/*----------------------------------------------------------------------------

   FLASH_Queue_set_num_threads

----------------------------------------------------------------------------*/
{
   FLA_Error e_val;

   // Verify that the number of threads is positive. 
   e_val = FLA_Check_num_threads( n_threads );
   FLA_Check_error_code( e_val );

   // Keep track of the number of threads internally.
   flash_queue_n_threads = n_threads;

#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP

   // No additional action is necessary to set the number of OpenMP threads
   // since setting the number of threads is handled at the parallel for loop
   // with a num_threads() clause. This gives the user more flexibility since
   // he can use the OMP_NUM_THREADS environment variable or the
   // omp_set_num_threads() function to set the global number of OpenMP threads
   // independently of the number of SuperMatrix threads.
   
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS

   // No additional action is necessary to set the number of pthreads
   // since setting the number of threads is handled entirely on our end.

#endif

   return;
}


unsigned int FLASH_Queue_get_num_threads( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_num_threads

----------------------------------------------------------------------------*/
{
   return flash_queue_n_threads;
}


#ifdef FLA_ENABLE_SUPERMATRIX


void FLASH_Queue_init( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_init

----------------------------------------------------------------------------*/
{
   // Exit early if we're already initialized.
   if ( flash_queue_initialized == TRUE )
      return;
   
   // Reset all the initial values.
   FLASH_Queue_reset();

   // Set the initialized flag.
   flash_queue_initialized = TRUE;

#ifdef FLA_ENABLE_GPU
   FLASH_Queue_init_gpu();
#endif

#ifdef FLA_ENABLE_HIP
   FLASH_Queue_init_hip();
#endif

   return;
}


void FLASH_Queue_finalize( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_finalize

----------------------------------------------------------------------------*/
{
   // Exit early if we're not already initialized.
   if ( flash_queue_initialized == FALSE )
      return;

   // Clear the initialized flag.
   flash_queue_initialized = FALSE;

#ifdef FLA_ENABLE_GPU
   FLASH_Queue_finalize_gpu();
#endif

#ifdef FLA_ENABLE_HIP
   FLASH_Queue_finalize_hip();
#endif

   return;
}


unsigned int FLASH_Queue_get_num_tasks( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_num_tasks

----------------------------------------------------------------------------*/
{
   return _tq.n_tasks;
}


void FLASH_Queue_set_verbose_output( FLASH_Verbose verbose )
/*----------------------------------------------------------------------------

   FLASH_Queue_set_verbose_output

----------------------------------------------------------------------------*/
{ 
   flash_queue_verbose = verbose;

   return;
}


FLASH_Verbose FLASH_Queue_get_verbose_output( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_verbose_output

----------------------------------------------------------------------------*/
{ 
   return flash_queue_verbose;
}


void FLASH_Queue_set_sorting( FLA_Bool sorting )
/*----------------------------------------------------------------------------

   FLASH_Queue_set_sorting

----------------------------------------------------------------------------*/
{ 
   flash_queue_sorting = sorting; 

   return;
}


FLA_Bool FLASH_Queue_get_sorting( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_sorting

----------------------------------------------------------------------------*/
{ 
   return flash_queue_sorting;
}


void FLASH_Queue_set_caching( FLA_Bool caching )
/*----------------------------------------------------------------------------

   FLASH_Queue_set_caching

----------------------------------------------------------------------------*/
{ 
   flash_queue_caching = caching; 

   return;
}


FLA_Bool FLASH_Queue_get_caching( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_caching

----------------------------------------------------------------------------*/
{ 
   return flash_queue_caching;
}


void FLASH_Queue_set_work_stealing( FLA_Bool work_stealing )
/*----------------------------------------------------------------------------

   FLASH_Queue_set_work_stealing

----------------------------------------------------------------------------*/
{
   flash_queue_work_stealing = work_stealing;

   return;
}


FLA_Bool FLASH_Queue_get_work_stealing( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_work_stealing

----------------------------------------------------------------------------*/
{
   return flash_queue_work_stealing;
}


void FLASH_Queue_set_data_affinity( FLASH_Data_aff data_affinity )
/*----------------------------------------------------------------------------

   FLASH_Queue_set_data_affinity

----------------------------------------------------------------------------*/
{ 
   flash_queue_data_affinity = data_affinity; 

   return;
}


FLASH_Data_aff FLASH_Queue_get_data_affinity( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_data_affinity

----------------------------------------------------------------------------*/
{ 
   return flash_queue_data_affinity;
}


double FLASH_Queue_get_total_time( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_total_time

----------------------------------------------------------------------------*/
{
   // Only return time if out of parallel region.
   if ( flash_queue_stack == 0 )
      return flash_queue_total_time;

   return 0.0;
}


double FLASH_Queue_get_parallel_time( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_parallel_time

----------------------------------------------------------------------------*/
{
   // Only return time if out of parallel region.
   if ( flash_queue_stack == 0 )
      return flash_queue_parallel_time;

   return 0.0;
}


// --- helper functions --- ===================================================


void FLASH_Queue_set_parallel_time( double dtime )
/*----------------------------------------------------------------------------

   FLASH_Queue_set_parallel_time

----------------------------------------------------------------------------*/
{
   flash_queue_parallel_time = dtime;

   return;
}


void FLASH_Queue_set_block_size( dim_t size )
/*----------------------------------------------------------------------------

   FLASH_Queue_set_block_size

----------------------------------------------------------------------------*/
{
   // Only adjust the block size if the new block is larger.
   if ( flash_queue_block_size < size )
      flash_queue_block_size = size;

   return;
}


dim_t FLASH_Queue_get_block_size( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_block_size

----------------------------------------------------------------------------*/
{
   return flash_queue_block_size;
}


void FLASH_Queue_set_cache_size( dim_t size )
/*----------------------------------------------------------------------------

   FLASH_Queue_set_cache_size

----------------------------------------------------------------------------*/
{
   flash_queue_cache_size = size;

   return;
}


dim_t FLASH_Queue_get_cache_size( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_cache_size

----------------------------------------------------------------------------*/
{
   return flash_queue_cache_size;
}


void FLASH_Queue_set_cache_line_size( dim_t size )
/*----------------------------------------------------------------------------

   FLASH_Queue_set_cache_line_size

----------------------------------------------------------------------------*/
{
   flash_queue_cache_line_size = size;

   return;
}


dim_t FLASH_Queue_get_cache_line_size( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_cache_line_size

----------------------------------------------------------------------------*/
{
   return flash_queue_cache_line_size;
}


void FLASH_Queue_set_cores_per_cache( int cores )
/*----------------------------------------------------------------------------

   FLASH_Queue_set_cores_per_cache

----------------------------------------------------------------------------*/
{
   flash_queue_cores_per_cache = cores;

   return;
}


int FLASH_Queue_get_cores_per_cache( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_cores_per_cache

----------------------------------------------------------------------------*/
{
   return flash_queue_cores_per_cache;
}


void FLASH_Queue_set_cores_per_queue( int cores )
/*----------------------------------------------------------------------------

   FLASH_Queue_set_cores_per_queue

----------------------------------------------------------------------------*/
{
   flash_queue_cores_per_queue = cores;

   return;
}


int FLASH_Queue_get_cores_per_queue( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_cores_per_queue

----------------------------------------------------------------------------*/
{
   return flash_queue_cores_per_queue;
}


void FLASH_Queue_reset( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_reset

----------------------------------------------------------------------------*/
{
   // Clear the other fields of the FLASH_Queue structure.
   _tq.n_tasks = 0;
   _tq.head    = NULL;
   _tq.tail    = NULL;

   // Reset the number of blocks.
   flash_queue_n_read_blocks  = 0;
   flash_queue_n_write_blocks = 0;

   return;
}


FLASH_Task* FLASH_Queue_get_head_task( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_head_task

----------------------------------------------------------------------------*/
{
   return _tq.head;
}


FLASH_Task* FLASH_Queue_get_tail_task( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_tail_task

----------------------------------------------------------------------------*/
{
   return _tq.tail;
}


void FLASH_Queue_push( void* func,
                       void* cntl,
                       char* name,
                       FLA_Bool enabled_gpu,
                       FLA_Bool enabled_hip,
                       int n_int_args,
                       int n_fla_args,
                       int n_input_args,
                       int n_output_args,
                       ... )
/*----------------------------------------------------------------------------

   FLASH_Queue_push

----------------------------------------------------------------------------*/
{
   int         i;
   va_list     var_arg_list;
   FLASH_Task* t;
   FLA_Obj     obj;

   // Allocate a new FLA_Task and populate its fields with appropriate values.
   t = FLASH_Task_alloc( func, cntl, name, enabled_gpu,
                         enabled_hip, n_int_args, n_fla_args,
                         n_input_args, n_output_args );
   
   // Initialize variable argument environment. In case you're wondering, the
   // second argument in this macro invocation of va_start() is supposed to be
   // the parameter that immediately preceeds the variable argument list
   // (ie: the ... above ).
   va_start( var_arg_list, n_output_args );

   // Extract the integer arguments.
   for ( i = 0; i < n_int_args; i++ )
      t->int_arg[i] = va_arg( var_arg_list, int );
   
   // Extract the FLA_Obj arguments.
   for ( i = 0; i < n_fla_args; i++ )
      t->fla_arg[i] = va_arg( var_arg_list, FLA_Obj );

   // Extract the input FLA_Obj arguments.
   for ( i = 0; i < n_input_args; i++ )
   {
      obj = va_arg( var_arg_list, FLA_Obj );
      t->input_arg[i] = obj;

      // Macroblock is used.
      if ( FLA_Obj_elemtype( obj ) == FLA_MATRIX )
      {
         dim_t    jj, kk;
         dim_t    m    = FLA_Obj_length( obj );
         dim_t    n    = FLA_Obj_width( obj );
         dim_t    cs   = FLA_Obj_col_stride( obj );
         FLA_Obj* buf  = FLASH_OBJ_PTR_AT( obj );

         // Dependence analysis for each input block in macroblock.
         for ( jj = 0; jj < n; jj++ )
            for ( kk = 0; kk < m; kk++ )
               FLASH_Queue_push_input( *( buf + jj * cs + kk ), t );

         // Set the number of blocks in the macroblock subtracted by one
         // since we do not want to recount an operand for each n_input_arg.
         t->n_macro_args += m * n - 1;
      }      
      else // Regular block.
      {
         // Dependence analysis for input operand.
         FLASH_Queue_push_input( obj, t );
      }
   }

   // Extract the output FLA_Obj arguments.
   for ( i = 0; i < n_output_args; i++ )
   {
      obj = va_arg( var_arg_list, FLA_Obj );
      t->output_arg[i] = obj;

      // Only assign data affinity to the first output block.
      if ( i == 0 )
      {
         FLA_Obj buf = obj;

         // Use the top left block of the macroblock.
         if ( FLA_Obj_elemtype( obj ) == FLA_MATRIX )
            buf = *FLASH_OBJ_PTR_AT( obj );
         
         if ( buf.base->write_task == NULL )
            t->queue = flash_queue_n_write_blocks;
         else 
            t->queue = buf.base->write_task->queue;
      }

      // Macroblock is used.
      if ( FLA_Obj_elemtype( obj ) == FLA_MATRIX )
      {
         dim_t    jj, kk;
         dim_t    m    = FLA_Obj_length( obj );
         dim_t    n    = FLA_Obj_width( obj );
         dim_t    cs   = FLA_Obj_col_stride( obj );
         FLA_Obj* buf  = FLASH_OBJ_PTR_AT( obj );

         // Dependence analysis for each output block in macroblock.
         for ( jj = 0; jj < n; jj++ )
            for ( kk = 0; kk < m; kk++ )
               FLASH_Queue_push_output( *( buf + jj * cs + kk ), t );

         // Set the number of blocks in the macroblock subtracted by one
         // since we do not want to recount an operand for each n_output_arg.
         t->n_macro_args += m * n - 1;
      }      
      else // Regular block.
      {
         // Dependence analysis for output operand.
         FLASH_Queue_push_output( obj, t );
      }
   }      

   // Finalize the variable argument environment.
   va_end( var_arg_list );
  
   // Add the task to the tail of the queue (and the head if queue is empty).
   if ( _tq.n_tasks == 0 )
   {
      _tq.head = t;
      _tq.tail = t;
   }
   else
   {
      t->prev_task = _tq.tail;
      _tq.tail->next_task = t;
      _tq.tail            = t;

      // Determine the index of the task in the task queue.
      t->order = t->prev_task->order + 1;
   }
   
   // Increment the number of tasks.
   _tq.n_tasks++;

   return;
}


void FLASH_Queue_push_input( FLA_Obj obj,
                             FLASH_Task* t )
/*----------------------------------------------------------------------------

   FLASH_Queue_push_input

----------------------------------------------------------------------------*/
{
   FLASH_Task* task;
   FLASH_Dep*  d;

   // Find dependence information.
   if ( obj.base->write_task == NULL )
   {
      t->n_ready--;
      
      // Add to number of blocks read if not written and not read before.
      if ( obj.base->n_read_tasks == 0 )
      {
         // Identify each read block with an id for freeing.
         obj.base->n_read_blocks = flash_queue_n_read_blocks;
         
         flash_queue_n_read_blocks++;            
      }
   }
   else
   { // Flow dependence.
      task = obj.base->write_task;
      
      d = (FLASH_Dep *) FLA_malloc( sizeof(FLASH_Dep) );
      
      d->task     = t;
      d->next_dep = NULL;
      
      if ( task->n_dep_args == 0 )
      {
         task->dep_arg_head = d;
         task->dep_arg_tail = d;
      }
      else
      {
         task->dep_arg_tail->next_dep = d;
         task->dep_arg_tail           = d;
      }
      
      task->n_dep_args++;
   }
   
   // Add task to the read task in the object if not already there.
   if ( obj.base->n_read_tasks == 0 ||
        obj.base->read_task_tail->task != t )
   { // Anti-dependence potentially.
      d = (FLASH_Dep *) FLA_malloc( sizeof(FLASH_Dep) );
      
      d->task     = t;
      d->next_dep = NULL;
      
      if ( obj.base->n_read_tasks == 0 )
      {
         obj.base->read_task_head = d;
         obj.base->read_task_tail = d;
      }
      else
      {
         obj.base->read_task_tail->next_dep = d;
         obj.base->read_task_tail           = d;
      }
      
      obj.base->n_read_tasks++;
   }      
   
   return;
}


void FLASH_Queue_push_output( FLA_Obj obj,
                              FLASH_Task* t )
/*----------------------------------------------------------------------------

   FLASH_Queue_push_output

----------------------------------------------------------------------------*/
{
   int         i;
   FLASH_Task* task;
   FLASH_Dep*  d;
   FLASH_Dep*  next_dep;

   // Assign tasks to threads with data affinity.
   if ( obj.base->write_task == NULL )
   {
      t->n_ready--;
      
      // Save index in which this output block is first encountered.
      obj.base->n_write_blocks = flash_queue_n_write_blocks;
      
      // Number of blocks written if not written before.
      flash_queue_n_write_blocks++;
      
      // Add to number of blocks read if not written or read before.
      if ( obj.base->n_read_tasks == 0 )
      {
         // Identify each read block with an id for freeing.
         obj.base->n_read_blocks = flash_queue_n_read_blocks;
         
         flash_queue_n_read_blocks++;
      }
   }
   else
   { // Flow dependence potentially.
      // The last task to overwrite this block is not itself.
      if ( obj.base->write_task != t )
      {
         // Create dependency from task that last wrote the block.
         task = obj.base->write_task;
         
         d = (FLASH_Dep *) FLA_malloc( sizeof(FLASH_Dep) );
         
         d->task     = t;
         d->next_dep = NULL;
         
         if ( task->n_dep_args == 0 )
         {
            task->dep_arg_head = d;
            task->dep_arg_tail = d;
         }
         else
         {
            task->dep_arg_tail->next_dep = d;
            task->dep_arg_tail           = d;
         }
         
         task->n_dep_args++;
      }
      else
      {
         // No need to notify task twice for output block already seen.
         t->n_ready--;
      }
   }
   
   // Clear read task for next set of reads and record the anti-dependence.
   d = obj.base->read_task_head;
   
   for ( i = 0; i < obj.base->n_read_tasks; i++ )
   {
      task     = d->task;
      next_dep = d->next_dep;
      
      // If the last task to read is not the current task, add dependence.
      if ( task != t )
      {
         d->task     = t;
         d->next_dep = NULL;
         
         if ( task->n_dep_args == 0 )
         {
            task->dep_arg_head = d;
            task->dep_arg_tail = d;
         }
         else
         {
            task->dep_arg_tail->next_dep = d;
            task->dep_arg_tail           = d;
         }
         
         task->n_dep_args++;
         
         t->n_war_args++;
      }  
      else
      {
         FLA_free( d );
      }
      
      d = next_dep;
   }
   
   obj.base->n_read_tasks   = 0;
   obj.base->read_task_head = NULL;
   obj.base->read_task_tail = NULL;
   
   // Record this task as the last to write to this block.
   obj.base->write_task = t;
   
   return;
}


FLASH_Task* FLASH_Task_alloc( void *func,
                              void *cntl,
                              char *name,
                              FLA_Bool enabled_gpu,
                              FLA_Bool enabled_hip,
                              int n_int_args,
                              int n_fla_args,
                              int n_input_args,
                              int n_output_args )
/*----------------------------------------------------------------------------

   FLASH_Task_alloc

----------------------------------------------------------------------------*/
{
   FLASH_Task* t;

   // Allocate space for the task structure t.
   t             = (FLASH_Task *) FLA_malloc( sizeof(FLASH_Task) );

   // Allocate space for the task's integer arguments.
   t->int_arg    = (int *) FLA_malloc( n_int_args * sizeof(int) );

   // Allocate space for the task's FLA_Obj arguments.
   t->fla_arg    = (FLA_Obj *) FLA_malloc( n_fla_args * sizeof(FLA_Obj) );

   // Allocate space for the task's input FLA_Obj arguments.
   t->input_arg  = (FLA_Obj *) FLA_malloc( n_input_args * sizeof(FLA_Obj) );

   // Allocate space for the task's output FLA_Obj arguments.
   t->output_arg = (FLA_Obj *) FLA_malloc( n_output_args * sizeof(FLA_Obj) );
   
   // Initialize other fields of the structure.
   t->n_ready       = 0;
   t->order         = 0;
   t->queue         = 0;
   t->height        = 0;
   t->thread        = 0;
   t->cache         = 0;
   t->hit           = FALSE;

   t->func          = func;
   t->cntl          = cntl;
   t->name          = name;
   t->enabled_gpu   = enabled_gpu;
   t->enabled_hip   = enabled_hip;
   t->n_int_args    = n_int_args;
   t->n_fla_args    = n_fla_args;
   t->n_input_args  = n_input_args;
   t->n_output_args = n_output_args;

   t->n_macro_args  = 0;
   t->n_war_args    = 0;
   t->n_dep_args    = 0;
   t->dep_arg_head  = NULL;
   t->dep_arg_tail  = NULL;
   t->prev_task     = NULL;
   t->next_task     = NULL;
   t->prev_wait     = NULL;
   t->next_wait     = NULL;
   
   // Return a pointer to the initialized structure.
   return t;
}


void FLASH_Task_free( FLASH_Task *t )
/*----------------------------------------------------------------------------

   FLASH_Task_free

----------------------------------------------------------------------------*/
{
   int        i, j, k;
   FLA_Obj    obj;
   FLASH_Dep* d;
   FLASH_Dep* next_dep;

   // Clearing the last write task in each output block.
   for ( i = 0; i < t->n_output_args; i++ )
   {
      obj = t->output_arg[i];
      
      // Macroblock is used.
      if ( FLA_Obj_elemtype( obj ) == FLA_MATRIX )
      {
         dim_t    jj, kk;
         dim_t    m    = FLA_Obj_length( obj );
         dim_t    n    = FLA_Obj_width( obj );
         dim_t    cs   = FLA_Obj_col_stride( obj );
         FLA_Obj* buf  = FLASH_OBJ_PTR_AT( obj );
         
         // Clear each block in macroblock.
         for ( jj = 0; jj < n; jj++ )
            for ( kk = 0; kk < m; kk++ )
               ( buf + jj * cs + kk )->base->write_task = NULL;
      }
      else // Clear regular block.
      {
         obj.base->write_task = NULL;
      }
   }
   
   // Cleaning the last read tasks in each input block.
   for ( i = 0; i < t->n_input_args; i++ )
   {
      obj = t->input_arg[i];

      // Macroblock is used.
      if ( FLA_Obj_elemtype( obj ) == FLA_MATRIX )
      {
         dim_t    jj, kk;
         dim_t    m    = FLA_Obj_length( obj );
         dim_t    n    = FLA_Obj_width( obj );
         dim_t    cs   = FLA_Obj_col_stride( obj );
         FLA_Obj* buf  = FLASH_OBJ_PTR_AT( obj );

         // Clear each block in macroblock.
         for ( jj = 0; jj < n; jj++ )
         {
            for ( kk = 0; kk < m; kk++ )
            {
               obj = *( buf + jj * cs + kk );
               
               k = obj.base->n_read_tasks;
               d = obj.base->read_task_head;

               obj.base->n_read_tasks   = 0;
               obj.base->read_task_head = NULL;
               obj.base->read_task_tail = NULL;

               for ( j = 0; j < k; j++ )
               {
                  next_dep = d->next_dep;
                  FLA_free( d );
                  d = next_dep;
               }
            }
         }
      }      
      else // Regular block.
      {     
         k = obj.base->n_read_tasks;
         d = obj.base->read_task_head;
         
         obj.base->n_read_tasks   = 0;
         obj.base->read_task_head = NULL;
         obj.base->read_task_tail = NULL;
         
         for ( j = 0; j < k; j++ )
         {
            next_dep = d->next_dep;
            FLA_free( d );
            d = next_dep;
         }
      }
   }

   // Free the dep_arg field of t.
   d = t->dep_arg_head;

   for ( i = 0; i < t->n_dep_args; i++ )
   {
      next_dep = d->next_dep;
      FLA_free( d );
      d = next_dep;
   }   

   // Free the int_arg field of t.
   FLA_free( t->int_arg );
   
   // Free the fla_arg field of t.
   FLA_free( t->fla_arg );

   // Free the input_arg field of t.
   FLA_free( t->input_arg );

   // Free the output_arg field of t.
   FLA_free( t->output_arg );

   // Finally, free the struct itself.
   FLA_free( t );

   return;
}


void FLASH_Queue_exec_task( FLASH_Task* t )
/*----------------------------------------------------------------------------

   FLASH_Queue_exec_task

----------------------------------------------------------------------------*/
{
   // Define local function pointer types.

   // LAPACK-level
   typedef FLA_Error(*flash_lu_piv_macro_p)(FLA_Obj A, FLA_Obj p, fla_lu_t* cntl );
   typedef FLA_Error(*flash_apply_pivots_macro_p)(FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl);
   typedef FLA_Error(*flash_lu_piv_p)(FLA_Obj A, FLA_Obj p, fla_lu_t* cntl);
   typedef FLA_Error(*flash_lu_piv_copy_p)(FLA_Obj A, FLA_Obj p, FLA_Obj U, fla_lu_t* cntl);
   typedef FLA_Error(*flash_trsm_piv_p)(FLA_Obj A, FLA_Obj C, FLA_Obj p, fla_trsm_t* cntl);
   typedef FLA_Error(*flash_sa_lu_p)(FLA_Obj U, FLA_Obj D, FLA_Obj p, FLA_Obj L, int nb_alg, fla_lu_t* cntl);
   typedef FLA_Error(*flash_sa_fs_p)(FLA_Obj L, FLA_Obj D, FLA_Obj p, FLA_Obj C, FLA_Obj E, int nb_alg, fla_gemm_t* cntl);
   typedef FLA_Error(*flash_lu_nopiv_p)(FLA_Obj A, fla_lu_t* cntl);
   typedef FLA_Error(*flash_trinv_p)(FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A, fla_trinv_t* cntl);
   typedef FLA_Error(*flash_ttmm_p)(FLA_Uplo uplo, FLA_Obj A, fla_ttmm_t* cntl);
   typedef FLA_Error(*flash_chol_p)(FLA_Uplo uplo, FLA_Obj A, fla_chol_t* cntl);
   typedef FLA_Error(*flash_sylv_p)(FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl);
   typedef FLA_Error(*flash_lyap_p)(FLA_Trans trans, FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale, fla_lyap_t* cntl);
   typedef FLA_Error(*flash_qrut_macro_p)(FLA_Obj A, FLA_Obj T, fla_qrut_t* cntl);
   typedef FLA_Error(*flash_qrut_p)(FLA_Obj A, FLA_Obj T, fla_qrut_t* cntl);
   typedef FLA_Error(*flash_qrutc_p)(FLA_Obj A, FLA_Obj T, FLA_Obj U, fla_qrut_t* cntl);
   typedef FLA_Error(*flash_qr2ut_p)(FLA_Obj B, FLA_Obj D, FLA_Obj T, fla_qr2ut_t* cntl);
   typedef FLA_Error(*flash_lqut_macro_p)(FLA_Obj A, FLA_Obj T, fla_lqut_t* cntl);
   typedef FLA_Error(*flash_caqr2ut_p)(FLA_Obj B, FLA_Obj D, FLA_Obj T, fla_caqr2ut_t* cntl);
   typedef FLA_Error(*flash_uddateut_p)(FLA_Obj R, FLA_Obj C, FLA_Obj D, FLA_Obj T, fla_uddateut_t* cntl);
   typedef FLA_Error(*flash_apqut_p)(FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl);
   typedef FLA_Error(*flash_apq2ut_p)(FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C, FLA_Obj E, fla_apq2ut_t* cntl);
   typedef FLA_Error(*flash_apcaq2ut_p)(FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C, FLA_Obj E, fla_apcaq2ut_t* cntl);
   typedef FLA_Error(*flash_apqudut_p)(FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj T, FLA_Obj W, FLA_Obj R, FLA_Obj U, FLA_Obj C, FLA_Obj V, FLA_Obj D, fla_apqudut_t* cntl);
   typedef FLA_Error(*flash_eig_gest_p)(FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl);

   // Level-3 BLAS
   typedef FLA_Error(*flash_gemm_p)(FLA_Trans transa, FLA_Trans transb, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl);
   typedef FLA_Error(*flash_hemm_p)(FLA_Side side, FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_hemm_t* cntl);
   typedef FLA_Error(*flash_herk_p)(FLA_Uplo uplo, FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C, fla_herk_t* cntl);
   typedef FLA_Error(*flash_her2k_p)(FLA_Uplo uplo, FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_her2k_t* cntl);
   typedef FLA_Error(*flash_symm_p)(FLA_Side side, FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_symm_t* cntl);
   typedef FLA_Error(*flash_syrk_p)(FLA_Uplo uplo, FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C, fla_syrk_t* cntl);
   typedef FLA_Error(*flash_syr2k_p)(FLA_Uplo uplo, FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_syr2k_t* cntl);
   typedef FLA_Error(*flash_trmm_p)(FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj C, fla_trmm_t* cntl);
   typedef FLA_Error(*flash_trsm_p)(FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj C, fla_trsm_t* cntl);

   // Level-2 BLAS
   typedef FLA_Error(*flash_gemv_p)(FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl);
   typedef FLA_Error(*flash_trsv_p)(FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl);

   // Level-1 BLAS
   typedef FLA_Error(*flash_axpy_p)(FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_axpy_t* cntl);
   typedef FLA_Error(*flash_axpyt_p)(FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_axpyt_t* cntl);
   typedef FLA_Error(*flash_copy_p)(FLA_Obj A, FLA_Obj B, fla_copy_t* cntl);
   typedef FLA_Error(*flash_copyt_p)(FLA_Trans trans, FLA_Obj A, FLA_Obj B, fla_copyt_t* cntl);
   typedef FLA_Error(*flash_copyr_p)(FLA_Uplo uplo, FLA_Obj A, FLA_Obj B, fla_copyr_t* cntl);
   typedef FLA_Error(*flash_scal_p)(FLA_Obj alpha, FLA_Obj A, fla_scal_t* cntl);
   typedef FLA_Error(*flash_scalr_p)(FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, fla_scalr_t* cntl);

   // Base
   typedef FLA_Error(*flash_obj_create_buffer_p)(dim_t rs, dim_t cs, FLA_Obj A, void* cntl);
   typedef FLA_Error(*flash_obj_free_buffer_p)(FLA_Obj A, void* cntl);

   // Only execute task if it is not NULL.
   if ( t == NULL )
      return;

   // Now "switch" between the various possible task functions.

   // FLA_LU_piv_macro
   if ( t->func == (void *) FLA_LU_piv_macro_task )
   {
      flash_lu_piv_macro_p func;
      func = (flash_lu_piv_macro_p) t->func;
      
      func(               t->output_arg[0],
                          t->output_arg[1],
            ( fla_lu_t* ) t->cntl );
   }
   // FLA_Apply_pivots_macro
   else if ( t->func == (void *) FLA_Apply_pivots_macro_task )
   {
      flash_apply_pivots_macro_p func;
      func = (flash_apply_pivots_macro_p) t->func;
      
      func( ( FLA_Side  )    t->int_arg[0],
            ( FLA_Trans )    t->int_arg[1],
                             t->input_arg[0],
                             t->output_arg[0],
            ( fla_appiv_t* ) t->cntl );
   }
   // FLA_LU_piv
   else if ( t->func == (void *) FLA_LU_piv_task )
   {
      flash_lu_piv_p func;
      func = (flash_lu_piv_p) t->func;
      
      func(               t->output_arg[0],
                          t->fla_arg[0],
            ( fla_lu_t* ) t->cntl );
   }
   // FLA_LU_piv_copy
   else if ( t->func == (void *) FLA_LU_piv_copy_task )
   {
      flash_lu_piv_copy_p func;
      func = (flash_lu_piv_copy_p) t->func;

      func(               t->output_arg[0],
                          t->fla_arg[0],
                          t->output_arg[1],
            ( fla_lu_t* ) t->cntl );
   }
   // FLA_Trsm_piv
   else if ( t->func == (void *) FLA_Trsm_piv_task )
   {
      flash_trsm_piv_p func;
      func = (flash_trsm_piv_p) t->func;

      func(                 t->input_arg[0],
                            t->output_arg[0],
                            t->fla_arg[0],
            ( fla_trsm_t* ) t->cntl );
   }
   // FLA_SA_LU
   else if ( t->func == (void *) FLA_SA_LU_task )
   {
      flash_sa_lu_p func;
      func = (flash_sa_lu_p) t->func;
      
      func(               t->output_arg[1],
                          t->output_arg[0],
                          t->fla_arg[0],
                          t->fla_arg[1],
                          t->int_arg[0],
            ( fla_lu_t* ) t->cntl );
   }
   // FLA_SA_FS
   else if ( t->func == (void *) FLA_SA_FS_task )
   {
      flash_sa_fs_p func;
      func = (flash_sa_fs_p) t->func;
         
      func(                 t->fla_arg[0],
                            t->input_arg[0],
                            t->fla_arg[1],                          
                            t->output_arg[1],
                            t->output_arg[0],
                            t->int_arg[0],
            ( fla_gemm_t* ) t->cntl );
   }
   // FLA_LU_nopiv
   else if ( t->func == (void *) FLA_LU_nopiv_task )
   {
      flash_lu_nopiv_p func;
      func = (flash_lu_nopiv_p) t->func;

      func(               t->output_arg[0],
            ( fla_lu_t* ) t->cntl );
   }
   // FLA_Trinv
   else if ( t->func == (void *) FLA_Trinv_task )
   {
      flash_trinv_p func;
      func = (flash_trinv_p) t->func;
      
      func( ( FLA_Uplo     ) t->int_arg[0],
            ( FLA_Diag     ) t->int_arg[1],
                             t->output_arg[0],
            ( fla_trinv_t* ) t->cntl );
   }
   // FLA_Ttmm
   else if ( t->func == (void *) FLA_Ttmm_task )
   {   
      flash_ttmm_p func;
      func = (flash_ttmm_p) t->func;
      
      func( ( FLA_Uplo    ) t->int_arg[0],
                            t->output_arg[0],
            ( fla_ttmm_t* ) t->cntl );
   }
   // FLA_Chol
   else if ( t->func == (void *) FLA_Chol_task )
   {   
      flash_chol_p func;
      func = (flash_chol_p) t->func;
      
      func( ( FLA_Uplo    ) t->int_arg[0],
                            t->output_arg[0],
            ( fla_chol_t* ) t->cntl );
   }
   // FLA_Sylv
   else if ( t->func == (void *) FLA_Sylv_task )
   {   
      flash_sylv_p func;
      func = (flash_sylv_p) t->func;
      
      func( ( FLA_Trans   ) t->int_arg[0],
            ( FLA_Trans   ) t->int_arg[1],
                            t->fla_arg[0],
                            t->input_arg[0],
                            t->input_arg[1],
                            t->output_arg[0],
                            t->fla_arg[1],
            ( fla_sylv_t* ) t->cntl );
   }
   // FLA_Lyap
   else if ( t->func == (void *) FLA_Lyap_task )
   {   
      flash_lyap_p func;
      func = (flash_lyap_p) t->func;
      
      func( ( FLA_Trans   ) t->int_arg[0],
                            t->fla_arg[0],
                            t->input_arg[0],
                            t->output_arg[0],
                            t->fla_arg[1],
            ( fla_lyap_t* ) t->cntl );
   }
   // FLA_QR_UT_macro
   else if ( t->func == (void *) FLA_QR_UT_macro_task )
   {   
      flash_qrut_macro_p func;
      func = (flash_qrut_macro_p) t->func;
      
      func(                 t->output_arg[0],
                            t->output_arg[1],
            ( fla_qrut_t* ) t->cntl );
   }
   // FLA_QR_UT
   else if ( t->func == (void *) FLA_QR_UT_task )
   {   
      flash_qrut_p func;
      func = (flash_qrut_p) t->func;
      
      func(                 t->output_arg[0],
                            t->fla_arg[0],
            ( fla_qrut_t* ) t->cntl );
   }
   // FLA_QR_UT_copy
   else if ( t->func == (void *) FLA_QR_UT_copy_task )
   {   
      flash_qrutc_p func;
      func = (flash_qrutc_p) t->func;

      func(                 t->output_arg[0],
                            t->fla_arg[0],
                            t->output_arg[1],
            ( fla_qrut_t* ) t->cntl );
   }
   // FLA_QR2_UT
   else if ( t->func == (void *) FLA_QR2_UT_task )
   {   
      flash_qr2ut_p func;
      func = (flash_qr2ut_p) t->func;
         
      func(                 t->output_arg[1],
                            t->output_arg[0],
                            t->fla_arg[0],
           ( fla_qr2ut_t* ) t->cntl );
   }
   // FLA_LQ_UT_macro
   else if ( t->func == (void *) FLA_LQ_UT_macro_task )
   {   
      flash_lqut_macro_p func;
      func = (flash_lqut_macro_p) t->func;
      
      func(                 t->output_arg[0],
                            t->output_arg[1],
            ( fla_lqut_t* ) t->cntl );
   }
   // FLA_CAQR2_UT
   else if ( t->func == (void *) FLA_CAQR2_UT_task )
   {   
      flash_caqr2ut_p func;
      func = (flash_caqr2ut_p) t->func;
         
      func(                 t->output_arg[1],
                            t->output_arg[0],
                            t->fla_arg[0],
         ( fla_caqr2ut_t* ) t->cntl );
   }
   // FLA_UDdate_UT
   else if ( t->func == (void *) FLA_UDdate_UT_task )
   {   
      flash_uddateut_p func;
      func = (flash_uddateut_p) t->func;
      
      func(                 t->output_arg[0],
                            t->output_arg[1],
                            t->output_arg[2],
                            t->output_arg[3],
        ( fla_uddateut_t* ) t->cntl );
   }
   // FLA_Apply_Q_UT
   else if ( t->func == (void *) FLA_Apply_Q_UT_task )
   {   
      flash_apqut_p func;
      func = (flash_apqut_p) t->func;
      
      func( ( FLA_Side     ) t->int_arg[0],
            ( FLA_Trans    ) t->int_arg[1],
            ( FLA_Direct   ) t->int_arg[2],
            ( FLA_Store    ) t->int_arg[3],
                             t->input_arg[0],
                             t->fla_arg[0],
                             t->output_arg[1],
                             t->output_arg[0],
            ( fla_apqut_t* ) t->cntl );
   }
   // FLA_Apply_Q2_UT
   else if ( t->func == (void *) FLA_Apply_Q2_UT_task )
   {   
      flash_apq2ut_p func;
      func = (flash_apq2ut_p) t->func;
      
      func( ( FLA_Side      ) t->int_arg[0],
            ( FLA_Trans     ) t->int_arg[1],
            ( FLA_Direct    ) t->int_arg[2],
            ( FLA_Store     ) t->int_arg[3],
                              t->input_arg[0],
                              t->fla_arg[0],
                              t->output_arg[2],
                              t->output_arg[1],
                              t->output_arg[0],
            ( fla_apq2ut_t* ) t->cntl );
   }
   // FLA_Apply_CAQ2_UT
   else if ( t->func == (void *) FLA_Apply_CAQ2_UT_task )
   {   
      flash_apcaq2ut_p func;
      func = (flash_apcaq2ut_p) t->func;
      
      func( ( FLA_Side      ) t->int_arg[0],
            ( FLA_Trans     ) t->int_arg[1],
            ( FLA_Direct    ) t->int_arg[2],
            ( FLA_Store     ) t->int_arg[3],
                              t->input_arg[0],
                              t->fla_arg[0],
                              t->output_arg[2],
                              t->output_arg[1],
                              t->output_arg[0],
          ( fla_apcaq2ut_t* ) t->cntl );
   }
   // FLA_Apply_QUD_UT
   else if ( t->func == (void *) FLA_Apply_QUD_UT_task )
   {   
      flash_apqudut_p func;
      func = (flash_apqudut_p) t->func;
      
      func( ( FLA_Side       ) t->int_arg[0],
            ( FLA_Trans      ) t->int_arg[1],
            ( FLA_Direct     ) t->int_arg[2],
            ( FLA_Store      ) t->int_arg[3],
                               t->input_arg[0],
                               t->output_arg[0],
                               t->output_arg[1],
                               t->input_arg[1],
                               t->output_arg[2],
                               t->input_arg[2],
                               t->output_arg[3],
            ( fla_apqudut_t* ) t->cntl );
   }
   // FLA_Eig_gest
   else if ( t->func == (void *) FLA_Eig_gest_task )
   {   
      flash_eig_gest_p func;
      func = (flash_eig_gest_p) t->func;
      
      func( ( FLA_Inv         ) t->int_arg[0],
            ( FLA_Uplo        ) t->int_arg[1],
                                t->output_arg[1],
                                t->output_arg[0],
                                t->input_arg[0],
            ( fla_eig_gest_t* ) t->cntl );
   }
   // FLA_Gemm
   else if ( t->func == (void *) FLA_Gemm_task )
   {
      flash_gemm_p func;
      func = (flash_gemm_p) t->func;
      
      func( ( FLA_Trans   ) t->int_arg[0],
            ( FLA_Trans   ) t->int_arg[1],
                            t->fla_arg[0],
                            t->input_arg[0],
                            t->input_arg[1],
                            t->fla_arg[1],
                            t->output_arg[0],
            ( fla_gemm_t* ) t->cntl );
   }
   // FLA_Hemm
   else if ( t->func == (void *) FLA_Hemm_task )
   {
      flash_hemm_p func;
      func = (flash_hemm_p) t->func;
      
      func( ( FLA_Side    ) t->int_arg[0],
            ( FLA_Uplo    ) t->int_arg[1],
                            t->fla_arg[0],
                            t->input_arg[0],
                            t->input_arg[1],
                            t->fla_arg[1],
                            t->output_arg[0],
            ( fla_hemm_t* ) t->cntl );
   }
   // FLA_Herk
   else if ( t->func == (void *) FLA_Herk_task )
   {
      flash_herk_p func;
      func = (flash_herk_p) t->func;
      
      func( ( FLA_Uplo    ) t->int_arg[0],
            ( FLA_Trans   ) t->int_arg[1],
                            t->fla_arg[0],
                            t->input_arg[0],
                            t->fla_arg[1],
                            t->output_arg[0],
            ( fla_herk_t* ) t->cntl );
   }
   // FLA_Her2k
   else if ( t->func == (void *) FLA_Her2k_task )
   {
      flash_her2k_p func;
      func = (flash_her2k_p) t->func;
      
      func( ( FLA_Uplo     ) t->int_arg[0],
            ( FLA_Trans    ) t->int_arg[1],
                             t->fla_arg[0],
                             t->input_arg[0],
                             t->input_arg[1],
                             t->fla_arg[1],
                             t->output_arg[0],
            ( fla_her2k_t* ) t->cntl );
   }
   // FLA_Symm
   else if ( t->func == (void *) FLA_Symm_task )
   {
      flash_symm_p func;
      func = (flash_symm_p) t->func;
      
      func( ( FLA_Side    ) t->int_arg[0],
            ( FLA_Uplo    ) t->int_arg[1],
                            t->fla_arg[0],
                            t->input_arg[0],
                            t->input_arg[1],
                            t->fla_arg[1],
                            t->output_arg[0],
            ( fla_symm_t* ) t->cntl );
   }
   // FLA_Syrk
   else if ( t->func == (void *) FLA_Syrk_task )
   {
      flash_syrk_p func;
      func = (flash_syrk_p) t->func;
      
      func( ( FLA_Uplo    ) t->int_arg[0],
            ( FLA_Trans   ) t->int_arg[1],
                            t->fla_arg[0],
                            t->input_arg[0],
                            t->fla_arg[1],
                            t->output_arg[0],
            ( fla_syrk_t* ) t->cntl );
   }
   // FLA_Syr2k
   else if ( t->func == (void *) FLA_Syr2k_task )
   {
      flash_syr2k_p func;
      func = (flash_syr2k_p) t->func;
      
      func( ( FLA_Uplo     ) t->int_arg[0],
            ( FLA_Trans    ) t->int_arg[1],
                             t->fla_arg[0],
                             t->input_arg[0],
                             t->input_arg[1],
                             t->fla_arg[1],
                             t->output_arg[0],
            ( fla_syr2k_t* ) t->cntl );
   }
   // FLA_Trmm
   else if ( t->func == (void *) FLA_Trmm_task )
   {
      flash_trmm_p func;
      func = (flash_trmm_p) t->func;
      
      func( ( FLA_Side    ) t->int_arg[0],
            ( FLA_Uplo    ) t->int_arg[1],
            ( FLA_Trans   ) t->int_arg[2],
            ( FLA_Diag    ) t->int_arg[3],
                            t->fla_arg[0],
                            t->input_arg[0],
                            t->output_arg[0],
            ( fla_trmm_t* ) t->cntl );
   }
   // FLA_Trsm
   else if ( t->func == (void *) FLA_Trsm_task )
   {
      flash_trsm_p func;
      func = (flash_trsm_p) t->func;
      
      func( ( FLA_Side    ) t->int_arg[0],
            ( FLA_Uplo    ) t->int_arg[1],
            ( FLA_Trans   ) t->int_arg[2],
            ( FLA_Diag    ) t->int_arg[3],
                            t->fla_arg[0],
                            t->input_arg[0],
                            t->output_arg[0],
            ( fla_trsm_t* ) t->cntl );
   }
   // FLA_Gemv
   else if ( t->func == (void *) FLA_Gemv_task )
   {
      flash_gemv_p func;
      func = (flash_gemv_p) t->func;
      
      func( ( FLA_Trans   ) t->int_arg[0],
                            t->fla_arg[0],
                            t->input_arg[0],
                            t->input_arg[1],
                            t->fla_arg[1],
                            t->output_arg[0],
            ( fla_gemv_t* ) t->cntl );
   }
   // FLA_Trsv
   else if ( t->func == (void *) FLA_Trsv_task )
   {
      flash_trsv_p func;
      func = (flash_trsv_p) t->func;
      
      func( ( FLA_Uplo    ) t->int_arg[0],
            ( FLA_Trans   ) t->int_arg[1],
            ( FLA_Diag    ) t->int_arg[2],
                            t->input_arg[0],
                            t->output_arg[0],
            ( fla_trsv_t* ) t->cntl );
   }
   // FLA_Axpy
   else if ( t->func == (void *) FLA_Axpy_task )
   {
      flash_axpy_p func;
      func = (flash_axpy_p) t->func;
         
      func(                 t->fla_arg[0],
                            t->input_arg[0],
                            t->output_arg[0],
            ( fla_axpy_t* ) t->cntl );
   }
   // FLA_Axpyt
   else if ( t->func == (void *) FLA_Axpyt_task )
   {
      flash_axpyt_p func;
      func = (flash_axpyt_p) t->func;
      
      func( ( FLA_Trans    ) t->int_arg[0],
                             t->fla_arg[0],
                             t->input_arg[0],
                             t->output_arg[0],
            ( fla_axpyt_t* ) t->cntl );
   }
   // FLA_Copy
   else if ( t->func == (void *) FLA_Copy_task )
   {
      flash_copy_p func;
      func = (flash_copy_p) t->func;
         
      func(                 t->input_arg[0],
                            t->output_arg[0],
            ( fla_copy_t* ) t->cntl );
   }
   // FLA_Copyt
   else if ( t->func == (void *) FLA_Copyt_task )
   {
      flash_copyt_p func;
      func = (flash_copyt_p) t->func;
      
      func( ( FLA_Trans    ) t->int_arg[0],
                             t->input_arg[0],
                             t->output_arg[0],
            ( fla_copyt_t* ) t->cntl );
   }
   // FLA_Copyr
   else if ( t->func == (void *) FLA_Copyr_task )
   {
      flash_copyr_p func;
      func = (flash_copyr_p) t->func;
      
      func( ( FLA_Uplo     ) t->int_arg[0],
                             t->input_arg[0],
                             t->output_arg[0],
            ( fla_copyr_t* ) t->cntl );
   }
   // FLA_Scal
   else if ( t->func == (void *) FLA_Scal_task )
   {
      flash_scal_p func;
      func = (flash_scal_p) t->func;

      func(                 t->fla_arg[0],
                            t->output_arg[0],
            ( fla_scal_t* ) t->cntl );
   }
   // FLA_Scalr
   else if ( t->func == (void *) FLA_Scalr_task )
   {
      flash_scalr_p func;
      func = (flash_scalr_p) t->func;

      func( ( FLA_Uplo     ) t->int_arg[0],
                             t->fla_arg[0],
                             t->output_arg[0],
            ( fla_scalr_t* ) t->cntl );
   }
   // FLA_Obj_create_buffer
   else if ( t->func == (void *) FLA_Obj_create_buffer_task )
   {
      flash_obj_create_buffer_p func;
      func = (flash_obj_create_buffer_p) t->func;

      func( ( dim_t       ) t->int_arg[0],
            ( dim_t       ) t->int_arg[1],
                            t->output_arg[0],
                            t->cntl );
   }
   // FLA_Obj_free_buffer
   else if ( t->func == (void *) FLA_Obj_free_buffer_task )
   {
      flash_obj_free_buffer_p func;
      func = (flash_obj_free_buffer_p) t->func;

      func(                 t->output_arg[0],
                            t->cntl );
   }
   else
   {
      FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
   }
   
   return;
}


void FLASH_Queue_verbose_output( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_verbose_output

----------------------------------------------------------------------------*/
{
   int           i, j, k;
   int           n_threads = FLASH_Queue_get_num_threads();
   int           n_tasks   = FLASH_Queue_get_num_tasks();
   FLASH_Verbose verbose   = FLASH_Queue_get_verbose_output();
   FLASH_Task*   t;
   FLASH_Dep*    d;

   // Grab the head of the task queue.
   t = FLASH_Queue_get_head_task();
   
   if ( verbose == FLASH_QUEUE_VERBOSE_READABLE )
   {
      // Iterate over linked list of tasks.
      for ( i = 0; i < n_tasks; i++ )
      {
         printf( "%d\t%s\t", t->order, t->name );
         
         for ( j = 0; j < t->n_output_args; j++ )
            printf( "%lu[%lu,%lu] ", t->output_arg[j].base->id,
                    t->output_arg[j].base->m_index, 
                    t->output_arg[j].base->n_index );
         
         printf( ":= " );
         
         for ( j = 0; j < t->n_output_args; j++ )
            printf( "%lu[%lu,%lu] ", t->output_arg[j].base->id,
                    t->output_arg[j].base->m_index, 
                    t->output_arg[j].base->n_index );
         
         for ( j = 0; j < t->n_input_args; j++ )
            printf( "%lu[%lu,%lu] ", t->input_arg[j].base->id,
                    t->input_arg[j].base->m_index, 
                    t->input_arg[j].base->n_index );

         printf( "\n" );
         
         // Go to the next task.
         t = t->next_task;
      }

      printf( "\n" );
   }
   else
   {
      printf( "digraph SuperMatrix {\n" );

      if ( FLASH_Queue_get_data_affinity() == FLASH_QUEUE_AFFINITY_NONE )
      {
         // Iterate over linked list of tasks.
         for ( i = 0; i < n_tasks; i++ )
         {
            printf( "%d [label=\"%s\"]; %d -> {", t->order, t->name, t->order);
            
            d = t->dep_arg_head;
            for ( j = 0; j < t->n_dep_args; j++ )
            {
               printf( "%d;", d->task->order );
               d = d->next_dep;
            }

            printf( "};\n" );
            
            // Go to the next task.
            t = t->next_task;
         }
      }
      else
      {
         // Iterate over all the threads.
         for ( k = 0; k < n_threads; k++ )
         {
            printf( "subgraph cluster%d {\nlabel=\"%d\"\n", k, k );

            // Iterate over linked list of tasks.
            for ( i = 0; i < n_tasks; i++ )
            {  
               if ( t->queue == k )
                  printf( "%d [label=\"%s\"];\n", t->order, t->name );
               
               // Go to the next task.
               t = t->next_task;               
            }

            printf( "}\n" );

            // Grab the head of the task queue.
            t = FLASH_Queue_get_head_task();               
         }

         // Iterate over linked list of tasks.
         for ( i = 0; i < n_tasks; i++ )
         {
            printf( "%d -> {", t->order );
            
            d = t->dep_arg_head;
            for ( j = 0; j < t->n_dep_args; j++ )
            {
               printf( "%d;", d->task->order );
               d = d->next_dep;
            }

            printf( "};\n" );
            
            // Go to the next task.
            t = t->next_task;
         }
      }

      printf( "}\n\n" );
   }
   
   return;
}


#endif // FLA_ENABLE_SUPERMATRIX
