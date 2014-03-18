/*
 Here is the usage for FLAME workqueuing:

   FLA_Queue_init();

   FLA_Part( ... );

   while( ... )
   {
     FLA_Repart( ... );

     ENQUEUE_FLA_Gemm( ... );

     FLA_Cont_with( ... );
   }

   FLA_Queue_exec();

   FLA_Queue_finalize();
*/


#include "FLAME.h"

#if FLA_MULTITHREADING_MODEL == FLA_OPENMP

#include "FLA_queue_omp.h"
#include "FLA_task_partitioning.h"
#include <assert.h>
#include <stdarg.h>
#include <omp.h>



FLA_Queue  tq;
static int fla_queue_n_threads   = 1;
static int fla_queue_initialized = 0;


void FLA_Queue_exec()
{
	// If the queue is empty, return early.
	if( tq.n_tasks == 0 )
		return;
	
	//FLA_queue_print_costs();
	FLA_queue_exec_sync();
}

void FLA_queue_print_costs()
{
	int i;
	
	for( i = 0; i < tq.n_tasks; ++i )
	{
		fprintf( stdout, "cost of task %2d: %e\n", i, tq.task_array[i]->cost ); fflush(stdout);
	}
}   


void FLA_queue_exec_sync()
{
	int        i;
	int        n_tasks;
	FLA_Task** task_array;
	FLA_Task*  t;
	
	// The queue is full, so we may now create an array index of each task.
	FLA_queue_create_task_array();
	
	// Sort the array if it is not sorted already
	FLA_queue_sort_task_array();
	//FLA_queue_print_costs();

	// Copy the task_array pointer and n_tasks integer locally to make the
	// OpenMP compiler happy.
	n_tasks    = tq.n_tasks;
	task_array = tq.task_array;
	
	// Iterate over the task queue using the random-access array.
	#pragma omp parallel for \
	        shared( task_array, n_tasks ) \
	        private( i, t ) \
	        schedule( dynamic, 1 )
	for( i = n_tasks - 1; i >= 0; --i )
	//for( i = 0; i < n_tasks; ++i )
	{
		t = task_array[i];
		FLA_queue_exec_task( t );
	}
	
	// Flush the queue. To do this, we walk the task_array and free() each
	// element.
	FLA_queue_flush();

	// Now that we're done with the task array, we can free it. 
	FLA_queue_free_task_array();
}

void FLA_queue_exec_task( FLA_Task* t )
{
	// Define local function pointer types.
	typedef FLA_Error(*fla_gemm_p)(FLA_Trans transa, FLA_Trans transb, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C);
	typedef FLA_Error(*fla_symm_p)(FLA_Side side, FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C);
	typedef FLA_Error(*fla_syrk_p)(FLA_Uplo uplo, FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C);
	typedef FLA_Error(*fla_syr2k_p)(FLA_Uplo uplo, FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C);
	typedef FLA_Error(*fla_trmm_p)(FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj C);
	typedef FLA_Error(*fla_trsm_p)(FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj C);

	
	// Now "switch" between the various possible task functions.
	
	// FLA_Gemm
	if     ( t->func == (void*)FLA_Gemm_external )
	{
		fla_gemm_p func;
		
		// Here we unpack our void* function pointer and typecast it appropriately.
		func = (fla_gemm_p) t->func;

		// Invoke FLA_Gemm()
		func( ( FLA_Trans ) t->int_arg[0],
		      ( FLA_Trans ) t->int_arg[1],
		      ( FLA_Obj   ) t->fla_arg[0],
		      ( FLA_Obj   ) t->fla_arg[1],
		      ( FLA_Obj   ) t->fla_arg[2],
		      ( FLA_Obj   ) t->fla_arg[3],
		      ( FLA_Obj   ) t->fla_arg[4] );
	}
	
	// FLA_Symm
	else if( t->func == (void*)FLA_Symm_external )
	{
		fla_symm_p func;
		
		// Here we unpack our void* function pointer and typecast it appropriately.
		func = (fla_symm_p) t->func;

		// Invoke FLA_Symm()
		func( ( FLA_Side  ) t->int_arg[0],
		      ( FLA_Uplo  ) t->int_arg[1],
		      ( FLA_Obj   ) t->fla_arg[0],
		      ( FLA_Obj   ) t->fla_arg[1],
		      ( FLA_Obj   ) t->fla_arg[2],
		      ( FLA_Obj   ) t->fla_arg[3],
		      ( FLA_Obj   ) t->fla_arg[4] );
	}

	// FLA_Syrk
	else if( t->func == (void*)FLA_Syrk_external )
	{
		fla_syrk_p func;
		
		// Here we unpack our void* function pointer and typecast it appropriately.
		func = (fla_syrk_p) t->func;

		// Invoke FLA_Syrk()
		func( ( FLA_Uplo  ) t->int_arg[0],
		      ( FLA_Trans ) t->int_arg[1],
		      ( FLA_Obj   ) t->fla_arg[0],
		      ( FLA_Obj   ) t->fla_arg[1],
		      ( FLA_Obj   ) t->fla_arg[2],
		      ( FLA_Obj   ) t->fla_arg[3] );
	}

	// FLA_Syr2k
	else if( t->func == (void*)FLA_Syr2k_external )
	{
		fla_syr2k_p func;
		
		// Here we unpack our void* function pointer and typecast it appropriately.
		func = (fla_syr2k_p) t->func;

		// Invoke FLA_Syr2k()
		func( ( FLA_Uplo  ) t->int_arg[0],
		      ( FLA_Trans ) t->int_arg[1],
		      ( FLA_Obj   ) t->fla_arg[0],
		      ( FLA_Obj   ) t->fla_arg[1],
		      ( FLA_Obj   ) t->fla_arg[2],
		      ( FLA_Obj   ) t->fla_arg[3],
		      ( FLA_Obj   ) t->fla_arg[4] );
	}

	// FLA_Trmm
	else if( t->func == (void*)FLA_Trmm_external )
	{
		fla_trmm_p func;
		
		// Here we unpack our void* function pointer and typecast it appropriately.
		func = (fla_trmm_p) t->func;

		// Invoke FLA_Trmm()
		func( ( FLA_Side  ) t->int_arg[0],
		      ( FLA_Uplo  ) t->int_arg[1],
		      ( FLA_Trans ) t->int_arg[2],
		      ( FLA_Diag  ) t->int_arg[3],
		      ( FLA_Obj   ) t->fla_arg[0],
		      ( FLA_Obj   ) t->fla_arg[1],
		      ( FLA_Obj   ) t->fla_arg[2] );
	}

	// FLA_Trsm
	else if( t->func == (void*)FLA_Trsm_external )
	{
		fla_trsm_p func;
		
		// Here we unpack our void* function pointer and typecast it appropriately.
		func = (fla_trsm_p) t->func;

		// Invoke FLA_Trsm()
		func( ( FLA_Side  ) t->int_arg[0],
		      ( FLA_Uplo  ) t->int_arg[1],
		      ( FLA_Trans ) t->int_arg[2],
		      ( FLA_Diag  ) t->int_arg[3],
		      ( FLA_Obj   ) t->fla_arg[0],
		      ( FLA_Obj   ) t->fla_arg[1],
		      ( FLA_Obj   ) t->fla_arg[2] );
	}

	else
	{
        FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}
}

void FLA_Queue_init()
{
	// Exit early if we're already initialized.
	if( fla_queue_initialized == 1 )
		return;
	
	// Initialze the basic fields of FLA_Queue.
	tq.n_tasks = 0;
	tq.head    = NULL;
	tq.tail    = NULL;
	
	// Set the initialized flag.
	fla_queue_initialized = 1;
}

void FLA_Queue_finalize()
{
	// Exit early if we're not already initialized.
	if( fla_queue_initialized == 0 )
		return;
	
	// Make sure the queue is empty.
	assert( tq.n_tasks == 0 );
	
	// Clear the initialized flag.
	fla_queue_initialized = 0;
}

void FLA_Queue_set_num_threads( int n_threads )
{
	char env_var_str[32];
	char nth_str[32];
	int  overwrite = 1;
	
	// Set the number of threads using the OpenMP interface.
	omp_set_num_threads( n_threads );
	
	// Also keep that value locally.
	fla_queue_n_threads = n_threads;
	
	// Also, try to set the OMP_NUM_THREADS environment variable.
	//sprintf( nth_str, "OMP_NUM_THREADS=%d", n_threads );
	//putenv( nth_str );
	sprintf( env_var_str, "OMP_NUM_THREADS" );
	sprintf( nth_str, "%d", n_threads );
	setenv( env_var_str, nth_str, overwrite );
}

int FLA_Queue_get_num_threads()
{
	return fla_queue_n_threads;
}

void FLA_Queue_push( void* func, double cost, int n_int_params, int n_fla_params, ... )
{
	FLA_Task* t;
	int       i;
	va_list   var_arg_list;
	
	// Make sure we've initialized the queue
	if( fla_queue_initialized == 0 )
		FLA_Queue_init();
	
	// Allocate a new FLA_Task and populate its fields with appropriate values.
	t = FLA_task_alloc_init( func, cost, n_int_params, n_fla_params );
	
	// Initialize variable argument environment. In case you're wondering, the
	// second argument in this macro invocation of va_start() is supposed to be
	// the parameter that immediately preceeds the variable argument list 
	// (ie: the ... above ).
	va_start( var_arg_list, n_fla_params );
	
	// Extract the integer arguments.
	for( i = 0; i < n_int_params; ++i )
	{
		t->int_arg[i] = va_arg( var_arg_list, int );
	}

	// Extract the FLA_Obj arguments.
	for( i = 0; i < n_fla_params; ++i )
	{
		t->fla_arg[i] = va_arg( var_arg_list, FLA_Obj );
	}
	
	// Finalize the variable argument environment.
	va_end( var_arg_list );
	
	// DON'T USE _push() if you need the task_array()! It doesn't exist yet!
	// Push a pointer to the task structure onto the task queue.
	FLA_queue_push_unsorted( t );
}

void FLA_queue_push_unsorted( FLA_Task* t )
{
	// Add the task to the tail of the queue (and the head if the queue is empty).
	if( tq.n_tasks == 0 )
	{
		tq.head = t;
		tq.tail = t;
	}
	else
	{
		tq.tail->next_task = t;
		tq.tail            = t;
	}
	
	// Increment the number of tasks.
	++tq.n_tasks;
}

void FLA_queue_insert_sorted( FLA_Task* t )
{
	FLA_Task* t_curr;
	FLA_Task* t_prev;
	
	// If the queue is empty, then add the task to the head of the queue and return
	// early. Notice that we don't need to maintain a tail for a sorted queue.
	if( tq.n_tasks == 0 )
	{
		tq.head = t;
	}
	// If the cost of the incoming task is greater than the cost of the head task,
	// then insert the task at the head of the queue, and return.
	else if( t->cost >= tq.head->cost )
	{
		t->next_task = tq.head;
		tq.head      = t;
	}
	else
	{
		// Initialize our pointers. Notice that we have to keep track of the current
		// task as well as the previous task since the linked list is only linked in 
		// one direction.
		t_prev = tq.head;
		t_curr = tq.head->next_task;
		
		// Iterate through the task queue. The only way to terminate is through inserting
		// the task, either somewhere in the middle, or at the tail of the linked list.
		while( 1 )
		{
			// If we reach the end of the linked list, then we insert the task at
			// the end and break out of the while loop.
			if( t_curr == NULL )
			{
				t_prev->next_task = t;
				break;
			}
			
			// If the cost of the incoming task is greater than the cost of the current
			// task, then insert the task before the current task (and after the previous
			// task), and break out of the while loop.
			if( t->cost >= t_curr->cost )
			{
				t->next_task      = t_curr;
				t_prev->next_task = t;
				break;
			}
			
			// Advance our pointers.
			t_prev = t_curr;
			t_curr = t_curr->next_task;
		}
	}
	
	// Increment the number of tasks.
	++tq.n_tasks;
}

void FLA_queue_create_task_array()
{
	int       i;
	FLA_Task* t_curr;
	
	// Allocate memory for the array of FLA_Task pointers.
	tq.task_array = (FLA_Task**)malloc( tq.n_tasks * sizeof(FLA_Task*) );
	assert( tq.task_array != NULL );
	
	// Initialize index and FLA_Task* variables.
	t_curr = tq.head;
	i      = 0;
	
	// Iterate over linked list of tasks.
	while( i < tq.n_tasks )
	{
		// Copy the current task pointer to the task array.
		tq.task_array[i] = t_curr;
		
		// Increment task array index and advance FLA_Task pointer.
		t_curr = t_curr->next_task;
		++i;
	}
}


void FLA_queue_sort_task_array()
{
	void*  task_array;
	size_t task_array_length;
	size_t task_pointer_size;

	// Typecast the task array to a void pointer. Notice that the task array
	// must already be created or else we're sorting a NULL array!
	task_array = (void*)tq.task_array;
	
	// Set the length of the task array and the size of each element.
	task_array_length = tq.n_tasks;
	task_pointer_size = sizeof(FLA_Task*);

	// Invoke quicksort from the standard C library.
	qsort( task_array, task_array_length, task_pointer_size, 
	       FLA_queue_task_cost_compare );
}

int  FLA_queue_task_cost_compare( const void* t0, const void* t1 )
{
	FLA_Task* task0;
	FLA_Task* task1;
	double    diff;
	int       r_val;
	
	task0 = *((FLA_Task**)t0);
	task1 = *((FLA_Task**)t1);
	diff  = task0->cost - task1->cost;
	
	r_val = (int)( diff/fabs(diff) );
	
	return r_val;
}

void FLA_queue_free_task_array()
{
	// Assert that the task array exists.
	assert( tq.task_array != NULL );

	// Free the task array.
	free( tq.task_array );
	
	// Clear the task_array field.
	tq.task_array = NULL;
}

void FLA_queue_flush()
{
	int i;
	
	for( i = 0; i < tq.n_tasks; ++i )
	{
		// Free the current task.
		FLA_task_free( tq.task_array[i] );
	}
	
	// Clear the other fields of the FLA_Queue structure.
	tq.n_tasks = 0;
	tq.head    = NULL;
	tq.tail    = NULL;
}

FLA_Task* FLA_task_alloc_init( void* func, double cost, int n_int_args, int n_fla_args )
{
	FLA_Task* t;
	
	// Allocate space for the task structure t.
	t             = (FLA_Task*)malloc( sizeof(FLA_Task) );
	assert( t != NULL );
	
	// Allocate space for the task's integer arguments.
	t->int_arg    =      (int*)malloc( n_int_args * sizeof(int) );
	assert( t->int_arg != NULL );
	
	// Allocate space for the task's FLA_Obj arguments.
	t->fla_arg    =  (FLA_Obj*)malloc( n_fla_args * sizeof(FLA_Obj) );
	assert( t->fla_arg != NULL );
	
	// Initialize other fields of the structure.
	t->cost       = cost;
	t->func       = func;
	t->n_int_args = n_int_args;
	t->n_fla_args = n_fla_args;
	t->next_task  = NULL;
	
	// Return a pointer to the initialized structure.
	return t;
}

void FLA_task_free( FLA_Task* t )
{
	// Free the int_arg field of t.
	free( t->int_arg );

	// Free the fla_arg field of t.
	free( t->fla_arg );
	
	// Finally, free the struct itself.
	free( t );
}

#endif
