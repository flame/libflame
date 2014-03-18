
/*

 Here is the usage for FLAME workqueuing:



 FLA_Queue_init();                -- optional

 FLA_Part( ... );

 while( ... )
 {
	FLA_Repart( ... );

	ENQUEUE_FLA_Gemm( ... );
	
	FLA_Cont_with( ... );
 }

 FLA_Queue_exec();

 FLA_Queue_finalize();            -- possibly optional

*/


#include "FLAME.h"
#include "FLA_queue.h"
#include <assert.h>
#include <pthread.h>





FLA_Queue  q;
int        fla_queue_n_threads = 1;
static int fla_queue_initialized = 0;




void FLA_Queue_exec( int mode )
{
	FLA_queue_exec_sync();
}

void FLA_queue_exec_sync()
{
	int i, r_val, n_threads;
	
	
	// Inquire the number of threads we are to use to dequeue
	// and execute tasks from the queue.
	n_threads = FLA_Queue_get_num_threads();
	
	// Spawn n_threads - 1 pthreads.
	for( i = 1; i < n_threads; ++i )
	{
		r_val = pthread_create( &q.thread[i], &q.thread_attr, FLA_queue_exec_thread_sync, (void*)i );
		assert( r_val == 0 );
	}
	
	// Now put the main thread to work.
	FLA_queue_exec_thread_sync( (void*)0 );
	
	// Synchronize. Join all pthreads.
	for( i = 1; i < n_threads; ++i )
	{
		//r_val = pthread_join( q.thread[i], (void**)&status );
		r_val = pthread_join( q.thread[i], NULL );
		assert( r_val == 0 );
	}
	
	// Cleanup the FLA_Queue.
	//FLA_Queue_finalize();
	
}

void* FLA_queue_exec_thread_sync( void* p )
{
	FLA_Task* t;
	
	// Attempt to pop a task off the queue.
	t = FLA_queue_pop();
	
	// Continue popping and executing tasks off the queue until
	// the queue is empty.
	while( t != NULL )
	{
		// Execute the current task.
		FLA_queue_exec_task( t );

		// Free the FLA_Task object now that the task has been executed.
		FLA_task_free( t );
		
		// Attempt to pop the next task off the queue. If this fails,
		// Then the while loop conditional will fail, and the thread
		// returns (and exits).
		t = FLA_queue_pop( t );
	} 
	
	return NULL;
}

void FLA_queue_exec_task( FLA_Task* t )
{
	unsigned long p_fla_gemm = (unsigned long)FLA_Gemm;

	if( (unsigned long)t->func == p_fla_gemm )
	{
		int (*func)(int transa, int transb, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C);
		
		// Here we unpack our void* function pointer and typecast it appropriately.
		func = t->func;

		// Invoke FLA_Gemm()
		func( t->int_arg[0], t->int_arg[1], t->fla_arg[0], t->fla_arg[1], t->fla_arg[2], t->fla_arg[3], t->fla_arg[4] );
	}
	else
	{
		FLA_Abort( "Unimplemented task function!", __LINE__, __FILE__ );
		exit(1);
	}
	
/*
	// Switch amongst the types of tasks anticipated.
	switch( (unsigned long)t->func )
	{
		//case FLA_Gemm:
		case p_fla_gemm:
		{
			int (*func)(int transa, int transb, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C);
			
			// Here we unpack our void* function pointer and typecast it appropriately.
			func = t->func;
			
			// Invoke FLA_Gemm()
			func( t->int_arg[0], t->int_arg[1], t->fla_arg[0], t->fla_arg[1], t->fla_arg[2], t->fla_arg[3], t->fla_arg[4] );
		}
		break;
		
		default:
			FLA_Abort( "Unimplemented task function!", __LINE__, __FILE__ );
			exit(1);
	}
*/
}



void FLA_Queue_init()
{
	int r_val, n_threads;
	
	
	// Exit early if we're already initialized.
	if( fla_queue_initialized == 1 )
		return;
	
	// Initialze the basic fields of FLA_Queue.
	q.n_tasks = 0;
	q.head    = NULL;
	q.tail    = NULL;
	
	// Inquire the number of threads we are to use to dequeue
	// and execute tasks from the queue.
	n_threads = FLA_Queue_get_num_threads();
	
	// Initialize the pthread mutex object. We request default mutex attributes
	// by passing NULL for the mutex attribute object.
	r_val = pthread_mutex_init( &q.lock, NULL );
	assert( r_val == 0 );
	
	// Initialize the pthread attribute object.
	r_val = pthread_attr_init( &q.thread_attr );
	assert( r_val == 0 );
	
	// Set the thread attribute scope to system-wide.
	r_val = pthread_attr_setscope( &q.thread_attr, PTHREAD_SCOPE_SYSTEM );
	assert( r_val == 0 );
	
	// Allocate memory for the pthread objects.
	q.thread = (pthread_t*)malloc( n_threads * sizeof(pthread_t) );
	assert( q.thread );
	
	// Set the initialized flag.
	fla_queue_initialized = 1;
}


void FLA_Queue_finalize()
{
	int r_val;
	
	
	// Exit early if we're not already initialized.
	if( fla_queue_initialized == 0 )
		return;
	
	// Make sure the queue is empty.
	assert( q.n_tasks == 0 );
	
	// Free the pthread objects.
	free( q.thread );
	
	// Free the pthread attribute object.
	r_val = pthread_attr_destroy( &q.thread_attr );
	assert( r_val == 0 );
	
	// Free the pthread mutex object that protected the queue.
	r_val = pthread_mutex_destroy( &q.lock );
	assert( r_val == 0 );
	
	// Clear the initialized flag.
	fla_queue_initialized = 0;
}


void FLA_Queue_set_num_threads( int n_threads )
{
	fla_queue_n_threads = n_threads;
}

int FLA_Queue_get_num_threads()
{
	return fla_queue_n_threads;
}


void FLA_Queue_push( void* func, int cost, int n_int_params, int n_fla_params, int param0, int param1,
                    FLA_Obj param2, FLA_Obj param3, FLA_Obj param4, FLA_Obj param5, FLA_Obj param6 )
{
	FLA_Task* t;
	
	
	// Make sure we've initialized the queue
	if( fla_queue_initialized == 0 )
		FLA_Queue_init();
	
	// Allocate a new FLA_Task and populate its fields with appropriate values.
	t = FLA_task_alloc_init( func, cost, n_int_params, n_fla_params );
	
	// Insert integer arguments into the task structure.
	t->int_arg[0] = param0;
	t->int_arg[1] = param1;
	
	// Insert FLA_Obj arguments into the task structure.
	t->fla_arg[0] = param2;
	t->fla_arg[1] = param3;
	t->fla_arg[2] = param4;
	t->fla_arg[3] = param5;
	t->fla_arg[4] = param6;
	
	// Push a pointer to the task structure onto the queue.
	// Notice that this does not maintain any sorted order.
	FLA_queue_push( t );
	
	//fprintf( stderr, "pushed\n" ); fflush(stderr);
}


void FLA_queue_push( FLA_Task* t )
{
	
	// Add the task to the tail of the queue (and the head if the queue is empty).
	if( q.n_tasks == 0 )
	{
		q.head = t;
		q.tail = t;
	}
	else
	{
		q.tail->next_task = t;
		q.tail            = t;
	}
	
	// Increment the number of tasks.
	q.n_tasks++;
}

void FLA_queue_flush()
{
	int i;
	FLA_Task* t;
	FLA_Task* t_next_task;
	
	// Start with the head of the queue. If the queue is empty, then
	// the while loop will not execute.
	t = q.head;
	
	while( t != NULL )
	{
		// Keep track of the next_task field, because we're about to free the whole
		// task and would otherwise loose track of it.
		t_next_task = t->next_task;
		
		// Free the current task t.
		FLA_task_free( t );
		
		// Make the next_task the current task. If t_next is NULL, then we're done.
		t = t_next_task;
	}
	
	// Clear the other fields of the FLA_Queue structure.
	q.n_tasks = 0;
	q.head    = NULL;
	q.tail    = NULL;
}

FLA_Task* FLA_queue_pop()
{
	FLA_Task* t;
	int       r_val;
	
	// Lock the mutex so we may enter the critical section.
	r_val = pthread_mutex_lock( &q.lock );
	assert( r_val == 0 );
	
	// Simply return NULL if queue is empty.
	if( q.n_tasks == 0 )
	{
		t = NULL;
	}
	else
	{
		// Get the pointer to the head task on the queue.
		t = q.head;
		
		// Change the head pointer to reflect popping the task.
		q.head = q.head->next_task;
		
		// Decrement the number of tasks.
		q.n_tasks--;
		
		// If the queue is now empty, then nullify the tail pointer.
		if( q.n_tasks == 0 )
		{
			q.tail = NULL;
		}
		
		// Nullify the next_task field of the task we just removed from the queue.
		t->next_task = NULL;
	}
	
	// Unlock the mutex. We're done accessing the critical data.
	r_val = pthread_mutex_unlock( &q.lock );
	assert( r_val == 0 );
	
	// Return a pointer to the separated task (or NULL if the queue was empty).
	return t;
}


FLA_Task* FLA_task_alloc_init( void* func, int cost, int n_int_args, int n_fla_args )
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

