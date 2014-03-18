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
#include "FLA_queue.h"
#include <assert.h>
#include <pthread.h>



static FLA_Queue        tq;
static int              fla_queue_n_threads           = 1;
static int              fla_queue_n_threads_awake;
static int              fla_queue_initialized         = FALSE;
static int              fla_queue_threads_should_live;
static pthread_mutex_t  mutex_of_cond_vars;
static pthread_cond_t   cond_var_wake_up_threads, 
                        cond_var_last_thread_to_sleep;


void FLA_Queue_exec()
{
	FLA_queue_exec_sync();
}

void FLA_queue_exec_sync()
{
	int n_threads;
	
	// Inquire the number of threads we are to use to dequeue
	// and execute tasks from the queue.
	n_threads = FLA_Queue_get_num_threads();
	
	// Wake up the worker threads now that tasks exist on the task queue.
	if( n_threads > 1 )
		FLA_queue_threads_awaken();
	
	// The pthreads now begin to dequeue and execute tasks from the task queue...
	
	// We want the main thread to do its share of work since we only spawned 
	// (n_threads - 1) additional pthreads.
	FLA_queue_engage_tasks();
	
	// Now that the main thread has detected an empty queue, we have to wait
	// for all other threads to signal that they are going back to sleep.
	if( n_threads > 1 )
		FLA_queue_threads_sleep_barrier();
	
}

void FLA_queue_threads_awaken()
{
	int r_val;
	
	// Lock the mutex.
	r_val = pthread_mutex_lock( &mutex_of_cond_vars );
	assert( r_val == 0 );
	
	// Send a "wake up" signal to all threads.
	r_val = pthread_cond_broadcast( &cond_var_wake_up_threads );
	assert( r_val == 0 );

	// Unlock the mutex.
	r_val = pthread_mutex_unlock( &mutex_of_cond_vars );
	assert( r_val == 0 );
}

void FLA_queue_threads_sleep_barrier()
{
	int r_val;
	
	// Lock the mutex associated with the condition variable.
	r_val = pthread_mutex_lock( &mutex_of_cond_vars );
	assert( r_val == 0 );
	
	// Wait for all threads to fall asleep.
	//if( fla_queue_n_threads_awake > 0 )
	while( fla_queue_n_threads_awake > 0 )
	{
		r_val = pthread_cond_wait( &cond_var_last_thread_to_sleep, &mutex_of_cond_vars );
		assert( r_val == 0 );
	}
	
	// When all threads are asleep, unlock the mutex, return, and exit.
	r_val = pthread_mutex_unlock( &mutex_of_cond_vars );
	assert( r_val == 0 );
}

void FLA_queue_threads_sleep_init()
{
	int r_val;
	
	// Lock the mutex.
	r_val = pthread_mutex_lock( &mutex_of_cond_vars );
	assert( r_val == 0 );
	
	// Register that the thread is falling asleep.
	--fla_queue_n_threads_awake;
	
	// Wait for the "wake up" signal.
	r_val = pthread_cond_wait( &cond_var_wake_up_threads, &mutex_of_cond_vars );
	assert( r_val == 0 );

	// Register that the thread is awake.
	++fla_queue_n_threads_awake;

	// Unlock the mutex.
	r_val = pthread_mutex_unlock( &mutex_of_cond_vars );
	assert( r_val == 0 );
}

void FLA_queue_threads_sleep()
{
	int r_val;
	
	// Lock the mutex.
	r_val = pthread_mutex_lock( &mutex_of_cond_vars );
	assert( r_val == 0 );
	
	// Register that the thread is falling asleep.
	--fla_queue_n_threads_awake;
	
	// Send signal to main thread: last pthread is about to fall asleep.
	if( fla_queue_n_threads_awake == 0 )
	{
		r_val = pthread_cond_signal( &cond_var_last_thread_to_sleep );
		assert( r_val == 0 );
	}
	
	// Wait for the "wake up" signal.
	r_val = pthread_cond_wait( &cond_var_wake_up_threads, &mutex_of_cond_vars );
	assert( r_val == 0 );

	// Register that the thread is awake.
	++fla_queue_n_threads_awake;

	// Unlock the mutex.
	r_val = pthread_mutex_unlock( &mutex_of_cond_vars );
	assert( r_val == 0 );
}

void* FLA_queue_exec_thread( void* p )
{
	// Fall asleep until we get a signal from the main thread that the queue
	// is ready. This function is different from FLA_queue_threads_sleep() in
	// that it does not signal the main thread when the last thread is falling
	// asleep.
	FLA_queue_threads_sleep_init();
	
	// Loop until the main thread calls FLA_Queue_finalize(). There, the flag
	// fla_queue_threads_should_live is set to FALSE and the threads are woken
	// up. (This global variable does not need a lock since only the main thread
	// writes to it, while the worker pthreads only read from it.)
	while( fla_queue_threads_should_live == TRUE )
	{
		// Dequeue and execute tasks from the queue until the queue is empty.
		FLA_queue_engage_tasks();
	
		// Fall back asleep and wait for more tasks.
		FLA_queue_threads_sleep();
	}
	
	// Return a value to make the compiler happy.
	return NULL;
}

void FLA_queue_engage_tasks()
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
	int i, r_val, n_threads;
	
	// Exit early if we're already initialized.
	if( fla_queue_initialized == 1 )
		return;
	
	// Initialze the basic fields of FLA_Queue.
	tq.n_tasks = 0;
	tq.head    = NULL;
	tq.tail    = NULL;
	
	// Inquire the number of threads we are to use to dequeue
	// and execute tasks from the queue.
	n_threads = FLA_Queue_get_num_threads();
	
	// Initialize the pthread mutex object. We request default mutex attributes
	// by passing NULL for the mutex attribute object.
	r_val = pthread_mutex_init( &tq.lock, NULL );
	assert( r_val == 0 );
	
	// Initialize the pthread attribute object.
	r_val = pthread_attr_init( &tq.thread_attr );
	assert( r_val == 0 );
	
	// Set the thread attribute scope to system-wide.
	r_val = pthread_attr_setscope( &tq.thread_attr, PTHREAD_SCOPE_SYSTEM );
	assert( r_val == 0 );
	
	// Initialize the condition variable mutex object.
	r_val = pthread_mutex_init( &mutex_of_cond_vars, NULL );
	assert( r_val == 0 );
	
	// Initialize the condition variable object.
	r_val = pthread_cond_init( &cond_var_wake_up_threads, NULL );
	assert( r_val == 0 );
	
	// Initialize the condition variable object.
	r_val = pthread_cond_init( &cond_var_last_thread_to_sleep, NULL );
	assert( r_val == 0 );
	
	// Allocate memory for the pthread objects.
	tq.thread = (pthread_t*)malloc( n_threads * sizeof(pthread_t) );
	assert( tq.thread );
	
	// Start with (n_threads - 1) pthreads awake. They'll fall asleep almost immediately, though.
	fla_queue_n_threads_awake = n_threads - 1;
	
	// Indicate that the threads should live until told to die.
	fla_queue_threads_should_live = TRUE;
	
	// Spawn n_threads - 1 pthreads.
	for( i = 1; i < n_threads; ++i )
	{
		r_val = pthread_create( &tq.thread[i], &tq.thread_attr, FLA_queue_exec_thread, (void*)i );
		assert( r_val == 0 );
	}
	
	// Set the initialized flag.
	fla_queue_initialized = 1;
}


void FLA_Queue_finalize()
{
	int i, r_val, n_threads;
	
	// Exit early if we're not already initialized.
	if( fla_queue_initialized == 0 )
		return;
	
	// Indicate that the threads, upon waking up, should exit gracefully rather
	// than attempt to dequeue and execute tasks from the task queue.
	fla_queue_threads_should_live = FALSE;
	
	// Wake up the threads. Time to go home.
	FLA_queue_threads_awaken();
	
	// Inquire the number of threads we are to use to dequeue
	// and execute tasks from the queue.
	n_threads = FLA_Queue_get_num_threads();
	
	// Now we Synchronize. Join all pthreads, just to be safe. (Probably not needed.)
	for( i = 1; i < n_threads; ++i )
	{
		//r_val = pthread_join( tq.thread[i], (void**)&status );
		r_val = pthread_join( tq.thread[i], NULL );
		assert( r_val == 0 );
	}
	
	// Make sure the queue is empty.
	assert( tq.n_tasks == 0 );
	
	// Free the pthread objects.
	free( tq.thread );
	
	// Free the condition variable mutex object.
	r_val = pthread_mutex_destroy( &mutex_of_cond_vars );
	assert( r_val == 0 );
	
	// Free the condition variable object.
	r_val = pthread_cond_destroy( &cond_var_wake_up_threads );
	assert( r_val == 0 );
	
	// Free the condition variable object.
	r_val = pthread_cond_destroy( &cond_var_last_thread_to_sleep );
	assert( r_val == 0 );
	
	// Free the pthread attribute object.
	r_val = pthread_attr_destroy( &tq.thread_attr );
	assert( r_val == 0 );
	
	// Free the pthread mutex object that protected the queue.
	r_val = pthread_mutex_destroy( &tq.lock );
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
	
	//fprintf( stderr, "initing %p %d %d %d\n", func, cost, n_int_params, n_fla_params ); fflush(stderr);
	// Allocate a new FLA_Task and populate its fields with appropriate values.
	t = FLA_task_alloc_init( func, cost, n_int_params, n_fla_params );
	//fprintf( stderr, "inited\n" ); fflush(stderr);
	
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
}


void FLA_queue_push( FLA_Task* t )
{
	int r_val;
	
	// Lock the mutex so we may modify the task queue.
	r_val = pthread_mutex_lock( &tq.lock );
	assert( r_val == 0 );
	
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
	tq.n_tasks++;
	
	// Unlock the mutex. We're done modifying the task queue.
	r_val = pthread_mutex_unlock( &tq.lock );
	assert( r_val == 0 );
}

void FLA_queue_flush()
{
	int       i;
	FLA_Task* t;
	FLA_Task* t_next_task;
	
	// Start with the head of the queue. If the queue is empty, then
	// the while loop will not execute.
	t = tq.head;
	
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
	tq.n_tasks = 0;
	tq.head    = NULL;
	tq.tail    = NULL;
}

FLA_Task* FLA_queue_pop()
{
	FLA_Task* t;
	int       r_val;
	
	// Lock the mutex so we may modify the task queue.
	r_val = pthread_mutex_lock( &tq.lock );
	assert( r_val == 0 );
	
	// Simply return NULL if queue is empty.
	if( tq.n_tasks == 0 )
	{
		t = NULL;
	}
	else
	{
		// Get the pointer to the head task on the queue.
		t = tq.head;
		
		// Change the head pointer to reflect popping the task.
		tq.head = tq.head->next_task;
		
		// Decrement the number of tasks.
		tq.n_tasks--;
		
		// If the queue is now empty, then nullify the tail pointer.
		if( tq.n_tasks == 0 )
		{
			tq.tail = NULL;
		}
		
		// Nullify the next_task field of the task we just removed from the queue.
		t->next_task = NULL;
	}
	
	// Unlock the mutex. We're done modifying the task queue.
	r_val = pthread_mutex_unlock( &tq.lock );
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

