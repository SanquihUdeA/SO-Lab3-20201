#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <pthread.h>
#include <sys/time.h>
#include <semaphore.h>

typedef struct Args{
	int curr_index;
	int end_index;
	int n_threads;
	int max_iters;
	int *n_iter;
	int it;
	int p;
	double scalar;
	double *X;
	double *Y;
	double *Y_avgs;
} Saxpy_args;

// Semphs. Def.
sem_t mutex;

// Sep_ args_ shared
Saxpy_args** create_args (int n_elements, int n_threads, int max_iters, int *n_iter, double scalar,  double *X, double *Y, double *Y_avgs) {
	Saxpy_args** parms = (Saxpy_args**) malloc(n_threads*sizeof(Saxpy_args*));

	int arr_chunk = n_elements / n_threads;	
	int n_leftovers = n_elements % n_threads;	
	int curr_index = 0;	
	int end_index = 0;	
	int n_threads_copy = n_threads;	
	int n; 	
	while (n_threads--) {	
		parms[n_threads] = (Saxpy_args*) malloc(sizeof(Saxpy_args));	
		if (n_leftovers) {	
			n = arr_chunk + 1;	
			n_leftovers--;	
		} else { n = arr_chunk; }
		end_index = curr_index + n - 1;	
		parms[n_threads] -> curr_index = curr_index;
		parms[n_threads] -> p = n_elements;
		parms[n_threads] -> end_index = end_index;	
		parms[n_threads] -> n_threads = n_threads_copy;	
		parms[n_threads] -> max_iters = max_iters;	
		parms[n_threads] -> n_iter = n_iter;	
		parms[n_threads] -> scalar = scalar;	
		parms[n_threads] -> X = X;	
		parms[n_threads] -> Y = Y;	
		parms[n_threads] -> Y_avgs = Y_avgs;	
		curr_index = end_index + 1;	
	}	
	return parms;
}

// Function to executed by threads...
void *SAXPYOp (void *args) {
	//aló
	Saxpy_args* parms = (Saxpy_args*) args;

	int p = parms->p;
	int it = parms->it;
	int end_index = parms->end_index;
	int max_iters = parms->max_iters;
	int curr_index = parms->curr_index;
	double scalar = parms->scalar;
	double* X = parms->X;
	double* Y = parms->Y;
	double* Y_avgs = parms->Y_avgs;

	int i;
	double result;

	#ifdef DEBUG
	  	printf("Thread values start = %d, end = %d, max_iters = %d, p = %d \n", ini, end, max_iters, p);
	#endif
	//SAXPY iterative SAXPY mfunction
	for (it = 0; it < max_iters; it++) {
		result = 0; //Variable local para optimizar el calculo del promedio
		for (i = curr_index; i < end_index; i++) {
			Y[i] = Y[i] + scalar * X[i];
			result += Y[i];
		}
		sem_wait(&mutex); 
		Y_avgs[it] += result / p; //Seccion critica protegida con semaforo binario
		sem_post(&mutex);
	}
  return NULL;
}

pthread_t *create_thread_arr (int n_threads, Saxpy_args **args) {
	pthread_t *threads = (pthread_t*) malloc(n_threads*sizeof(pthread_t));	
	while (n_threads--) {
		pthread_create(&threads[n_threads], NULL, SAXPYOp, args[n_threads]);
  }
	return threads;	
}

void wait_for_threads (pthread_t *threads, int n_threads) {	
	while (n_threads--) {	
		pthread_join(threads[n_threads], NULL);	
	}	
}

int main(int argc, char* argv[]) {
	// Variables to obtain command line parameters
	unsigned int seed = 1;
	int p = 10000000;
	int n_threads = 2;
	int max_iters = 1000;
	// Variables to perform SAXPY operation
	double* X;
	double a;
	double* Y;
	double* Y_avgs;
	int i, it;
	// Variables to get execution time
	struct timeval t_start, t_end;
	double exec_time;

	// Getting input values
	int opt;
	while ((opt = getopt(argc, argv, ":p:s:n:i:")) != -1) {  
		switch(opt){  
			case 'p':  
				printf("vector size: %s\n", optarg);
				p = strtol(optarg, NULL, 10);
				assert(p > 0 && p <= 2147483647);
				break;  
			case 's':  
				printf("seed: %s\n", optarg);
				seed = strtol(optarg, NULL, 10);
				break;
			case 'n':  
				printf("threads number: %s\n", optarg);
				n_threads = strtol(optarg, NULL, 10);
				break;  
			case 'i':  
				printf("max. iterations: %s\n", optarg);
				max_iters = strtol(optarg, NULL, 10);
				break;  
			case ':':  
				printf("option -%c needs a value\n", optopt);  
				break;  
			case '?':  
				fprintf(stderr, "Usage: %s [-p <vector size>] [-s <seed>] [-n <threads number>]\n", argv[0]);
				exit(EXIT_FAILURE);
		}  
	}  
	srand(seed);

	printf("p = %d, seed = %d, n_threads = %d, max_iters = %d\n", \
	p, seed, n_threads, max_iters);	

	// initializing data
	X = (double*) malloc(sizeof(double) * p);
	Y = (double*) malloc(sizeof(double) * p);
	Y_avgs = (double*) malloc(sizeof(double) * max_iters);

	for(i = 0; i < p; i++){
		X[i] = (double)rand() / RAND_MAX;
		Y[i] = (double)rand() / RAND_MAX;
	}
	for(i = 0; i < max_iters; i++){
		Y_avgs[i] = 0.0;
	}
	a = (double)rand() / RAND_MAX;

	#ifdef DEBUG
		printf("vector X= [ ");
		for(i = 0; i < p-1; i++){
			printf("%f, ",X[i]);
		}
		printf("%f ]\n",X[p-1]);

		printf("vector Y= [ ");
		for(i = 0; i < p-1; i++){
			printf("%f, ", Y[i]);
		}
		printf("%f ]\n", Y[p-1]);

		printf("a= %f \n", a);	
	#endif

	/**
	*
	* Function to parallelize 
	*
	**/
	gettimeofday(&t_start, NULL);
  	sem_init(&mutex,0,1); //Inicializacion de semaforo binario
	//SAXPY iterative SAXPY mfunction

  	it = 0;
	Saxpy_args **args = create_args(p, n_threads, max_iters, &it, a,  X, Y, Y_avgs);
	pthread_t *threads;
	threads = create_thread_arr(n_threads, args); // n_threads are created
	wait_for_threads(threads, n_threads);

	gettimeofday(&t_end, NULL);

	#ifdef DEBUG
		printf("RES: final vector Y= [ ");
		for(i = 0; i < p-1; i++){
			printf("%f, ", Y[i]);
		}
		printf("%f ]\n", Y[p-1]);
	#endif
	
	// Computing execution time
	exec_time = (t_end.tv_sec - t_start.tv_sec) * 1000.0;  // sec to ms
	exec_time += (t_end.tv_usec - t_start.tv_usec) / 1000.0; // us to ms
	printf("Execution time: %f ms \n", exec_time);
	printf("Last 3 values of Y: %f, %f, %f \n", Y[p-3], Y[p-2], Y[p-1]);
	printf("Last 3 values of Y_avgs: %f, %f, %f \n", Y_avgs[max_iters-3], Y_avgs[max_iters-2], Y_avgs[max_iters-1]);
	free(threads);
	return 0;
}