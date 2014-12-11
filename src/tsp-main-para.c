#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <assert.h>
#include <complex.h>
#include <stdbool.h>
#include <unistd.h>
#include <pthread.h>
#include <semaphore.h>

#include "tsp-types.h"
#include "tsp-job.h"
#include "tsp-genmap.h"
#include "tsp-print.h"
#include "tsp-tsp.h"


/* macro de mesure de temps, retourne une valeur en nanosecondes */
#define TIME_DIFF(t1, t2) \
  ((t2.tv_sec - t1.tv_sec) * 1000000000ll + (long long int) (t2.tv_nsec - t1.tv_nsec))


/* tableau des distances */
tsp_distance_matrix_t distance ={};

/** Paramètres **/

/* nombre de villes */
int nb_towns=10;
/* graine */
long int myseed= 0;
/* nombre de threads */
int nb_threads=1;

/* affichage SVG */
bool affiche_sol= false;


// Mutex (tous utilises dans des fonctions d'autres modules)
extern pthread_mutex_t mutex_sol;		// Declare dans tsp-tsp.c
extern pthread_mutex_t mutex_sol_len;	// Declare dans tsp-tsp.c
extern pthread_mutex_t mutex_cuts;		// Declare dans tsp-tsp.c
extern pthread_mutex_t mutex_jobs;		// Declare dans tsp-job.c

extern pthread_mutex_t mutex_print;		// Pour des raisons esthetiques
										// Declare dans tsp-print.c

// Semaphore
sem_t sem_threads_number;

// Structures pour les threads
struct threads_args {
	int * threads_finished_table;
	int thread_index;
	int jumps;
	int length;
	long long int * cuttings;
	tsp_path_t * t_path;
	tsp_path_t * t_solution;
	int * solution_length;
	struct tsp_queue * jobs_list;
};

void * status;
void * PTHREAD_TERMINATED = (void *)123456789L;

/** Fin des paramètres **/


static void generate_tsp_jobs (struct tsp_queue *q, int hops, int len, tsp_path_t path, long long int *cuts, tsp_path_t sol, int *sol_len, int depth)
{
    if (len >= minimum) {
        (*cuts)++ ;
        return;
    }

    if (hops == depth) {
        /* On enregistre du travail à faire plus tard... */
        add_job (q, path, hops, len);
    } else {
        int me = path [hops - 1];
        for (int i = 0; i < nb_towns; i++) {
            if (!present (i, hops, path)) {
                path[hops] = i;
                int dist = distance[me][i];
                generate_tsp_jobs (q, hops + 1, len + dist, path, cuts, sol, sol_len, depth);
            }
        }
    }
}

static void usage(const char *name) {
  fprintf (stderr, "Usage: %s [-s] <ncities> <seed> <nthreads>\n", name);
  exit (3);
}


static void *threads_loop(void * arg) {
	struct threads_args * conversion = (struct threads_args *) arg;

// 	printf("%lX appelle get_job()\n", pthread_self());
	get_job(conversion->jobs_list, *(conversion->t_path), &(conversion->jumps), &(conversion->length));
// 	printf("%lX appelle tsp()\n", pthread_self());
	tsp(conversion->jumps, conversion->length, *(conversion->t_path), conversion->cuttings, *(conversion->t_solution), conversion->solution_length);

	// Le thread signale qu'il a fini (test dans le main)
	*((conversion->threads_finished_table)+conversion->thread_index) = 1;

	// Le thread met a jour le semaphore
	sem_post(&sem_threads_number);

	// printf("Je suis le thread %lX : j'ai termine mon job\n", pthread_self());
	return PTHREAD_TERMINATED;
}


int main (int argc, char **argv)
{
    unsigned long long perf;
    tsp_path_t path;
    tsp_path_t sol;	// Accès en exclusion mutuelle
    int sol_len;	// Accès en exclusion mutuelle
    long long int cuts = 0;	// Accès en exclusion mutuelle (moins critique)
    struct tsp_queue q;
	struct timespec t1, t2;

	/* Tableau de threads, et tableau de signalisation des threads finis */
	pthread_t *threads_table = NULL;
	int *threads_finished = NULL;
	// et tableau de structures d'arguments pour les threads
	struct threads_args * arguments_for_threads = NULL;
	// Indice pour la creation de threads
	int next_thread_index = 0;
	int next_thread_index_buffer = 0;
	/*  */

    /* lire les arguments */
    int opt;
    while ((opt = getopt(argc, argv, "s")) != -1) {
      switch (opt) {
      case 's':
	affiche_sol = true;
	break;
      default:
	usage(argv[0]);
	break;
      }
    }

    if (optind != argc-3)
      usage(argv[0]);

    nb_towns = atoi(argv[optind]);
    myseed = atol(argv[optind+1]);
    nb_threads = atoi(argv[optind+2]);
    assert(nb_towns > 0);
    assert(nb_threads > 0);

	/* Generation des threads */
	threads_table = (pthread_t *) calloc(nb_threads, sizeof(pthread_t));
	threads_finished = (int *) calloc(nb_threads, sizeof(int));
	arguments_for_threads = (struct threads_args *) calloc(nb_threads, sizeof(struct threads_args));

	// Initialisation du semaphore
	sem_init(&sem_threads_number, 0, nb_threads);
	/*  */

    minimum = INT_MAX;

    /* generer la carte et la matrice de distance */
    fprintf (stderr, "ncities = %3d\n", nb_towns);
    genmap ();

    init_queue (&q);

    clock_gettime (CLOCK_REALTIME, &t1);

    memset (path, -1, MAX_TOWNS * sizeof (int));
    path[0] = 0;

    /* mettre les travaux dans la file d'attente */
    generate_tsp_jobs (&q, 1, 0, path, &cuts, sol, & sol_len, 3);
    no_more_jobs (&q);

    /* calculer chacun des travaux */
    tsp_path_t solution;
    memset (solution, -1, MAX_TOWNS * sizeof (int));
    solution[0] = 0;
	while (!empty_queue (&q)) {
		sem_wait(&sem_threads_number);

		// On gère le fait qu'il y ait plus de jobs que de threads
		if (next_thread_index >= nb_threads) {
			next_thread_index_buffer = next_thread_index;
			next_thread_index = 0;
			while ((next_thread_index < nb_threads) && (threads_finished[next_thread_index] == 0)) {
				next_thread_index++;
			}
			if (threads_finished[next_thread_index] == 1) {
				pthread_join(threads_table[next_thread_index], &status);	// On attend un seul thread a la fois et seulement si besoin est
			}
			threads_finished[next_thread_index] = 0;
		}

        int hops = 0, len = 0;

		// Initialisation de la structure d'arguments
		(arguments_for_threads+next_thread_index)->threads_finished_table = threads_finished;
		(arguments_for_threads+next_thread_index)->thread_index = next_thread_index;
		(arguments_for_threads+next_thread_index)->jumps = hops;
		(arguments_for_threads+next_thread_index)->length = len;
		(arguments_for_threads+next_thread_index)->cuttings = &cuts;
		(arguments_for_threads+next_thread_index)->t_path = &solution;
		(arguments_for_threads+next_thread_index)->t_solution = &sol;
		(arguments_for_threads+next_thread_index)->solution_length = &sol_len;
		(arguments_for_threads+next_thread_index)->jobs_list = &q;

		// Creation d'un thread
		pthread_create(threads_table + next_thread_index, NULL, threads_loop, (void *) (arguments_for_threads + next_thread_index));

		// Correctif pour l'indice du thread a creer
		next_thread_index = next_thread_index_buffer;
    }

    fprintf(stderr, "On attend les threads.\n");
	int j = 0;
    for (j = 0; j < nb_threads; j++) {
		pthread_join(threads_table[j], &status);
	}

	fprintf(stderr, "On libere les structures d'exclusion mutuelle.\n");
    sem_destroy(&sem_threads_number);
	pthread_mutex_destroy(&mutex_sol);
	pthread_mutex_destroy(&mutex_sol_len);
	pthread_mutex_destroy(&mutex_cuts);
	pthread_mutex_destroy(&mutex_jobs);
	pthread_mutex_destroy(&mutex_print);

	fprintf(stderr, "On libere les structures liees aux threads.\n");
    free(arguments_for_threads);
	free(threads_table);

    clock_gettime (CLOCK_REALTIME, &t2);

    if (affiche_sol)
      print_solution_svg (sol, sol_len);

    perf = TIME_DIFF (t1,t2);
    printf("<!-- # = %d seed = %ld len = %d threads = %d time = %lld.%03lld ms ( %lld coupures ) -->\n",
	   nb_towns, myseed, sol_len, nb_threads,
	   perf/1000000ll, perf%1000000ll, cuts);

    return 0 ;
}
