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


// Mutex
extern pthread_mutex_t mutex_sol;
extern pthread_cond_t cond_sol;

extern pthread_mutex_t mutex_sol_len;
extern pthread_cond_t cond_sol_len;

extern pthread_mutex_t mutex_cuts;
extern pthread_cond_t cond_cuts;

extern pthread_mutex_t mutex_jobs;
extern pthread_cond_t cond_jobs;

extern pthread_mutex_t mutex_print;	// Pour des raisons esthétiques
extern pthread_cond_t cond_print;

// Sémaphore
sem_t sem_threads_number;


struct threads_args {
	int * number_of_threads;
	int * jumps;
	int * length;
	long long int * cuttings;
	tsp_path_t * t_path;
	tsp_path_t * t_solution;
	int * solution_length;
	struct tsp_queue * jobs_list;
};

void * status;
void * PTHREAD_TERMINATED = (void *)123456789L;


static void generate_tsp_jobs (struct tsp_queue *q, int hops, int len, tsp_path_t path, long long int *cuts, tsp_path_t sol, int *sol_len, int depth)
{
    if (len >= minimum) {
        (*cuts)++ ;
        return;
    }

    pthread_cond_init(&cond_jobs, NULL);

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
	get_job(conversion->jobs_list, *(conversion->t_path), conversion->jumps, conversion->length);
// 	printf("%lX appelle tsp()\n", pthread_self());
	tsp(*(conversion->jumps), *(conversion->length), *(conversion->t_path), conversion->cuttings, *(conversion->t_solution), conversion->solution_length);
	(*(conversion->number_of_threads))--;
// 	printf("Je suis le thread %lX : j'ai terminé mon job\n", pthread_self());
// 	pthread_cond_signal(&cond_threads_number);
	sem_post(&sem_threads_number);

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

	// Tableau de TID
	pthread_t *threads_table = NULL;
	struct threads_args * arguments_for_threads = NULL;
	int i = 0;

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

	/* Génération des threads */
	threads_table = (pthread_t *) calloc(nb_threads, sizeof(pthread_t));
	arguments_for_threads = (struct threads_args *) calloc(nb_threads, sizeof(struct threads_args));

	sem_init(&sem_threads_number, 0, nb_threads);
	pthread_cond_init(&cond_sol, NULL);
	pthread_cond_init(&cond_sol_len, NULL);
	pthread_cond_init(&cond_cuts, NULL);
	pthread_cond_init(&cond_print, NULL);
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

        int hops = 0, len = 0;

		(arguments_for_threads+i)->number_of_threads = &i;
		(arguments_for_threads+i)->jumps = &hops;
		(arguments_for_threads+i)->length = &len;
		(arguments_for_threads+i)->cuttings = &cuts;
		(arguments_for_threads+i)->t_path = &solution;
		(arguments_for_threads+i)->t_solution = &sol;
		(arguments_for_threads+i)->solution_length = &sol_len;
		(arguments_for_threads+i)->jobs_list = &q;
		pthread_create(threads_table+i, NULL, threads_loop, (void *) (arguments_for_threads+i));
		i++;
// 		printf("%d threads créés\n", i);
    }

	int j = 0;
    for (j = 0; j < nb_threads; j++) {
		pthread_join(threads_table[i], &status);
	}

    sem_destroy(&sem_threads_number);
	pthread_mutex_destroy(&mutex_sol);
	pthread_mutex_destroy(&mutex_sol_len);
	pthread_mutex_destroy(&mutex_jobs);
	pthread_mutex_destroy(&mutex_cuts);
	pthread_mutex_destroy(&mutex_print);
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
