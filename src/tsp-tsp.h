#ifndef TSP_TSP_H
#define TSP_TSP_H

/* dernier minimum trouv� */
extern int minimum;

struct tsp_threads {
	int jumps;
	int length;
	tsp_path_t * path;
	long long int * cuttings;
	tsp_path_t * solution;
	int * solution_length;
};

int present (int city, int hops, tsp_path_t path);
void tsp (int hops, int len, tsp_path_t path, long long int *cuts, tsp_path_t sol, int *sol_len);

#endif
