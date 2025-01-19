#ifndef EPO_TSP_H
#define EPO_TSP_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <mpi.h>

#define NUM_CITIES 5000
#define NUM_PENGUINS 50
#define MAX_ITERATIONS 300

/*
    Structure to represent a city in x and y coordinates
*/
typedef struct {
    double x, y;
} City;


/*
    Structure to represent a penguin with its path (path) and total cost (cost)
*/
typedef struct {
    int path[NUM_CITIES];
    double cost;
} Penguin;

double distance(City a, City b);
double calculate_path_length(City cities[], Penguin* penguin, int num_cities);
void generate_random_path(Penguin* penguin, int num_cities);
int contains(int* path, int city, int num_cities);
void enhanced_crossover(Penguin* follower, Penguin* leader, int num_cities);
void mutate(Penguin* penguin, int num_cities);
void emperor_penguin_optimization(City cities[], int num_cities, int num_penguins, int max_iterations, int rank, int size);
void read_cities_from_file(City cities[], const char* filename, int num_cities);

#endif 
