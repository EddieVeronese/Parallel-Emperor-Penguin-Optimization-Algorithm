#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAX_CITIES 16000    // Maximum number of cities in the problem
#define MAX_PENGUINS 500    // Maximum number of penguin agents

typedef struct {
    int id;      // Unique identifier for the city (unused)
    double x, y; // Coordinates of the city on a 2D plane
} City;

/**
 * Calculate the Euclidean distance between two cities.
 */
double calculate_distance(City a, City b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return sqrt(dx*dx + dy*dy);
}

/**
 * Compute the total length of a closed path that visits all cities exactly once and returns to start.
 */
double calculate_path_length(City path[], int num_cities) {
    double total = 0.0;
    for (int i = 0; i < num_cities - 1; i++) {
        total += calculate_distance(path[i], path[i+1]);
    }
    total += calculate_distance(path[num_cities-1], path[0]);
    return total;
}

/**
 * Load city coordinates from a file.
 */
int read_cities(const char *filename, City cities[], int num_cities) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Unable to open file '%s'\n", filename);
        return 0;
    }
    for (int i = 0; i < num_cities; i++) {
        if (fscanf(file, "%lf %lf", &cities[i].x, &cities[i].y) != 2) {
            fprintf(stderr, "Error: City %d not read correctly\n", i);
            fclose(file);
            return 0;
        }
    }
    fclose(file);
    return 1;
}

/**
 * Randomly permute the order of cities in a path.
 */
void shuffle_path(City path[], int num_cities) {
    for (int i = 0; i < num_cities * 2; i++) {
        int a = rand() % num_cities;
        int b = rand() % num_cities;
        City tmp = path[a];
        path[a] = path[b];
        path[b] = tmp;
    }
}

/**
 * Copy an 80% segment from the leader path into a follower,
 * then fill remaining positions with the unused cities.
 */
void copy_leader_segment(City leader[], City follower[], int num_cities, double start_frac) {
    int seg_size = (int)(num_cities * 0.8);
    int start = (int)(num_cities * start_frac);
    int end   = start + seg_size;
    if (end > num_cities) end = num_cities;

    int used[MAX_CITIES] = {0};

    // copy the segment
    for (int i = start; i < end; i++) {
        follower[i] = leader[i];
        // mark as used by matching coords
        for (int j = 0; j < num_cities; j++) {
            if (leader[i].x == leader[j].x && leader[i].y == leader[j].y) {
                used[j] = 1;
                break;
            }
        }
    }

    // fill the rest
    for (int i = 0; i < num_cities; i++) {
        if (i >= start && i < end) continue;
        while (1) {
            int r = rand() % num_cities;
            if (!used[r]) {
                follower[i] = leader[r];
                used[r] = 1;
                break;
            }
        }
    }
}

/**
 * Mutate a path by swapping positions according to a cooling schedule.
 */
void mutate_path(City path[], int num_cities, int iter, int max_iter) {
    double cool = 1.0 - ((double)iter / max_iter);
    int swaps = (int)(num_cities * 0.2 * cool);
    if (swaps < 2) swaps = 2;
    for (int i = 0; i < swaps; i++) {
        int a = rand() % num_cities;
        int b = rand() % num_cities;
        City tmp = path[a];
        path[a] = path[b];
        path[b] = tmp;
    }
}

/**
 * Serial Evolutionary Penguin Optimization (EPO).
 */
void run_epo(City cities[], int num_cities, int num_penguins,
             int num_iterations, City best_path[]) {
    printf(">> Start EPO with %d penguins, %d cities, %d iterations\n",
           num_penguins, num_cities, num_iterations);

    // Log file for analysis
    FILE *logf = fopen("epo_log.csv", "w");
    if (!logf) {
        fprintf(stderr, "Error: Unable to create log file\n");
        return;
    }
    fprintf(logf, "iteration,penguin,length\n");

    static City population[MAX_PENGUINS][MAX_CITIES];

    // Initialize random paths
    for (int i = 0; i < num_penguins; i++) {
        for (int j = 0; j < num_cities; j++) {
            population[i][j] = cities[j];
        }
        shuffle_path(population[i], num_cities);
    }
    printf(">> Generated initial random paths for penguins\n");

    double global_best = INFINITY;

    // Main loop
    for (int iter = 0; iter < num_iterations; iter++) {
        printf(">> Iteration %d/%d\n", iter+1, num_iterations);

        // find best local
        int leader_idx = 0;
        double leader_len = calculate_path_length(population[0], num_cities);

        for (int i = 0; i < num_penguins; i++) {
            double len = calculate_path_length(population[i], num_cities);
            fprintf(logf, "%d,%d,%.2f\n", iter+1, i, len);
            if (len < leader_len) {
                leader_len = len;
                leader_idx = i;
            }
        }

        printf("   >> Leader: Penguin %d with length %.2f\n", leader_idx, leader_len);

        // update global best
        if (leader_len < global_best) {
            global_best = leader_len;
            for (int j = 0; j < num_cities; j++) {
                best_path[j] = population[leader_idx][j];
            }
        }

        // evolve population
        for (int i = 0; i < num_penguins; i++) {
            if (i == leader_idx) continue;
            double start_frac = (i % 5) * 0.2;
            copy_leader_segment(population[leader_idx],
                                population[i],
                                num_cities,
                                start_frac);
            mutate_path(population[i],
                        num_cities,
                        iter,
                        num_iterations);
        }
    }

    fclose(logf);
    printf(">> End of EPO\n");
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <city_file> <num_cities> <num_iterations>\n", argv[0]);
        return EXIT_FAILURE;
    }

    srand((unsigned)time(NULL));

    const char *filename    = argv[1];
    int num_cities          = atoi(argv[2]);
    int num_iterations      = atoi(argv[3]);
    int num_penguins        = num_cities / 100;
    if (num_penguins < 5)   num_penguins = 5;
    if (num_penguins > MAX_PENGUINS) {
        fprintf(stderr, "Error: too many penguins (max %d)\n", MAX_PENGUINS);
        return EXIT_FAILURE;
    }
    if (num_cities > MAX_CITIES) {
        fprintf(stderr, "Error: too many cities (max %d)\n", MAX_CITIES);
        return EXIT_FAILURE;
    }

    City *cities = malloc(sizeof(City) * num_cities);
    if (!read_cities(filename, cities, num_cities)) {
        free(cities);
        return EXIT_FAILURE;
    }

    City *best_path = malloc(sizeof(City) * num_cities);

    double t0 = (double)clock() / CLOCKS_PER_SEC;
    run_epo(cities, num_cities, num_penguins, num_iterations, best_path);
    double t1 = (double)clock() / CLOCKS_PER_SEC;

    double best_len = calculate_path_length(best_path, num_cities);
    printf("\n>> Best path found (length: %.2f)\n", best_len);
    printf(">> Execution time: %.3f seconds\n", t1 - t0);

    free(cities);
    free(best_path);
    return EXIT_SUCCESS;
}
