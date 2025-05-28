#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>  

#define MAX_CITIES 16000 // Maximum number of cities supported
#define MAX_PENGUINS 500 // Maximum number of cities supported

typedef struct {
    int id; // Unique city identifier
    double x, y; // Coordinates of city on 2D plane
} City;

/**
 * Compute Euclidean distance between two cities.
 *
 * @param a First city
 * @param b Second city
 * @return Straight-line distance between cities a and b
 */
double calculate_distance(City a, City b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return sqrt(dx * dx + dy * dy);
}

/**
 * Calculate total length of a closed tour visiting all cities once.
 *
 * @param path Array of cities in visit order
 * @param num_cities Number of cities in the path
 * @return Total tour length, including return to start
 */
double calculate_path_length(City path[], int num_cities) {
    double total = 0.0;
    for (int i = 0; i < num_cities - 1; i++) {
        total += calculate_distance(path[i], path[i + 1]);
    }
    // Close the loop: last city back to first
    total += calculate_distance(path[num_cities - 1], path[0]);
    return total;
}

/**
 * Read city coordinates from a file, one pair per line.
 *
 * @param filename Path to input file: each line "x y"
 * @param cities Preallocated array to fill with city data
 * @param num_cities Number of cities expected
 * @return 1 on success, 0 on error
 */
int read_cities(const char *filename, City cities[], int num_cities) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Error: Unable to open file '%s'\n", filename);
        return 0;
    }

    for (int i = 0; i < num_cities; i++) {
        if (fscanf(file, "%lf %lf", &cities[i].x, &cities[i].y) != 2) {
            printf("Error: City %d not read correctly\n", i);
            fclose(file);
            return 0;
        }
    }
    fclose(file);
    printf(">> Read %d cities from '%s'\n", num_cities, filename);
    return 1;
}

/**
 * Randomly shuffle the order of an array of cities.
 *
 * @param path Array of cities to shuffle
 * @param num_cities Number of cities in the array
 */
void shuffle_path(City path[], int num_cities) {
    // Perform 2*num_cities random swaps for thorough mixing
    for (int i = 0; i < num_cities * 2; i++) {
        int a = rand() % num_cities;
        int b = rand() % num_cities;
        City temp = path[a];
        path[a] = path[b];
        path[b] = temp;
    }
}

/**
 * Copy a contiguous 80% segment of the leader's path into a follower,
 * then fill remaining positions from unused cities in random order.
 *
 * @param leader Source path of the best penguin
 * @param follower Destination path to be updated
 * @param num_cities Total number of cities
 * @param start_percent Fractional start index (0.0 to 0.2)
 */
void copy_leader_segment(City leader[], City follower[], int num_cities, double start_percent) {
    int segment_size = num_cities * 0.8;
    int start = (int)(num_cities * start_percent);
    int end = start + segment_size;
    if (end > num_cities) end = num_cities;

    int used[MAX_CITIES] = {0}; // Flags for cities already copied

    // Copy the central segment from leader to follower
    for (int i = start; i < end; i++) {
        follower[i] = leader[i];
        // Mark that city as used by matching coordinates
        for (int j = 0; j < num_cities; j++) {
            if (leader[i].x == leader[j].x && leader[i].y == leader[j].y) {
                used[j] = 1;
                break;
            }
        }
    }

    // Fill remaining positions with unused cities
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
 * Mutate a path by swapping city positions based on a cooling schedule.
 * Number of swaps decreases as iterations progress.
 *
 * @param path Array of cities to mutate
 * @param num_cities Number of cities
 * @param iter Current iteration index
 * @param max_iter Total iterations for cooling schedule
 */
void mutate_path(City path[], int num_cities, int iter, int max_iter) {
    double cooling_factor = 1.0 - ((double)iter / max_iter);
    int swaps = (int)(num_cities * 0.2 * cooling_factor);
    if (swaps < 2) swaps = 2;

    for (int i = 0; i < swaps; i++) {
        int a = rand() % num_cities;
        int b = rand() % num_cities;
        City temp = path[a];
        path[a] = path[b];
        path[b] = temp;
    }
}

/**
 * Core Evolutionary Penguin Optimization (EPO) loop with OpenMP parallelism.
 * Each thread evaluates a subset of penguin paths, finds the best local path,
 * then updates all penguins based on the global leader.
 *
 * @param cities Array of all city coordinates
 * @param num_cities Number of cities in the problem
 * @param num_penguins Total penguin agents
 * @param num_iterations EPO iterations to perform
 * @param best_path Output array to store the best discovered path
 */
void run_epo(City cities[], int num_cities, int num_penguins, int num_iterations, City best_path[]) {
    printf(">> Start EPO with %d penguins, %d cities, %d iterations\n", num_penguins, num_cities, num_iterations);

    // Open CSV log file to track iteration progress
    FILE *log_file = fopen("epo_log.csv", "w");
    if (!log_file) {
        printf("Error: Unable to create log file.\n");
        return;
    }
    fprintf(log_file, "iteration,penguin,length\n");
 
    static City population[MAX_PENGUINS][MAX_CITIES]; // Local population of paths

    // Initialize random paths in parallel
    #pragma omp parallel for
    for (int i = 0; i < num_penguins; i++) {
        for (int j = 0; j < num_cities; j++) {
            population[i][j] = cities[j];
        }
        shuffle_path(population[i], num_cities);
    }
    printf(">> Generate random initial paths for penguins\n");

    double best_length = INFINITY;
    int last_improve = 0;

    // Main optimization loop
    for (int iter = 0; iter < num_iterations; iter++) {

        // Stop criterion
        if (iter - last_improve >= 100) {
            printf(">> Stop: no improvement after %d iterations, exiting at step %d\n",
                iter - last_improve, iter);
            break;
        }

        printf(">> Iteration %d/%d\n", iter + 1, num_iterations);
        int leader_index = 0;
        double leader_length = INFINITY;

        // Parallel evaluation of each penguin's path
        #pragma omp parallel
        {
            int local_leader = -1;
            double local_best = INFINITY;

            // Distribute work among threads
            #pragma omp for private(current_length) nowait
            for (int i = 0; i < num_penguins; i++) {
                double current_length = calculate_path_length(population[i], num_cities);

                // Thread-safe logging
                #pragma omp critical
                fprintf(log_file, "%d,%d,%.2f\n", iter + 1, i, current_length);

                // Track best in this thread
                if (current_length < local_best) {
                    local_best = current_length;
                    local_leader = i;
                }
            }

            // Update global leader under critical section
            #pragma omp critical
            {
                if (local_best < leader_length) {
                    leader_length = local_best;
                    leader_index = local_leader;
                }
            }
        }

        printf("   >> Leader: Penguin %d with length %.2f\n", leader_index, leader_length);

        // Update overall best path if improved
        if (leader_length < best_length) {
            best_length = leader_length;
            last_improve = iter;
            for (int i = 0; i < num_cities; i++) {
                best_path[i] = population[leader_index][i];
            }
        }

        // Copy-and-mutate each penguin based on leader
        #pragma omp parallel for
        for (int i = 0; i < num_penguins; i++) {
            if (i == leader_index) continue;
            double start_percent = (i % 5) * 0.2;
            copy_leader_segment(population[leader_index], population[i], num_cities, start_percent);
            mutate_path(population[i], num_cities, iter, num_iterations);
        }
    }

    fclose(log_file);
    printf(">> End of EPO\n");
}

/**
 * Program entry point: sets up environment, reads input, and executes EPO.
 *
 * Usage: %s <city_file> <num_cities> <num_iterations>
 */
int main(int argc, char *argv[]) {
    // Display number of OpenMP threads available
    printf(">> OpenMP: %d threads available\n", omp_get_max_threads());

    if (argc != 4) {
        printf("Usage: %s <city_file> <city_number> <iteration_number>\n", argv[0]);
        return 1;
    }

    // Seed random number generator once
    srand(time(NULL));

    const char *filename = argv[1];
    int num_cities = atoi(argv[2]);
    int num_iterations = atoi(argv[3]);

    // Determine number of penguins (agents)
    int num_penguins = num_cities / 100;
    if (num_penguins < 5) num_penguins = 5;

    if (num_cities > MAX_CITIES) {
        printf("Error: City number too large. Max allowed: %d\n", MAX_CITIES);
        return 1;
    }
    if (num_penguins > MAX_PENGUINS) {
        printf("Error: Number of penguins too large. Max allowed: %d\n", MAX_PENGUINS);
        return 1;
    }

    // Read city coordinates from file
    City cities[MAX_CITIES];
    if (!read_cities(filename, cities, num_cities)) {
        printf("Error reading cities.\n");
        return 1;
    }

    // Array to store the best path found
    City best_path[MAX_CITIES];

    // Start timer
    double start_time = omp_get_wtime();

    // Run the EPO algorithm
    run_epo(cities, num_cities, num_penguins, num_iterations, best_path);

    // End timer
    double end_time = omp_get_wtime();

    printf("\n>> Best route calculation completed\n");
    double best_length = calculate_path_length(best_path, num_cities);
    printf("\nBest path found (length: %.2f)\n", best_length);
    printf("\n>> Execution time: %.3f seconds\n", end_time - start_time);

    return 0;
}
