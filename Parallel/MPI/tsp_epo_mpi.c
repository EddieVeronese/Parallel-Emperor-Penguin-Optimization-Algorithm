#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define MAX_CITIES 16000 // Maximum number of cities in the problem
#define MAX_PENGUINS 500  // Maximum number of penguin agents

typedef struct {
    int id; // Unique identifier for the city
    double x, y;  // Coordinates of the city on a 2D plane
} City;

/**
 * Calculate the Euclidean distance between two cities.
 *
 * @param a First city
 * @param b Second city
 * @return The straight-line distance between city a and city b
 */
double calculate_distance(City a, City b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return sqrt(dx * dx + dy * dy);
}

/**
 * Compute the total length of a closed path that visits all cities exactly once and returns to start.
 *
 * @param path Array of cities in visit order
 * @param num_cities Number of cities in the path
 * @return Total circumference of the tour
 */
double calculate_path_length(City path[], int num_cities) {
    double total = 0.0;
    for (int i = 0; i < num_cities - 1; i++) {
        total += calculate_distance(path[i], path[i + 1]);
    }
    total += calculate_distance(path[num_cities - 1], path[0]);
    return total;
}

/**
 * Load city coordinates from a file.
 *
 * @param filename Path to text file containing one city per line: "x y"
 * @param cities Preallocated array to store city data
 * @param num_cities Expected number of cities to read
 * @return True if successful, false on I/O error
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
    return 1;
}

/**
 * Randomly permute the order of cities in a path.
 *
 * @param path Array of cities to shuffle
 * @param num_cities Number of cities in the path
 */
void shuffle_path(City path[], int num_cities) {
    // Perform 2*num_cities random swaps for sufficient mixing
    for (int i = 0; i < num_cities * 2; i++) {
        int a = rand() % num_cities;
        int b = rand() % num_cities;
        City temp = path[a];
        path[a] = path[b];
        path[b] = temp;
    }
}

/**
 * Copy a central segment (80%) from the leader path into a follower path,
 * then fill the remaining positions with the leftover cities in random order.
 *
 * @param leader Best-known path (global leader)
 * @param follower Path to receive the segment
 * @param num_cities Number of cities in each path
 * @param start_percent Fractional start position (0.0 to 0.2) indicating where to copy
 */
void copy_leader_segment(City leader[], City follower[], int num_cities, double start_percent) {
    int segment_size = num_cities * 0.8;
    int start = (int)(num_cities * start_percent);
    int end = start + segment_size;
    if (end > num_cities) end = num_cities;

    int used[MAX_CITIES] = {0}; // Tracks which cities have been placed

    // Copy the segment from the leader into the same positions in the follower
    for (int i = start; i < end; i++) {
        follower[i] = leader[i];
        // Mark the city as used by matching coordinates
        for (int j = 0; j < num_cities; j++) {
            if (leader[i].x == leader[j].x && leader[i].y == leader[j].y) {
                used[j] = 1;
                break;
            }
        }
    }

    // Fill all other positions with the unused cities
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
 * Apply a mutation (swapping) to a path based on a cooling schedule.
 * The number of swaps decreases as iterations increase.
 *
 * @param path Path to mutate
 * @param num_cities Number of cities in the path
 * @param iter Current iteration index (0-based)
 * @param max_iter Total number of iterations for cooling schedule
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
 * Core EPO (Evolutionary Penguin Optimization) algorithm driver.
 * Each MPI process handles a subset of "penguin" paths,
 * discovers its local best, then participates in a global reduction to find the global leader.
 * Followers then copy and mutate based on that leader.
 */
void run_epo(City cities[], int num_cities, int num_penguins, int num_iterations, City best_path[], int rank, int size) {
    if (rank == 0)
        printf(">> Start EPO with %d penguins, %d cities, %d iterations\n", num_penguins, num_cities, num_iterations);

    FILE *log_file = NULL;
    if (rank == 0) {
        // Root process writes progress logs to CSV
        log_file = fopen("epo_log.csv", "w");
        if (!log_file) {
            printf("Error: Unable to create log file.\n");
            return;
        }
        fprintf(log_file, "iteration,penguin,length\n");
    }

    // Distribute penguins evenly across MPI ranks
    int local_penguins = num_penguins / size;
    if (rank < num_penguins % size) local_penguins++;

    // Allocate local population of paths
    static City population[MAX_PENGUINS][MAX_CITIES];

    // Initialize random paths for each local penguin
    for (int i = 0; i < local_penguins; i++) {
        for (int j = 0; j < num_cities; j++) {
            population[i][j] = cities[j];
        }
        shuffle_path(population[i], num_cities);
    }

    if (rank == 0)
        printf(">> Generate random initial paths for penguins\n");

    double best_length = INFINITY;
    City leader_global[MAX_CITIES]; // Holds the current global best path

    // Main optimization iterations
    for (int iter = 0; iter < num_iterations; iter++) {
        if (rank == 0)
            printf(">> Iteration %d/%d\n", iter + 1, num_iterations);

        // Find the best local penguin and its length
        int local_leader_index = 0;
        double local_best = calculate_path_length(population[0], num_cities);

        for (int i = 0; i < local_penguins; i++) {
            double len = calculate_path_length(population[i], num_cities);
            if (rank == 0 && log_file)
                fprintf(log_file, "%d,%d,%.2f\n", iter + 1, i + rank * local_penguins, len);
            if (len < local_best) {
                local_best = len;
                local_leader_index = i;
            }
        }

        // Copy local leader path
        City local_leader[MAX_CITIES];
        for (int i = 0; i < num_cities; i++) {
            local_leader[i] = population[local_leader_index][i];
        }

        // Structure to perform min reduction on path lengths with MPI
        struct {
            double len;
            int rank;
        } in = { local_best, rank }, out;

        // Global reduction: find smallest length and its rank
        MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

        // The rank owning the best path becomes the global leader
        if (rank == out.rank) {
            for (int i = 0; i < num_cities; i++)
                leader_global[i] = local_leader[i];
        }

        // Broadcast the global leader path to all processes
        MPI_Bcast(leader_global, sizeof(City) * num_cities, MPI_BYTE, out.rank, MPI_COMM_WORLD);

        if (rank == 0)
            printf("   >> Leader: Penguin %d (from process %d) with length %.2f\n",
                   local_leader_index + out.rank * local_penguins, out.rank, out.len);

        // Update the overall best if improved
        if (out.len < best_length) {
            best_length = out.len;
            for (int i = 0; i < num_cities; i++)
                best_path[i] = leader_global[i];
        }

        // Each penguin (except the leader) copies and mutates based on global leader
        for (int i = 0; i < local_penguins; i++) {
            if (i == local_leader_index && rank == out.rank) continue;
            double start_percent = (i % 5) * 0.2;
            copy_leader_segment(leader_global, population[i], num_cities, start_percent);
            mutate_path(population[i], num_cities, iter, num_iterations);
        }
    }

    // Clean up log file and announce end
    if (rank == 0 && log_file) fclose(log_file);
    if (rank == 0)
        printf(">> End of EPO\n");
}

/**
 * Main entry point: initializes MPI, distributes work, and invokes EPO.
 *
 * Usage: ./program <city_file> <num_cities> <num_iterations>
 */
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Check command-line arguments
    if (argc != 4) {
        if (rank == 0)
            printf("Usage: %s <city_file> <city_number> <iteration_number>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    // Seed random generator uniquely per process
    srand(time(NULL) + rank);

    const char *filename = argv[1];
    int num_cities = atoi(argv[2]);
    int num_iterations = atoi(argv[3]);

    // Determine number of penguins based on problem size
    int num_penguins = num_cities / 100;
    if (num_penguins < 5) num_penguins = 5;
    if (num_penguins > MAX_PENGUINS) {
        if (rank == 0)
            printf("Error: Number of penguins too large. Max allowed: %d\n", MAX_PENGUINS);
        MPI_Finalize();
        return 1;
    }

     // Read city data on root and broadcast to all
    City cities[MAX_CITIES];
    if (rank == 0 && !read_cities(filename, cities, num_cities)) {
        printf("Error reading cities.\n");
        MPI_Finalize();
        return 1;
    }

    MPI_Bcast(cities, sizeof(City) * num_cities, MPI_BYTE, 0, MPI_COMM_WORLD);

    // Allocate storage for the final best path
    City best_path[MAX_CITIES];

    // Measure execution time
    double start_time = MPI_Wtime();
    run_epo(cities, num_cities, num_penguins, num_iterations, best_path, rank, size);
    double end_time = MPI_Wtime();

    if (rank == 0) {
        double best_length = calculate_path_length(best_path, num_cities);
        printf("\n>> Best route calculation completed\n");
        printf("\nBest path found (length: %.2f)\n", best_length);
        printf("\n>> Execution time: %.3f seconds\n", end_time - start_time);
    }

    MPI_Finalize();
    return 0;
}
