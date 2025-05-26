#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#define MAX_CITIES 16000 // Maximum number of cities
#define MAX_PENGUINS 500 // Maximum number of penguin agents

typedef struct {
    int id;  // City identifier
    double x, y; // Coordinates of the city
} City;

/**
 * Compute the Euclidean distance between two cities.
 *
 * @param a First city
 * @param b Second city
 * @return Straight-line distance between a and b
 */
double calculate_distance(City a, City b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return sqrt(dx * dx + dy * dy);
}

/**
 * Compute the total length of a closed tour through all cities.
 * Adds distance from last back to first.
 *
 * @param path Array of cities in visitation order
 * @param num_cities Number of cities
 * @return Total tour length
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
 * Load city coordinates from a text file.
 *
 * @param filename Input file, each line: "x y"
 * @param cities Preallocated array to fill
 * @param num_cities Expected number of cities
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
    return 1;
}

/**
 * Randomly shuffle the order of cities for initial path generation.
 *
 * @param path Array of cities
 * @param num_cities Length of array
 */
void shuffle_path(City path[], int num_cities) {
    // Perform 2*num_cities random swaps
    for (int i = 0; i < num_cities * 2; i++) {
        int a = rand() % num_cities;
        int b = rand() % num_cities;
        City temp = path[a];
        path[a] = path[b];
        path[b] = temp;
    }
}

/**
 * Copy a central 80% segment from the leader's tour into a follower,
 * then fill remaining positions with unused cities in random order.
 *
 * @param leader Global best path
 * @param follower Path to update
 * @param num_cities Number of cities
 * @param start_percent Fractional start index for the segment (0.0 to 0.2)
 */
void copy_leader_segment(City leader[], City follower[], int num_cities, double start_percent) {
    int segment_size = num_cities * 0.8;
    int start = (int)(num_cities * start_percent);
    int end = start + segment_size;
    if (end > num_cities) end = num_cities;

    int used[MAX_CITIES] = {0}; // Flags for cities already placed

    // Copy segment from leader
    for (int i = start; i < end; i++) {
        follower[i] = leader[i];
        // Mark this city used by matching coordinates
        for (int j = 0; j < num_cities; j++) {
            if (leader[i].x == leader[j].x && leader[i].y == leader[j].y) {
                used[j] = 1;
                break;
            }
        }
    }

    // Fill rest with unused cities
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
 * Early iterations have more swaps; later ones fewer.
 *
 * @param path Array to mutate
 * @param num_cities Number of cities
 * @param iter Current iteration index
 * @param max_iter Total iterations
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
 * Hybrid EPO driver: each MPI rank handles a subset of penguins,
 * and within each rank OpenMP threads evaluate and evolve paths.
 *
 * @param cities Global array of city coordinates (broadcasted to all ranks)
 * @param num_cities Number of cities in problem
 * @param num_penguins Total agents across all ranks
 * @param num_iterations Optimization iterations
 * @param best_path Output best path found globally
 * @param rank MPI rank of this process
 * @param size Total MPI ranks
 */
void run_epo(City cities[], int num_cities, int num_penguins, int num_iterations, City best_path[], int rank, int size) {
    if (rank == 0)
        printf(">> Start EPO with %d penguins, %d cities, %d iterations\n", num_penguins, num_cities, num_iterations);

    FILE *log_file = NULL;
    if (rank == 0) {
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

    static City population[MAX_PENGUINS][MAX_CITIES]; // Local population

    // Initialize random paths for local penguins
    for (int i = 0; i < local_penguins; i++) {
        for (int j = 0; j < num_cities; j++) {
            population[i][j] = cities[j];
        }
        shuffle_path(population[i], num_cities);
    }

    if (rank == 0)
        printf(">> Generate random initial paths for penguins\n");

    double best_length = INFINITY;
    int last_improve = 0;
    City leader_global[MAX_CITIES];

    // Main iteration loop
    for (int iter = 0; iter < num_iterations; iter++) {

        // Stop criterion
        if (iter - last_improve >= 100) {
            if (rank == 0)
                printf(">> Stop: no improvement after %d iterations, exiting at step %d\n",
                       iter - last_improve, iter);
            break;
        }


        if (rank == 0)
            printf(">> Iteration %d/%d\n", iter + 1, num_iterations);

        // Local best search using OpenMP
        int local_leader_index = 0;
        double local_best = INFINITY;

        #pragma omp parallel
        {
            int thread_leader = -1;
            double thread_best = INFINITY;

            #pragma omp for
            for (int i = 0; i < local_penguins; i++) {
                double len = calculate_path_length(population[i], num_cities);
                #pragma omp critical
                {
                    if (rank == 0 && log_file)
                        fprintf(log_file, "%d,%d,%.2f\n", iter + 1, i + rank * local_penguins, len);
                }

                if (len < thread_best) {
                    thread_best = len;
                    thread_leader = i;
                }
            }

            #pragma omp critical
            {
                if (thread_best < local_best) {
                    local_best = thread_best;
                    local_leader_index = thread_leader;
                }
            }
        }

        // Copy local leader's path
        City local_leader[MAX_CITIES];
        for (int i = 0; i < num_cities; i++)
            local_leader[i] = population[local_leader_index][i];

        // MPI reduction to find global best
        struct {
            double len;
            int rank;
        } in = { local_best, rank }, out;

        MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

        // Root of out rank shares its leader path
        if (rank == out.rank) {
            for (int i = 0; i < num_cities; i++)
                leader_global[i] = local_leader[i];
        }

        MPI_Bcast(leader_global, sizeof(City) * num_cities, MPI_BYTE, out.rank, MPI_COMM_WORLD);

        if (rank == 0)
            printf("   >> Leader: Penguin %d (from process %d) with length %.2f\n",
                   local_leader_index + out.rank * local_penguins, out.rank, out.len);

        // Update best_path if improved
        if (out.len < best_length) {
            best_length = out.len;
            last_improve = iter;
            for (int i = 0; i < num_cities; i++)
                best_path[i] = leader_global[i];
        }

        // Followers copy & mutate in parallel
        #pragma omp parallel for
        for (int i = 0; i < local_penguins; i++) {
            if (i == local_leader_index && rank == out.rank) continue;
            double start_percent = (i % 5) * 0.2;
            copy_leader_segment(leader_global, population[i], num_cities, start_percent);
            mutate_path(population[i], num_cities, iter, num_iterations);
        }
    }

    if (rank == 0 && log_file) fclose(log_file);
    if (rank == 0)
        printf(">> End of EPO\n");
}

/**
 * Entry point: initialize MPI, read input, and launch hybrid EPO.
 *
 * Usage: %s <city_file> <num_cities> <num_iterations>
 */
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 4) {
        if (rank == 0)
            printf("Usage: %s <city_file> <city_number> <iteration_number>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    // Display MPI and OpenMP environment
    printf(">> [MPI %d/%d] OpenMP: %d threads available\n", rank, size, omp_get_max_threads());

    srand(time(NULL) + rank);

    const char *filename = argv[1];
    int num_cities = atoi(argv[2]);
    int num_iterations = atoi(argv[3]);
    int num_penguins = num_cities / 100;
    if (num_penguins < 5) num_penguins = 5;
    if (num_penguins > MAX_PENGUINS) {
        if (rank == 0)
            printf("Error: Number of penguins too large. Max allowed: %d\n", MAX_PENGUINS);
        MPI_Finalize();
        return 1;
    }

    // Read city data on root and broadcast
    City cities[MAX_CITIES];
    if (rank == 0 && !read_cities(filename, cities, num_cities)) {
        printf("Error reading cities.\n");
        MPI_Finalize();
        return 1;
    }

    MPI_Bcast(cities, sizeof(City) * num_cities, MPI_BYTE, 0, MPI_COMM_WORLD);

    // Storage for best global path
    City best_path[MAX_CITIES];

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
