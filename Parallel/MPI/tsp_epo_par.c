#include "tsp_epo_par.h"

/*
    Calculate the Euclidean distance between two cities
    Parameters:
    - a: First city
    - b: Second city
    Return:
    - The Euclidean distance between the two points
*/
double distance(City a, City b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

/*
    Calculate the total length of a penguin's path
    Parameters:
    - cities: Array of cities
    - penguin: Pointer to the penguin
    - num_cities: Number of cities in the route
    Return:
    - The total length of the route
*/
double calculate_path_length(City cities[], Penguin* penguin, int num_cities) {
    double total_distance = 0.0;
    for (int i = 0; i < num_cities - 1; i++) {
        total_distance += distance(cities[penguin->path[i]], cities[penguin->path[i + 1]]);
    }
    // Distance between the last city visited and the city of departure
    total_distance += distance(cities[penguin->path[num_cities - 1]], cities[penguin->path[0]]);
    return total_distance;
}

/*
    Generate a random path for a penguin
    Parameters:
    - penguin: Pointer to the penguin
    - num_cities: Total number of cities
*/
void generate_random_path(Penguin* penguin, int num_cities) {
    for (int i = 1; i < num_cities; i++) {
        penguin->path[i] = i;
    }
    for (int i = 1; i < num_cities; i++) {
        int j = 1 + rand() % (num_cities - 1); 
        int temp = penguin->path[i];
        penguin->path[i] = penguin->path[j];
        penguin->path[j] = temp;
    }
}

/*
    Check if a city is already present in a route
    Parameters:
    - path: Array of the path
    - city: City to check
    - num_cities: Total number of cities
    Return:
    - 1 if the city is present, 0 otherwise
*/
int contains(int* path, int city, int num_cities) {
    for (int i = 0; i < num_cities; i++) {
        if (path[i] == city) {
            return 1;
        }
    }
    return 0; 
}

/*
    Performs improved crossover between a leader and a follower
    Parameters:
    - follower: Pointer to the follower penguin
    - leader: Pointer to the leader penguin
    - num_cities: Total number of cities
*/
void enhanced_crossover(Penguin* follower, Penguin* leader, int num_cities) {
    int segments = 3; 
    for (int s = 0; s < segments; s++) {
        int start = rand() % num_cities;
        int end = start + rand() % (num_cities - start);

        for (int i = start; i <= end; i++) {
            follower->path[i] = leader->path[i];
        }
    }

    int index = 0;
    for (int i = 0; i < num_cities; i++) {
        if (!contains(follower->path, leader->path[i], num_cities)) {
            while (contains(follower->path, follower->path[index], num_cities)) {
                index++;
            }
            follower->path[index] = leader->path[i];
        }
    }
}


/*
    Performs a random mutation on a penguin
    Parameters:
    - penguin: Pointer to the penguin
    - num_cities: Total number of cities
*/
void mutate(Penguin* penguin, int num_cities) {
    double mutation_rate = 0.2; 
    if ((rand() / (double)RAND_MAX) < mutation_rate) {
        int i = rand() % num_cities; 
        int j = rand() % num_cities; 
        int temp = penguin->path[i];
        penguin->path[i] = penguin->path[j];
        penguin->path[j] = temp;
    }
}

/*
    Optimize the route using the EPO algorithm
    Parameters:
    - cities: Array of cities
    - num_cities: Total number of cities
    - num_penguins: Total number of penguins
    - max_iterations: Maximum number of iterations
    - rank: Rank of the MPI process
    - size: Total number of MPI processes
*/
void emperor_penguin_optimization(City cities[], int num_cities, int num_penguins, int max_iterations, int rank, int size) {

    FILE* file = NULL;
    if (rank == 0) {
        file = fopen("optimal_paths.csv", "w");
        if (!file) {
            printf("Error opening file for writing.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    Penguin penguins[NUM_PENGUINS / size];    
    double best_fitness = INFINITY;     
    Penguin best_penguin;             

    //Generate a complete path for each penguin
    for (int i = 0; i < num_penguins / size; i++) {
        generate_random_path(&penguins[i], num_cities);
        //Calculate the length of each path
        penguins[i].cost = calculate_path_length(cities, &penguins[i], num_cities);
        //Find the best path and copy it into best_penguin
        if (penguins[i].cost < best_fitness) {
            best_fitness = penguins[i].cost;
            best_penguin = penguins[i];
        }
    }

    //Counter for iterations without improvement
    int iterations_without_improvement = 0;  

    //Execute for each iteration
    for (int iter = 0; iter < max_iterations; iter++) {
        //Recalculate the best penguin
        int best_found_in_iteration = 0; // Flag to track if improvement is found in the current iteration
        for (int i = 0; i < num_penguins / size; i++) {
            if (penguins[i].cost < best_fitness) {
                best_fitness = penguins[i].cost;
                best_penguin = penguins[i];
                best_found_in_iteration = 1;  // Mark that an improvement was found
            }
        }

        //Update followers
         for (int i = 0; i < num_penguins / size; i++) {
            // If the length is greater than the optimal one, change part of the path
            if (penguins[i].cost > best_fitness) {
                int start = rand() % num_cities;
                int end = start + rand() % (num_cities - start); 
                //crossover(&penguins[i], &best_penguin, num_cities, start, end);
                enhanced_crossover(&penguins[i], &best_penguin, num_cities);
                mutate(&penguins[i], num_cities);

                //Calculate the new fitness
                penguins[i].cost = calculate_path_length(cities, &penguins[i], num_cities);
                //If the length is less than optimal, update the best route
                if (penguins[i].cost < best_fitness) {
                    best_fitness = penguins[i].cost;
                    best_penguin = penguins[i];
                    best_found_in_iteration = 1;  
                }
            }
        }

        //Increment the counter if no improvement was found in this iteration
        if (!best_found_in_iteration) {
            iterations_without_improvement++;
        } else {
            iterations_without_improvement = 0;  
        }

        //If no improvement has been made in 30 iterations, stop the algorithm
        
        if (iterations_without_improvement >= 50) {
            printf("No improvement after 50 iterations, stopping...\n");
            break;
        }
        
        

        //Print the best route at the end of the iteration
        /*
        if (rank == 0) {
            printf("Iteration %d: Best Path Length = %.2f\n", iter + 1, best_fitness);
        } */    

        if (rank == 0) {
            fprintf(file, "%d,%.2f\n", iter + 1, best_fitness);
        }



        MPI_Barrier(MPI_COMM_WORLD);
        double global_best_fitness;
        Penguin global_best_penguin;
        MPI_Allreduce(&best_fitness, &global_best_fitness, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

        if (best_fitness == global_best_fitness) {
            global_best_penguin = best_penguin;
        }

        MPI_Bcast(&global_best_penguin, sizeof(Penguin), MPI_BYTE, 0, MPI_COMM_WORLD);
        best_fitness = global_best_fitness;
        best_penguin = global_best_penguin;
    }

    //Print the best route at the end of all the iterations
    if (rank == 0) {
        printf("Optimal Path Length: %.2f\n", best_fitness);
        fclose(file);
    }
}



/*
    Reads cities from a file
    Parameters:
    - cities: Array of cities
    - filename: Name of the file
    - num_cities: Total number of cities to read
*/
void read_cities_from_file(City cities[], const char* filename, int num_cities) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }
    for (int i = 0; i < num_cities; i++) {
        fscanf(file, "%lf %lf", &cities[i].x, &cities[i].y);
    }
    fclose(file);
}

int main(int argc, char** argv) {
    srand(time(NULL));

    MPI_Init(&argc, &argv);
    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //Read cities from file
    City cities[NUM_CITIES];
    if (rank == 0) {
        read_cities_from_file(cities, "cities.txt", NUM_CITIES);
    }

    MPI_Bcast(cities, NUM_CITIES * sizeof(City), MPI_BYTE, 0, MPI_COMM_WORLD);

    double start_time, end_time;

    start_time = MPI_Wtime();

    //Apply the EPO algorithm to the TSP problem
    emperor_penguin_optimization(cities, NUM_CITIES, NUM_PENGUINS, MAX_ITERATIONS, rank, size);

    end_time = MPI_Wtime();

    double elapsed_time = end_time - start_time;
    double max_time;

    MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Local Execution Time (Rank 0): %.3f seconds\n", elapsed_time);
        printf("Maximum Execution Time: %.3f seconds\n", max_time);
    }

    MPI_Finalize();

    return 0;
}
