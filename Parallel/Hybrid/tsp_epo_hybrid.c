#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#define MAX_CITIES 16000
#define MAX_PENGUINS 500

typedef struct {
    int id;
    double x, y;
} City;

double calculate_distance(City a, City b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return sqrt(dx * dx + dy * dy);
}

double calculate_path_length(City path[], int num_cities) {
    double total = 0.0;
    for (int i = 0; i < num_cities - 1; i++) {
        total += calculate_distance(path[i], path[i + 1]);
    }
    total += calculate_distance(path[num_cities - 1], path[0]);
    return total;
}

int read_cities(const char *filename, City cities[], int num_cities) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Errore: impossibile aprire il file '%s'\n", filename);
        return 0;
    }

    for (int i = 0; i < num_cities; i++) {
        if (fscanf(file, "%lf %lf", &cities[i].x, &cities[i].y) != 2) {
            printf("Errore: città %d non letta correttamente\n", i);
            fclose(file);
            return 0;
        }
    }
    fclose(file);
    return 1;
}

void shuffle_path(City path[], int num_cities) {
    for (int i = 0; i < num_cities * 2; i++) {
        int a = rand() % num_cities;
        int b = rand() % num_cities;
        City temp = path[a];
        path[a] = path[b];
        path[b] = temp;
    }
}

void copy_leader_segment(City leader[], City follower[], int num_cities, double start_percent) {
    int segment_size = num_cities * 0.8;
    int start = (int)(num_cities * start_percent);
    int end = start + segment_size;
    if (end > num_cities) end = num_cities;

    int used[MAX_CITIES] = {0};

    for (int i = start; i < end; i++) {
        follower[i] = leader[i];
        for (int j = 0; j < num_cities; j++) {
            if (leader[i].x == leader[j].x && leader[i].y == leader[j].y) {
                used[j] = 1;
                break;
            }
        }
    }

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

void run_epo(City cities[], int num_cities, int num_penguins, int num_iterations, City best_path[], int rank, int size) {
    if (rank == 0)
        printf(">> Inizio EPO con %d pinguini, %d città, %d iterazioni\n", num_penguins, num_cities, num_iterations);

    FILE *log_file = NULL;
    if (rank == 0) {
        log_file = fopen("epo_log.csv", "w");
        if (!log_file) {
            printf("Errore: impossibile creare il file di log.\n");
            return;
        }
        fprintf(log_file, "iterazione,pinguino,lunghezza\n");
    }

    int local_penguins = num_penguins / size;
    if (rank < num_penguins % size) local_penguins++;

    static City population[MAX_PENGUINS][MAX_CITIES];

    for (int i = 0; i < local_penguins; i++) {
        for (int j = 0; j < num_cities; j++) {
            population[i][j] = cities[j];
        }
        shuffle_path(population[i], num_cities);
    }

    if (rank == 0)
        printf(">> Generati percorsi iniziali casuali per i pinguini\n");

    double best_length = INFINITY;
    City leader_global[MAX_CITIES];

    for (int iter = 0; iter < num_iterations; iter++) {
        if (rank == 0)
            printf(">> Iterazione %d/%d\n", iter + 1, num_iterations);

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

        City local_leader[MAX_CITIES];
        for (int i = 0; i < num_cities; i++)
            local_leader[i] = population[local_leader_index][i];

        struct {
            double len;
            int rank;
        } in = { local_best, rank }, out;

        MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

        if (rank == out.rank) {
            for (int i = 0; i < num_cities; i++)
                leader_global[i] = local_leader[i];
        }

        MPI_Bcast(leader_global, sizeof(City) * num_cities, MPI_BYTE, out.rank, MPI_COMM_WORLD);

        if (rank == 0)
            printf("   >> Leader: Pinguino %d (da processo %d) con lunghezza %.2f\n",
                   local_leader_index + out.rank * local_penguins, out.rank, out.len);

        if (out.len < best_length) {
            best_length = out.len;
            for (int i = 0; i < num_cities; i++)
                best_path[i] = leader_global[i];
        }

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
        printf(">> Fine EPO\n");
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 4) {
        if (rank == 0)
            printf("Uso: %s <file_citta> <numero_citta> <numero_iterazioni>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    printf(">> [MPI %d/%d] OpenMP: %d thread disponibili\n", rank, size, omp_get_max_threads());

    srand(time(NULL) + rank);

    const char *filename = argv[1];
    int num_cities = atoi(argv[2]);
    int num_iterations = atoi(argv[3]);
    int num_penguins = num_cities / 100;
    if (num_penguins < 5) num_penguins = 5;
    if (num_penguins > MAX_PENGUINS) {
        if (rank == 0)
            printf("Errore: numero pinguini troppo grande. Max consentito: %d\n", MAX_PENGUINS);
        MPI_Finalize();
        return 1;
    }

    City cities[MAX_CITIES];
    if (rank == 0 && !read_cities(filename, cities, num_cities)) {
        printf("Errore nella lettura delle città.\n");
        MPI_Finalize();
        return 1;
    }

    MPI_Bcast(cities, sizeof(City) * num_cities, MPI_BYTE, 0, MPI_COMM_WORLD);

    City best_path[MAX_CITIES];

    double start_time = MPI_Wtime();
    run_epo(cities, num_cities, num_penguins, num_iterations, best_path, rank, size);
    double end_time = MPI_Wtime();

    if (rank == 0) {
        double best_length = calculate_path_length(best_path, num_cities);
        printf("\n>> Calcolo percorso migliore terminato\n");
        printf("\nMiglior percorso trovato (lunghezza: %.2f)\n", best_length);
        printf("\n>> Tempo di esecuzione: %.3f secondi\n", end_time - start_time);
    }

    MPI_Finalize();
    return 0;
}
