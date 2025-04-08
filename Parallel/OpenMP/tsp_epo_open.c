#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>  // <-- OpenMP

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
    printf(">> Letti %d città da '%s'\n", num_cities, filename);
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

void run_epo(City cities[], int num_cities, int num_penguins, int num_iterations, City best_path[]) {
    printf(">> Inizio EPO con %d pinguini, %d città, %d iterazioni\n", num_penguins, num_cities, num_iterations);

    FILE *log_file = fopen("epo_log.csv", "w");
    if (!log_file) {
        printf("Errore: impossibile creare il file di log.\n");
        return;
    }
    fprintf(log_file, "iterazione,pinguino,lunghezza\n");

    static City population[MAX_PENGUINS][MAX_CITIES];

    // Parallel init
    #pragma omp parallel for
    for (int i = 0; i < num_penguins; i++) {
        for (int j = 0; j < num_cities; j++) {
            population[i][j] = cities[j];
        }
        shuffle_path(population[i], num_cities);
    }
    printf(">> Generati percorsi iniziali casuali per i pinguini\n");

    double best_length = INFINITY;

    for (int iter = 0; iter < num_iterations; iter++) {
        printf(">> Iterazione %d/%d\n", iter + 1, num_iterations);
        int leader_index = 0;
        double leader_length = INFINITY;

        #pragma omp parallel
        {
            int local_leader = -1;
            double local_best = INFINITY;

            #pragma omp for nowait
            for (int i = 0; i < num_penguins; i++) {
                double current_length = calculate_path_length(population[i], num_cities);

                #pragma omp critical
                fprintf(log_file, "%d,%d,%.2f\n", iter + 1, i, current_length);

                if (current_length < local_best) {
                    local_best = current_length;
                    local_leader = i;
                }
            }

            #pragma omp critical
            {
                if (local_best < leader_length) {
                    leader_length = local_best;
                    leader_index = local_leader;
                }
            }
        }

        printf("   >> Leader: Pinguino %d con lunghezza %.2f\n", leader_index, leader_length);

        if (leader_length < best_length) {
            best_length = leader_length;
            for (int i = 0; i < num_cities; i++) {
                best_path[i] = population[leader_index][i];
            }
        }

        #pragma omp parallel for
        for (int i = 0; i < num_penguins; i++) {
            if (i == leader_index) continue;
            double start_percent = (i % 5) * 0.2;
            copy_leader_segment(population[leader_index], population[i], num_cities, start_percent);
            mutate_path(population[i], num_cities, iter, num_iterations);
        }
    }

    fclose(log_file);
    printf(">> Fine EPO\n");
}

int main(int argc, char *argv[]) {
    printf(">> OpenMP: %d thread disponibili\n", omp_get_max_threads());

    if (argc != 4) {
        printf("Uso: %s <file_citta> <numero_citta> <numero_iterazioni>\n", argv[0]);
        return 1;
    }

    srand(time(NULL));

    const char *filename = argv[1];
    int num_cities = atoi(argv[2]);
    int num_iterations = atoi(argv[3]);
    int num_penguins = num_cities / 100;
    if (num_penguins < 5) num_penguins = 5;

    if (num_cities > MAX_CITIES) {
        printf("Errore: numero città troppo grande. Max consentito: %d\n", MAX_CITIES);
        return 1;
    }
    if (num_penguins > MAX_PENGUINS) {
        printf("Errore: numero pinguini troppo grande. Max consentito: %d\n", MAX_PENGUINS);
        return 1;
    }

    City cities[MAX_CITIES];
    if (!read_cities(filename, cities, num_cities)) {
        printf("Errore nella lettura delle città.\n");
        return 1;
    }

    City best_path[MAX_CITIES];

    // TIMER START
    double start_time = omp_get_wtime();

    run_epo(cities, num_cities, num_penguins, num_iterations, best_path);

    // TIMER END
    double end_time = omp_get_wtime();

    printf("\n>> Calcolo percorso migliore terminato\n");
    double best_length = calculate_path_length(best_path, num_cities);
    printf("\nMiglior percorso trovato (lunghezza: %.2f)\n", best_length);

    // STAMPA TEMPO
    printf("\n>> Tempo di esecuzione: %.3f secondi\n", end_time - start_time);

    return 0;
}
