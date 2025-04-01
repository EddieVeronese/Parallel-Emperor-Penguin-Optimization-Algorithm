#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <time.h>

#define MAX_CITIES 10000  // Numero massimo di città nel dataset (leggi dal file)
#define MAX_ITERATIONS 300  // Numero massimo di iterazioni
#define NUM_PENGUINS 100  // Numero di pinguini nella popolazione

// Dati globali
int num_cities;
double distances[MAX_CITIES][MAX_CITIES];
int path[MAX_CITIES];
int best_path[MAX_CITIES];
double best_distance = DBL_MAX;

// Funzione per calcolare la distanza euclidea tra due città
double euclidean_distance(int i, int j, double cities[][2]) {
    return sqrt(pow(cities[i][0] - cities[j][0], 2) + pow(cities[i][1] - cities[j][1], 2));
}

// Funzione per calcolare la lunghezza totale di un percorso
double path_length(int *path) {
    double total_distance = 0.0;
    for (int i = 0; i < num_cities - 1; i++) {
        total_distance += distances[path[i]][path[i + 1]];
    }
    total_distance += distances[path[num_cities - 1]][path[0]];  // ritorno al punto di partenza
    return total_distance;
}

// Funzione per calcolare le distanze tra tutte le città (pre-calcolo)
void calculate_distances(double cities[][2]) {
    for (int i = 0; i < num_cities; i++) {
        for (int j = i + 1; j < num_cities; j++) {
            distances[i][j] = euclidean_distance(i, j, cities);
            distances[j][i] = distances[i][j];
        }
    }
}

// Funzione di ottimizzazione locale (2-opt)
void local_search(int *path) {
    int i, j;
    double best_len = path_length(path);
    for (i = 0; i < num_cities - 1; i++) {
        for (j = i + 1; j < num_cities; j++) {
            if (i != j) {
                // Scambia i due segmenti
                int temp = path[i];
                path[i] = path[j];
                path[j] = temp;

                double new_len = path_length(path);
                if (new_len < best_len) {
                    best_len = new_len;
                    for (int k = 0; k < num_cities; k++) {
                        best_path[k] = path[k];
                    }
                } else {
                    // Ripristina l'ordine
                    path[j] = path[i];
                    path[i] = temp;
                }
            }
        }
    }
}

// Funzione per leggere le città da un file
int read_cities_from_file(const char *filename, double cities[][2], int max_cities) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error opening file %s\n", filename);
        return -1;
    }

    int i = 0;
    while (fscanf(file, "%lf %lf", &cities[i][0], &cities[i][1]) == 2 && i < max_cities) {
        i++;
    }
    fclose(file);
    return i;
}

// Funzione per inizializzare il percorso (ordine casuale delle città)
void initialize_path(int *path) {
    for (int i = 0; i < num_cities; i++) {
        path[i] = i;
    }

    // Shuffle per ottenere un percorso casuale
    for (int i = 0; i < num_cities; i++) {
        int j = rand() % num_cities;
        int temp = path[i];
        path[i] = path[j];
        path[j] = temp;
    }
}

// Funzione per raccogliere i migliori risultati con MPI
void gather_results(int rank, double *best_distance, int *best_path) {
    double global_best_distance;
    int global_best_path[MAX_CITIES];

    MPI_Reduce(best_distance, &global_best_distance, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        *best_distance = global_best_distance;
        MPI_Gather(best_path, num_cities, MPI_INT, global_best_path, num_cities, MPI_INT, 0, MPI_COMM_WORLD);
    }
}

// Funzione principale per il calcolo parallelo del TSP
void emperor_penguin_optimization(int rank, int size, double cities[][2]) {
    // Distribuzione iniziale del lavoro (ogni processo lavora su un sottoinsieme dei pinguini)
    int local_best_path[MAX_CITIES];
    double local_best_distance = DBL_MAX;
    int global_best_path[MAX_CITIES];
    double global_best_distance;

    // Inizializzazione dei pinguini
    initialize_path(path);
    local_best_distance = path_length(path);
    for (int i = 0; i < num_cities; i++) {
        local_best_path[i] = path[i];
    }

    // Esegui la ricerca locale (2-opt)
    local_search(local_best_path);
    local_best_distance = path_length(local_best_path);

    // Sincronizzazione e raccolta dei migliori risultati
    gather_results(rank, &local_best_distance, local_best_path);

    // Aggiornamento globale
    if (rank == 0) {
        global_best_distance = local_best_distance;
        for (int i = 0; i < num_cities; i++) {
            global_best_path[i] = local_best_path[i];
        }
    }
    MPI_Bcast(&global_best_distance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(global_best_path, num_cities, MPI_INT, 0, MPI_COMM_WORLD);

    // Stampa del miglior percorso trovato
    if (rank == 0) {
        printf("Best path found by EPO (Length: %.2f):\n", global_best_distance);
        for (int i = 0; i < num_cities; i++) {
            printf("%d -> ", global_best_path[i]);
        }
        printf("%d\n", global_best_path[0]);  // Ritorno alla città iniziale
    }
}

// Funzione per inizializzare MPI e avviare l'algoritmo
int main(int argc, char **argv) {
    double cities[MAX_CITIES][2];
    int rank, size;
    
    // Inizializzazione di MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Lettura del numero di città da analizzare (parametro di input)
    if (argc < 2) {
        printf("Usage: %s <number_of_cities_to_analyze>\n", argv[0]);
        MPI_Finalize();
        return -1;
    }

    num_cities = atoi(argv[1]); // Numero di città da analizzare

    // Assicurati che il numero di città non ecceda il massimo disponibile
    if (num_cities > MAX_CITIES) {
        printf("Error: Cannot analyze more than %d cities.\n", MAX_CITIES);
        MPI_Finalize();
        return -1;
    }

    // Genera città da file
    if (read_cities_from_file("cities.txt", cities, MAX_CITIES) < 0) {
        MPI_Finalize();
        return -1;
    }

    // Calcola le distanze tra le città
    calculate_distances(cities);

    // Esegui l'Emperor Penguin Optimization (EPO)
    emperor_penguin_optimization(rank, size, cities);

    // Finalizzazione di MPI
    MPI_Finalize();

    return 0;
}
