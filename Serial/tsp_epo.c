#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define NUM_CITIES 4000
#define NUM_PENGUINS 50
#define MAX_ITERATIONS 500

// struttura città
typedef struct {
    double x, y;
} City;

// distanza tra due città
double distance(City a, City b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

double calculate_path_length(City cities[], int path[], int num_cities) {
    double total_distance = 0.0;
    for (int i = 0; i < num_cities - 1; i++) {
        total_distance += distance(cities[path[i]], cities[path[i + 1]]);
    }
    total_distance += distance(cities[path[num_cities - 1]], cities[path[0]]);
    return total_distance;
}

// genera percorso casuale
void generate_random_path(int path[], int num_cities) {
    for (int i = 0; i < num_cities; i++) {
        path[i] = i;
    }
    for (int i = 0; i < num_cities; i++) {
        int j = rand() % num_cities;
        int temp = path[i];
        path[i] = path[j];
        path[j] = temp;
    }
}

void copy_path(int dest[], int src[], int num_cities) {
    for (int i = 0; i < num_cities; i++) {
        dest[i] = src[i];
    }
}

// Funzione per verificare se una città è già nel percorso
int contains(int* path, int city, int num_cities) {
    for (int i = 0; i < num_cities; i++) {
        if (path[i] == city) {
            return 1;
        }
    }
    return 0; 
}

// Crossover tra follower e leader
void crossover(int* follower_path, int* leader_path, int num_cities, int start, int end) {
    int* new_path = (int*)malloc(num_cities * sizeof(int));
    int i, index = 0;
    memset(new_path, -1, num_cities * sizeof(int)); // Inizializza con -1 (città non assegnata)

    // Copia il segmento dal percorso del leader
    for (i = start; i <= end; i++) {
        new_path[i] = leader_path[i];
    }

    // Riempi le posizioni rimanenti del follower con le città mancanti
    for (i = 0; i < num_cities; i++) {
        if (new_path[i] == -1) {
            while (contains(new_path, follower_path[index], num_cities)) {
                index++;
            }
            new_path[i] = follower_path[index];
            index++;
        }
    }

    copy_path(follower_path, new_path, num_cities);
    free(new_path);
}

// Emperor Penguin Optimization
void emperor_penguin_optimization(City cities[], int num_cities, int num_penguins, int max_iterations) {
    int paths[NUM_PENGUINS][NUM_CITIES]; 
    double fitness[NUM_PENGUINS];        
    int best_path[NUM_CITIES];           
    double best_fitness = INFINITY;      

    // Inizializza la popolazione di pinguini
    for (int i = 0; i < num_penguins; i++) {
        generate_random_path(paths[i], num_cities);
        fitness[i] = calculate_path_length(cities, paths[i], num_cities);
        if (fitness[i] < best_fitness) {
            best_fitness = fitness[i];
            copy_path(best_path, paths[i], num_cities);
        }
    }

    for (int iter = 0; iter < max_iterations; iter++) {
        // Dividi i pinguini in migliori e follower
        for (int i = 0; i < num_penguins; i++) {
            if (fitness[i] < best_fitness) {
                best_fitness = fitness[i];
                copy_path(best_path, paths[i], num_cities);
            }
        }

        // Aggiorna i follower
        for (int i = 0; i < num_penguins; i++) {
            if (fitness[i] > best_fitness) {
                // Follower: prendi una parte del percorso del leader e mescolalo con il proprio
                int start = rand() % num_cities;
                int end = start + rand() % (num_cities - start); 
                crossover(paths[i], best_path, num_cities, start, end);
                
                // Calcola la nuova fitness
                fitness[i] = calculate_path_length(cities, paths[i], num_cities);
                if (fitness[i] < best_fitness) {
                    best_fitness = fitness[i];
                    copy_path(best_path, paths[i], num_cities);
                }
            }
        }
         // Stampa il miglior percorso trovato ad ogni iterazione
        printf("Iterazione %d: Miglior Percorso Lunghezza = %.2f\n", iter + 1, best_fitness);
    }

    // Stampa il risultato finale
    printf("Lunghezza del Percorso Ottimale: %.2f\n", best_fitness);
}

int main() {
    srand(time(NULL));

    // Genera un insieme di città casuali
    City cities[NUM_CITIES];
    for (int i = 0; i < NUM_CITIES; i++) {
        cities[i].x = rand() % 100;
        cities[i].y = rand() % 100;
    }

    // Misura il tempo di esecuzione
    clock_t start = clock();

    // Applica l'algoritmo EPO al problema del TSP
    emperor_penguin_optimization(cities, NUM_CITIES, NUM_PENGUINS, MAX_ITERATIONS);

    clock_t end = clock();
    printf("Tempo di Esecuzione: %.2f secondi\n", (double)(end - start) / CLOCKS_PER_SEC);

    return 0;
}
