#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void generate_cities(int num_cities, const char* filename) {
    
    FILE *file = fopen(filename, "w");
    if (!file) {
        printf("Errore nell'aprire il file %s\n", filename);
        exit(1);
    }

    
    for (int i = 0; i < num_cities; i++) {
        double x = (rand() % 10000) / 100.0; 
        double y = (rand() % 10000) / 100.0;
        fprintf(file, "%.2f %.2f\n", x, y); 
    }

    fclose(file);
}

int main() {
    srand(time(NULL));  
    int num_cities;
    printf("Enter the number of cities to generate: ");
    scanf("%d", &num_cities); 

    const char* filename = "cities.txt";  
    generate_cities(num_cities, filename);

    printf("Cities have been generated and saved in %s\n", filename);

    return 0;
}
