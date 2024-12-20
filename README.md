# Parallel-Emperor-Penguin-Optimization-Algorithm

Cities are generated via the `generate_cities.c` file, which creates a `cities.txt` file with one city per line.

To generate cities, run the following commands:
```bash
gcc generate_cities.c -o generate_cities
./generate_cities.exe
```

Then you can use EPO to fix the TSP problem with the following commands:
```bash
gcc tsp_epo.c -o tsp_epo
.\tsp_epo.exe
```

## Parallelizzazione
Inizialmente usiamo un numero di città fisse (4000) e di pinguini fissi (50), cerchiamo di parallelizzare usando MPI e dividendo la popolazione di pinguini tra i vari processi.
Successivamente implementaremo anche OpenMP all'inerno di ogni processo per parallelizzare i calcoli di aggiornamento dei pinguini tramite thread.
Infine useremo dataset di varie dimensioni (numero di città) per vedere le performance con vari numeri di processi
