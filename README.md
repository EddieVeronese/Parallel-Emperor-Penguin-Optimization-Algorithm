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
