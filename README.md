# Parallel Emperor Penguin Optimization Algorithm

## Description
This project implements the **Emperor Penguin Optimization (EPO)** algorithm for the **Traveling Salesman Problem (TSP)** using three parallelization strategies:
- **MPI**: process-level parallelism
- **OpenMP**: thread-level parallelism
- **Hybrid (MPI + OpenMP)**

It also includes a city generator and two Python scripts for plotting EPO execution logs.

## Repository Structure
```
├── generate_cities.c              # generates a cities.txt file with N random cities
├── plot_epo.py                    # EPO log plotting script
├── plot_epo2.py                   # EPO log plotting script
├── Parallel
│   ├── MPI
│   │   ├── cities.txt             # input city file
│   │   ├── epo_log.csv            # output log file
│   │   └── tsp_epo_mpi.c          # MPI C implementation
│   │   └── tsp_epo_mpi.sh         # test script
│   ├── OpenMP
│   │   ├── cities.txt             # input city file
│   │   ├── epo_log.csv            # output log file
│   │   ├── tsp_epo_mpi.c          # OpenMP C implementation
|   |   └── tsp_epo_open.sh        # test script
│   └── Hybrid
│       ├── cities.txt             # input city file
│       ├── epo_log.csv            # output log file
│       ├── tsp_epo_hybrid.c       # MPI + OpenMP C implementation
│       └── tsp_epo_hybrid.sh      # test script
└── README.md                      
```

## Installation and Usage

### 1. Compile and run OpenMP version
1. **Go to folder**  
   ```bash
   cd Parallel/MPI
   ```
2. **Load library**  
   ```bash
   module load mpich-3.2
   ```
3. **Compile**  
   ```bash
   mpicc -g -Wall -std=c99 -o tsp_epo_mpi tsp_epo_mpi.c -lm
   ```
4. **Run**  
   ```bash
   qsub tsp_epo_mpi.sh
   ```

### 2. Compile and run OpenMP version
1. **Go to folder**  
   ```bash
   cd Parallel/OpenMP
   ```
2. **Compile**  
   ```bash
   mpicc -g -Wall -std=c99 -fopenmp -o tsp_epo_open tsp_epo_open.c -lm
   ```
3. **Run**  
   ```bash
   qsub tsp_epo_open.sh
   ```

### 3. Compile and run Hybrid (MPI + OpenMP) version
1. **Go to folder**  
   ```bash
   cd Parallel/Hybrid
   ```
2. **Load library**  
   ```bash
   module load mpich-3.2
   ```
3. **Compile**  
   ```bash
   mpicc -g -Wall -std=c99 -fopenmp -o tsp_epo_hybrid tsp_epo_hybrid.c -lm
   ```
4. **Run**  
   ```bash
   qsub tsp_epo_hybrid.sh
   ```

## Additional Tools

### City Generator
```bash
gcc generate_cities.c -o generate_cities
./generate_cities
# Enter the number of cities when prompted (e.g., 100000)
```

### Plotting Results
```bash
# Interactive graph of fitness distribution per iteration
python3 plot_epo.py
# Graph of the best fitness per iteration
python3 plot_epo2.py
```
