import pandas as pd
import matplotlib.pyplot as plt

# Carica il CSV
df = pd.read_csv("epo_log.csv")

# Trova il valore minimo di 'lunghezza' per ogni 'iterazione'
best_per_iteration = df.groupby("iterazione")["lunghezza"].min().reset_index()

# Plot
plt.figure(figsize=(10, 6))
plt.plot(best_per_iteration["iterazione"], best_per_iteration["lunghezza"], marker='o')
plt.title("Evoluzione della lunghezza del percorso migliore")
plt.xlabel("Iterazione")
plt.ylabel("Lunghezza migliore")
plt.grid(True)
plt.tight_layout()
plt.show()
