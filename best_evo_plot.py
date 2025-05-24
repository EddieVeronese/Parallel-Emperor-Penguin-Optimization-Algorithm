import pandas as pd
import matplotlib.pyplot as plt

# Load CSV
df = pd.read_csv("epo_log.csv")

# Find the minimum value of 'length' for each 'iteration'
best_per_iteration = df.groupby("iteration")["length"].min().reset_index()

# Plot
plt.figure(figsize=(10, 6))
plt.plot(best_per_iteration["iteration"], best_per_iteration["length"], marker='o')
plt.title("Evolution of the length of the best path")
plt.xlabel("Iteration")
plt.ylabel("Better length")
plt.grid(True)
plt.tight_layout()
plt.show()
