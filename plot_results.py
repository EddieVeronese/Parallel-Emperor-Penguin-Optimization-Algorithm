import pandas as pd
import matplotlib.pyplot as plt

# read the csv file
df = pd.read_csv('optimization_results.csv', header=None, names=['Iteration', 'Best Path Length'])

print(df.head()) 

# create the graph 
plt.plot(df['Iteration'], df['Best Path Length'], marker='o', linestyle='-', color='b')
plt.title('Path improvement in EPO optimization')
plt.xlabel('Iteration')
plt.ylabel('Length of the best route')
plt.grid(True)
plt.show()
