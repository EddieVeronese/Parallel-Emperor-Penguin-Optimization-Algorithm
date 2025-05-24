import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# Load data
df = pd.read_csv("epo_log.csv")
df.columns = df.columns.str.strip()  

# Extremes
max_iter = df["iteration"].max()
max_penguin = df["penguin"].max()

# Figure and axes
fig, ax = plt.subplots(figsize=(10, 6))
plt.subplots_adjust(bottom=0.2)

# First subset
initial_iter = df["iteration"].min()
subset = df[df["iteration"] == initial_iter]

# Colors: all blue, one red
colors = ["red" if val == subset["length"].min() else "blue" for val in subset["length"]]
sc = ax.scatter(subset["penguin"], subset["length"], c=colors, edgecolor="black")

ax.set_xlim(0, max_penguin + 1)
ax.set_ylim(df["length"].min() * 0.98, df["length"].max() * 1.02)
ax.set_xlabel("Penguin")
ax.set_ylabel("Fitness (path length)")
title = ax.set_title(f"Fitness Distribution - Iteration {initial_iter}")
ax.grid(True)

# Slider
ax_slider = plt.axes([0.2, 0.05, 0.6, 0.04])
slider = Slider(ax_slider, "Iteration", valmin=df["iteration"].min(), valmax=max_iter, valinit=initial_iter, valstep=1)

# Update chart
def update(val):
    iter_selected = int(slider.val)
    subset = df[df["iteration"] == iter_selected]

    # Updated colors
    colors = ["red" if val == subset["length"].min() else "blue" for val in subset["length"]]
    sc.set_offsets(list(zip(subset["penguin"], subset["length"])))
    sc.set_color(colors)
    title.set_text(f"Fitness Distribution - Iteration {iter_selected}")
    fig.canvas.draw_idle()

slider.on_changed(update)

plt.show()
