import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# Carica i dati
df = pd.read_csv("epo_log.csv")
df.columns = df.columns.str.strip()  # Rimuove spazi nei nomi colonne

# Estremi
max_iter = df["iterazione"].max()
max_penguin = df["pinguino"].max()

# Figura e assi
fig, ax = plt.subplots(figsize=(10, 6))
plt.subplots_adjust(bottom=0.2)

# Primo subset
initial_iter = df["iterazione"].min()
subset = df[df["iterazione"] == initial_iter]

# Colori: tutti blu, uno rosso
colors = ["red" if val == subset["lunghezza"].min() else "blue" for val in subset["lunghezza"]]
sc = ax.scatter(subset["pinguino"], subset["lunghezza"], c=colors, edgecolor="black")

ax.set_xlim(0, max_penguin + 1)
ax.set_ylim(df["lunghezza"].min() * 0.98, df["lunghezza"].max() * 1.02)
ax.set_xlabel("Pinguino")
ax.set_ylabel("Fitness (lunghezza percorso)")
title = ax.set_title(f"Distribuzione fitness - Iterazione {initial_iter}")
ax.grid(True)

# Slider
ax_slider = plt.axes([0.2, 0.05, 0.6, 0.04])
slider = Slider(ax_slider, "Iterazione", valmin=df["iterazione"].min(), valmax=max_iter, valinit=initial_iter, valstep=1)

# Aggiorna grafico
def update(val):
    iter_selected = int(slider.val)
    subset = df[df["iterazione"] == iter_selected]

    # Colori aggiornati
    colors = ["red" if val == subset["lunghezza"].min() else "blue" for val in subset["lunghezza"]]
    sc.set_offsets(list(zip(subset["pinguino"], subset["lunghezza"])))
    sc.set_color(colors)
    title.set_text(f"Distribuzione fitness - Iterazione {iter_selected}")
    fig.canvas.draw_idle()

slider.on_changed(update)

plt.show()
