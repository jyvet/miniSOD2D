import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Read the data from the file
file_path = 'performance_results_cartesian.txt'
df = pd.read_csv(file_path, delim_whitespace=True)
df.loc[df['Speedup'] > 4, 'Speedup'] /= 4
df.loc[df['Speedup'] > 2, 'Speedup'] /= 2
df.loc[df['Speedup'] > 1.5, 'Speedup'] = 1
print(df)


# 1. 3D plot with Tile_X, Tile_Y, and Speedup
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
sc = ax.scatter(df['Tile_X'], df['Tile_Y'], df['Speedup'], c=df['Speedup'], cmap='viridis')
ax.set_xlabel('Tile_X')
ax.set_ylabel('Tile_Y')
ax.set_zlabel('Speedup')
ax.set_title('Speedup vs Tile Dimensions')
plt.colorbar(sc, label='Speedup')
plt.show()

# 2. 2D plot with speedup as y and tile dimension product as x
df['Tile_Product'] = df['Tile_X'] * df['Tile_Y']
groups = df.groupby('Tile_X')

plt.figure(figsize=(12, 6))
for name, group in groups:
    plt.scatter(group['Tile_Product'], group['Speedup'], label=f'Tile_X={name}')
plt.xlabel('Tile_X * Tile_Y')
plt.ylabel('Speedup')
plt.title('Speedup vs Tile Product')
plt.legend(title='Tile_X')
plt.grid(True)
plt.show()
