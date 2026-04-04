import numpy as np
import matplotlib.pyplot as plt

filename = "output_files/f16_states.csv"

# Load data safely (handles headers + missing values)
data = np.genfromtxt(filename, delimiter=',', skip_header=1)

# Remove rows with NaNs
data = data[~np.isnan(data).any(axis=1)]

# Extract x, y, z (columns 8–10 → indices 7,8,9)
x = data[:, 7]
y = data[:, 8]
z = data[:, 9] * -1.0

# Plot
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.plot(x, y, z)

# Add start point (green) and end point (red)
ax.scatter(x[0], y[0], z[0], color='green', s=50, label='Start')
ax.scatter(x[-1], y[-1], z[-1], color='red', s=50, label='End')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('F16 Trajectory')

ax.set_box_aspect([1,1,1])


plt.show()

