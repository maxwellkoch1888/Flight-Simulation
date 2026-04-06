import numpy as np
import matplotlib.pyplot as plt

filename = "output_files/f16_states.csv"

# Load data
data = np.genfromtxt(filename, delimiter=',', skip_header=1)

# Remove rows with NaNs
data = data[~np.isnan(data).any(axis=1)]

# Extract pqr
t = data[:, 0]
z = data[:, 9] * -1.0

# Plot
plt.plot(t,z) 
plt.xlabel('Time [sec]')
plt.ylabel('Altitude [ft]')
plt.title('F16 Altitude')


plt.show()

