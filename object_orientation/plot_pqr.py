import numpy as np
import matplotlib.pyplot as plt

filename = "output_files/f16_states.csv"

# Load data
data = np.genfromtxt(filename, delimiter=',', skip_header=1)

# Remove rows with NaNs
data = data[~np.isnan(data).any(axis=1)]

# Extract pqr
t = data[:, 0]
p = data[:, 4] * 180.0 / np.pi 
q = data[:, 5] * 180.0 / np.pi 
r = data[:, 6] * 180.0 / np.pi 

# Plot
plt.plot(t,p)
plt.plot(t,q)
plt.plot(t,r) 
plt.xlabel('Time [sec]')
plt.ylabel('Rotation Rates [deg]')
plt.title('F16 Rotation Rates')
plt.legend(['Roll Rate', 'Pitch Rate', 'Yaw Rate'])


plt.show()

