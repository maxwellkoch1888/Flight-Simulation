import numpy as np
import matplotlib.pyplot as plt

filename = "output_files/f16_states.csv"

# Load data
data = np.genfromtxt(filename, delimiter=',', skip_header=1)

# Remove rows with NaNs
data = data[~np.isnan(data).any(axis=1)]

# Extract pqr
t = data[:, 0]
da = data[:, 13]
de = data[:, 14]
dr = data[:, 15]
tau = data[:, 17]

plt.figure()

# Control surfaces
plt.subplot(2,1,1)
plt.plot(t, da)
plt.plot(t, de)
plt.plot(t, dr)
plt.ylabel('Deflection [deg]')
plt.title('F16 Control Surfaces')
plt.legend(['Aileron', 'Elevator', 'Rudder'])

# Throttle
plt.subplot(2,1,2)
plt.plot(t, tau)
plt.xlabel('Time [sec]')
plt.ylabel('Throttle')
plt.title('Throttle Command')

plt.tight_layout()
plt.show()

