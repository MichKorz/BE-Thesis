import numpy as np

# Wczytanie wszystkich kolumn z pliku
data = np.loadtxt("Speck_L2_16k-hz.dat")

# Obliczenie amplitudy z dwóch ostatnich kolumn
amplitude = np.hypot(data[:, -2], data[:, -1])

# Zastąpienie części rzeczywistej amplitudą oraz wyzerowanie części urojonej
data[:, -2] = amplitude
data[:, -1] = 0.0

# Zapis do nowego pliku
np.savetxt("Speck_L2_16k-hz-tl.dat", data, fmt="%.7g")
