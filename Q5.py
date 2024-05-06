#!/usr/bin/env python
# coding: utf-8

# In[18]:


import numpy as np
import time

# Function to compute DFT of a sequence x
def compute_dft(x):
    N = len(x)
    n = np.arange(N)
    k = n.reshape((N, 1))
    W = np.exp(-2j * np.pi * k * n / N)/np.sqrt(N)
    return np.dot(W, x)

# Number of elements
n = 2**4

# Generate random numbers
x = np.random.random(n)

# Compute DFT using direct computation
start_time = time.time()
dft_direct = compute_dft(x)
direct_time = time.time() - start_time

# Compute DFT using numpy's FFT
start_time = time.time()
dft_fft = np.fft.fft(x, norm="ortho")
fft_time = time.time() - start_time

print(f"Time taken for direct computation: {direct_time} seconds")
print(f"Time taken for numpy's FFT: {fft_time} seconds")
print(x)
print(dft_direct)
print(dft_fft)


# In[27]:


import numpy as np
import time
import matplotlib.pyplot as plt

# Function to compute DFT of a sequence x
def compute_dft(x):
    N = len(x)
    n = np.arange(N)
    k = n.reshape((N, 1))
    W = np.exp(-2j * np.pi * k * n / N)/np.sqrt(N)
    return np.dot(W, x)

# Initialize lists to store time taken by each method
direct_times = []
fft_times = []

# Range of n values
n_values = range(4, 1024)

for n in n_values:
    # Generate random numbers
    x = 10*np.random.random(n)
    
    # Compute DFT using direct computation
    start_time = time.time()
    dft_direct = compute_dft(x)
    direct_time = time.time() - start_time
    direct_times.append(direct_time)

    # Compute DFT using numpy's FFT
    start_time = time.time()
    dft_fft = np.fft.fft(x)/np.sqrt(n)
    fft_time = time.time() - start_time
    fft_times.append(fft_time)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(n_values, direct_times, label='Direct Computation')
plt.plot(n_values, fft_times, label='numpy FFT')
plt.xlabel('Number of Elements (n)')
plt.ylabel('Time (seconds)')
plt.title('Time taken for DFT computation vs. Number of Elements')
plt.legend()
plt.grid(True)
plt.show()


# In[28]:


import numpy as np
import time
import matplotlib.pyplot as plt

# Function to compute DFT of a sequence x
def compute_dft(x):
    N = len(x)
    n = np.arange(N)
    k = n.reshape((N, 1))
    W = np.exp(-2j * np.pi * k * n / N)/np.sqrt(N)
    return np.dot(W, x)

# Initialize lists to store time taken by each method
direct_times = []
fft_times = []

# Range of n values
n_values = range(4, 100)

for n in n_values:
    # Generate random numbers
    x = 10*np.random.random(n)
    
    # Compute DFT using direct computation
    start_time = time.time()
    dft_direct = compute_dft(x)
    direct_time = time.time() - start_time
    direct_times.append(direct_time)

    # Compute DFT using numpy's FFT
    start_time = time.time()
    dft_fft = np.fft.fft(x)/np.sqrt(n)
    fft_time = time.time() - start_time
    fft_times.append(fft_time)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(n_values, direct_times, label='Direct Computation')
plt.plot(n_values, fft_times, label='numpy FFT')
plt.xlabel('Number of Elements (n)')
plt.ylabel('Time (seconds)')
plt.title('Time taken for DFT computation vs. Number of Elements')
plt.legend()
plt.grid(True)
plt.show()


# In[29]:


import numpy as np
import time
import matplotlib.pyplot as plt

# Function to compute DFT of a sequence x
def compute_dft(x):
    N = len(x)
    n = np.arange(N)
    k = n.reshape((N, 1))
    W = np.exp(-2j * np.pi * k * n / N)/np.sqrt(N)
    return np.dot(W, x)

# Initialize lists to store time taken by each method
direct_times = []
fft_times = []

# Range of n values
n_values = range(1000, 1100)

for n in n_values:
    # Generate random numbers
    x = 10*np.random.random(n)
    
    # Compute DFT using direct computation
    start_time = time.time()
    dft_direct = compute_dft(x)
    direct_time = time.time() - start_time
    direct_times.append(direct_time)

    # Compute DFT using numpy's FFT
    start_time = time.time()
    dft_fft = np.fft.fft(x)/np.sqrt(n)
    fft_time = time.time() - start_time
    fft_times.append(fft_time)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(n_values, direct_times, label='Direct Computation')
plt.plot(n_values, fft_times, label='numpy FFT')
plt.xlabel('Number of Elements (n)')
plt.ylabel('Time (seconds)')
plt.title('Time taken for DFT computation vs. Number of Elements')
plt.legend()
plt.grid(True)
plt.show()





