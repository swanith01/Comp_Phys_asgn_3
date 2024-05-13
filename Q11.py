import numpy as np
import matplotlib.pyplot as plt

# Define the box function
def box_function(x):
    return np.where(np.logical_and(-1 < x, x < 1), 1, 0)

# Generate x values
x = np.linspace(-3, 3, 1000)

# Compute the box function for each x
f = box_function(x)

# Compute the convolution of the box function with itself
convolution = np.convolve(f, f, mode='same') / len(x)  # Normalize the convolution

# Plotting
plt.figure(figsize=(10, 6))

# Plot the original box function
plt.plot(x, f, label='Box Function')

# Plot the convolution result
plt.plot(x, convolution, label='Convolution')

plt.title('Convolution of Box Function with Itself')
plt.xlabel('x')
plt.ylabel('Amplitude')
plt.legend()
plt.grid(True)
plt.show()

