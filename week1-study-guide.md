# Week 1 Study Guide: Python Basics, NumPy, and Visualization

## Overview

This study guide covers the essential Python concepts needed for Week 1 of ENGR 240. The activities, worksheets, and programming assignments for week 1 are designed to provide practice with these concepts and techniques.

## 1. Variables, Data Types, and Math Expressions

### Key Concepts

- **Variables**: Names that store values
- **Basic Data Types**:
  - `int`: Integer numbers (e.g., 1, 42, -10)
  - `float`: Floating-point numbers (e.g., 3.14, -2.5)
  - `str`: Text strings (e.g., "Hello", 'Engineering')
  - `bool`: Boolean values (True, False)
- **Math Operators**:
  - Addition: `+`
  - Subtraction: `-` 
  - Multiplication: `*`
  - Division: `/`
  - Integer division: `//`
  - Modulus (remainder): `%`
  - Exponentiation: `**`
- **Math Functions**:
  - From the `math` module: `sqrt()`, `sin()`, `cos()`, `exp()`, `log()`, etc.

### Example Code

```python
# Import necessary libraries
import numpy as np
import math

# Calculate factorials for MacLaurin series terms
fact3 = 3 * 2 * 1  # 3!
fact5 = 5 * 4 * 3 * 2 * 1  # 5!
fact7 = 7 * 6 * 5 * 4 * 3 * 2 * 1  # 7!

# Define x values for sine approximation
x = np.linspace(-2*np.pi, 2*np.pi, 1000)

# MacLaurin series approximations for sine
approx1 = x  # First term approximation
approx2 = x - (x**3) / fact3  # Second term approximation
approx3 = x - (x**3) / fact3 + (x**5) / fact5  # Third term approximation
```

### Resources

- **Berkeley Python Numerical Methods**:
  - https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter02.00-Variables-and-Basic-Data-Structures.html
 
- **Sullivan Numerical Methods YouTube Playlist**:
  - https://youtu.be/V5MPjDIw_sE?si=LVjZ-EndlCCY5kGt
  - https://youtu.be/yCHK9jrkUAU?si=aHS20aJMuShU_U-o 

- **Python 4 Everybody Videos (ENGR 151 Resource)**:
  - https://youtu.be/7KHdV6FSpo8?si=TnSQxAo1oN8A07-s
  - https://youtu.be/kefrGMAglGs?si=l_oD5O6IBhqdz6Wx

## 2. Overview of Libraries

### NumPy

NumPy is the fundamental package for scientific computing in Python. It provides:

- N-dimensional array objects
- Sophisticated mathematical functions
- Tools for linear algebra, random number generation, and Fourier transforms

### Matplotlib

Matplotlib is a plotting library that produces publication-quality figures. It's used for:

- Creating static, animated, and interactive visualizations
- Making basic plots (line plots, scatter plots, bar charts, histograms, etc.)
- Customizing plots with titles, labels, legends, and more

### Example Code

```python
# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt

# Define parameter values for contaminant decay model
A0 = 200  # Initial concentration of contaminant A (ppm)
B0 = 100  # Initial concentration of contaminant B (ppm)
kA = 1    # Decay rate of contaminant A (day^-1)
kB = 0.5  # Decay rate of contaminant B (day^-1)

# Time domain for simulation
t = np.linspace(0, 6, 100)  # 100 time points from 0 to 6 days

# Calculate contaminant concentration
p = A0 * np.exp(-kA * t) + B0 * np.exp(-kB * t)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(t, p, linewidth=2)
plt.title('Contaminant Decay Model')
plt.xlabel('Time (days)')
plt.ylabel('Concentration (ppm)')
plt.grid(True)
plt.show()
```

## 3. Basic NumPy Array Creation and Manipulation

### Creating Arrays

- `np.array()`: Create an array from a Python list
- `np.zeros()`: Create an array filled with zeros
- `np.ones()`: Create an array filled with ones
- `np.linspace()`: Create an array with evenly spaced values (inclusive)
- `np.arange()`: Create an array with evenly spaced values (exclusive of end)

### Array Operations

- Element-wise operations: `+`, `-`, `*`, `/`, `**`
- Functions: `np.sin()`, `np.cos()`, `np.exp()`, `np.log()`, etc.
- Array statistics: `np.mean()`, `np.max()`, `np.min()`, `np.sum()`
- Indexing and slicing: `array[0]`, `array[1:5]`, etc.

### Example Code

```python
# From the Contaminant Decay Sensitivity Analysis example

# Define decay rates for sensitivity analysis
kB_values = [0.5, 1.0, 2.0, 5.0]  # Different decay rates to compare
t = np.linspace(0, 6, 100)  # Time points from 0 to 6 days
standard = 10  # Regulatory standard (ppm)

# Calculate concentrations for different decay rates
A0 = 200  # Initial concentration of contaminant A
B0 = 100  # Initial concentration of contaminant B
kA = 1    # Decay rate of contaminant A

# Calculate half-life for each decay rate
half_lives = np.log(2) / np.array(kB_values)

# Display half-lives with formatting
for i, kB in enumerate(kB_values):
    half_life = half_lives[i]
    print(f"When kB = {kB}: Half-life = {half_life:.2f} days")
```

### Resources

- **Berkeley Python Numerical Methods**:
  - https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter02.07-Introducing_numpy_arrays.html

- **Sullivan Numerical Methods YouTube Playlist**:
  - https://youtu.be/0dPSwArrgsQ?si=TmMaQn_GO1V-lgck
  - https://youtu.be/4k7elXtKLUo?si=Qn7Jy0l8cP9eyNe0

## 4. Visualization with Matplotlib

### Basic Plotting

- `plt.figure()`: Create a new figure
- `plt.plot()`: Plot y versus x as lines and/or markers
- `plt.scatter()`: Create a scatter plot
- `plt.xlabel()`, `plt.ylabel()`: Add axis labels
- `plt.title()`: Add a title
- `plt.legend()`: Add a legend
- `plt.grid()`: Add a grid
- `plt.show()`: Display the figure

### Plot Customization

- Line styles: `'-'`, `'--'`, `'-.'`, `':'`
- Markers: `'o'`, `'s'`, `'^'`, `'*'`
- Colors: `'r'`, `'g'`, `'b'`, `'k'`, or hex values like `'#FF5733'`
- Figure size: `plt.figure(figsize=(width, height))`

### Example Code

```python
# From the decaying oscillation worksheet

# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt

# Define parameter values for decaying oscillation
A = 10          # Initial amplitude
k = 0.5         # Decay rate
omega = np.pi/2 # Angular frequency (rad/s)
phi = -np.pi/4  # Phase shift (radians)

# Set up time vector
T = 2 * np.pi / omega  # Period of oscillation
t_end = 3 * T          # End time (3 periods)
t = np.linspace(0, t_end, 200)  # Time vector with 200 points

# Generate the data for the decaying cosine model
y = A * np.exp(-k * t) * np.cos(omega * t + phi)

# Create and format the plot
plt.figure(figsize=(10, 6))
plt.plot(t, y, 'b-', linewidth=2)
plt.title('Decaying Oscillation Example')
plt.xlabel('t (seconds)')
plt.ylabel('y')
plt.grid(True)

# Add a text box with the equation for reference
equation = f'$y = {A}e^{{-{k}t}}\\cos({omega:.2f}t {phi:+.2f})$'
plt.text(0.6*t_end, 0.8*A, equation, fontsize=12, 
         bbox=dict(facecolor='white', alpha=0.5))

plt.show()
```

### Resources

- **Berkeley Python Numerical Methods**:
  - https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter12.00-Visualization-and-Plotting.html

- **Sullivan Numerical Methods YouTube Playlist**:
  - https://youtu.be/-e5e4AZFwjg?si=GNwOydeuKZ1Kg1pp
  - https://youtu.be/v2KMuXmIBZg?si=O6YZbY_LQfxaIkuN
  - https://youtu.be/LABgO9lIIuo?si=2FgJW3wXPjxiHpJ3

## Tips for Week 1 Assignments

1. **Start with basic calculations**: Make sure you understand how to use Python for basic math operations before moving to more complex tasks.

2. **Break down the problems**: Each programming assignment can be broken into smaller steps - tackle them one at a time.

3. **Use vectorized operations**: Use NumPy's vectorized operations on arrays for better performance and cleaner code.

4. **Focus on visualization quality**: Pay attention to proper labeling, titles, and legend entries in your plots.

5. **Comment your code**: Document your code with plenty of comments explaining your approach and what each section of code does.

6. **Test your code**: Test your code at intermediate steps to verify it works as expected.
