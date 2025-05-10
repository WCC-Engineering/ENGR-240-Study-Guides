# Week 6 Study Guide: Interpolation and Splines

## Overview

This study guide covers interpolation techniques for constructing new data points within the range of a discrete set of known data points. These methods are essential in engineering for approximating functions, estimating intermediate values, and creating smooth curves through data points. We'll explore various interpolation approaches including linear interpolation, polynomial interpolation, Lagrange polynomials, and spline interpolation.

## 1. Introduction to Interpolation

### Key Concepts

- **Definition**: Interpolation is the process of finding a function that passes exactly through a given set of data points
- **Purpose**:
  - Estimate values between measured data points
  - Create smooth curves connecting discrete data
  - Reconstruct missing data values
  - Approximate complex functions with simpler representations
- **Common Engineering Applications**:
  - Digital signal processing
  - Computer graphics and animation
  - Analysis of experimental data
  - Finite element analysis
  - Control systems design
  - Numerical integration and differentiation

### Mathematical Basis

For a set of data points $(x_i, y_i)$ where $i = 0, 1, 2, ..., n$, an interpolating function $f(x)$ satisfies:

$$f(x_i) = y_i \text{ for all } i = 0, 1, 2, ..., n$$

This means the function passes exactly through every data point.

### Key Properties of Interpolation

1. **Accuracy at Data Points**: By definition, interpolating functions pass exactly through the given data points
2. **Behavior Between Points**: Different interpolation methods produce different behaviors between data points
3. **Error Characteristics**: The error in interpolation depends on:
   - The interpolation method used
   - The spacing of data points
   - The smoothness of the underlying function
4. **Extrapolation Limitations**: Interpolation functions often perform poorly when extrapolating beyond the data range

### Types of Interpolation Methods

1. **Linear Interpolation**: Connects data points with straight lines
2. **Polynomial Interpolation**: Uses a single polynomial to pass through all data points
3. **Lagrange Interpolation**: A specific form of polynomial interpolation
4. **Piecewise Interpolation**: Uses different functions for different intervals between data points
5. **Spline Interpolation**: Uses piecewise polynomials with smoothness constraints at the knots

### Resources

- **Davishahl Numerical Methods Videos**:
  - [Introduction to Interpolation] <!-- Video link to be added -->

- **Python Numerical Methods**:
  - [Interpolation Basics](https://pythonnumericalmethods.berkeley.edu/notebooks/chapter17.01-Interpolation.html)

## 2. Linear Interpolation

### Key Concepts

- **Definition**: Linear interpolation connects data points with straight line segments
- **Applications**:
  - Simple estimation of intermediate values
  - First approximation for many engineering problems
  - Computer graphics (e.g., color gradients)
  - Real-time calculations where speed is critical
- **Advantages**:
  - Simplicity and computational efficiency
  - Stability (no oscillations)
  - Easy to understand and implement
- **Limitations**:
  - Limited accuracy for curved data
  - Discontinuous first derivative at data points
  - Not suitable for smooth function approximation

### Mathematical Formulation

For two known points $(x_0, y_0)$ and $(x_1, y_1)$, the linear interpolant for a value $x$ where $x_0 \leq x \leq x_1$ is:

$$y = y_0 + \frac{y_1 - y_0}{x_1 - x_0}(x - x_0)$$

This can also be expressed as a weighted average:

$$y = \frac{x_1 - x}{x_1 - x_0}y_0 + \frac{x - x_0}{x_1 - x_0}y_1$$

### Implementation in Python

{% raw %}
```python
import numpy as np
import matplotlib.pyplot as plt

# Example of manual linear interpolation
def linear_interp(x, x0, y0, x1, y1):
    """
    Linear interpolation between two points (x0,y0) and (x1,y1).
    
    Parameters:
    -----------
    x : float
        Point at which to interpolate
    x0, y0 : float
        First data point
    x1, y1 : float
        Second data point
        
    Returns:
    --------
    y : float
        Interpolated value at x
    """
    # Check that x is within the range [x0, x1]
    if not (x0 <= x <= x1):
        raise ValueError(f"x={x} is outside the range [{x0}, {x1}]")
    
    # Calculate interpolation
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)

# Sample data points
x_data = np.array([0, 2, 4, 6, 8, 10])
y_data = np.array([1, 3, 8, 15, 20, 22])

# Create a dense x-grid for plotting the interpolation
x_interp = np.linspace(min(x_data), max(x_data), 100)
y_interp = np.zeros_like(x_interp)

# Perform linear interpolation
for i, x in enumerate(x_interp):
    # Find the appropriate segment
    for j in range(len(x_data)-1):
        if x_data[j] <= x <= x_data[j+1]:
            y_interp[i] = linear_interp(x, x_data[j], y_data[j], 
                                        x_data[j+1], y_data[j+1])
            break

# Using NumPy's built-in interp function
y_interp_numpy = np.interp(x_interp, x_data, y_data)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(x_data, y_data, 'o', label='Data points')
plt.plot(x_interp, y_interp, '-', label='Manual linear interpolation')
plt.plot(x_interp, y_interp_numpy, '--', label='NumPy interp')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Linear Interpolation Example')
plt.legend()
plt.grid(True)
plt.show()
```
{% endraw %}

### Applications and Considerations

1. **When to Use Linear Interpolation**:
   - Quick estimates with minimal computation
   - When data appears to change linearly between points
   - When higher-order methods might introduce unwanted oscillations
   - As a first approximation before trying more complex methods

2. **Error Analysis**:
   - If the true function is continuously differentiable, the error is proportional to the square of the distance between data points
   - Maximum error occurs where the second derivative of the true function has its maximum

3. **Python Implementation Options**:
   - `numpy.interp()`: Efficient 1D linear interpolation
   - `scipy.interpolate.interp1d`: More flexible with additional features
   - Manual implementation: For educational purposes or custom requirements

### Resources

- **Davishahl Numerical Methods Videos**:
  - [Linear Interpolation] <!-- Video link to be added -->

- **NumPy Documentation**:
  - [numpy.interp](https://numpy.org/doc/stable/reference/generated/numpy.interp.html)
