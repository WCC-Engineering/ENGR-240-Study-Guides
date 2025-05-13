# Using numpy.trapz
x = np.linspace(a, b, 100)
y = f(x)
trapz_result = np.trapz(y, x)
print(f"numpy.trapz: {trapz_result}, Error: {abs(trapz_result-exact)}")

# Using scipy.integrate.trapz (same as numpy.trapz)
trapz_scipy = integrate.trapz(y, x)
print(f"scipy.integrate.trapz: {trapz_scipy}, Error: {abs(trapz_scipy-exact)}")

# Using scipy.integrate.simpson (Simpson's rule)
simpson_result = integrate.simpson(y, x)
print(f"scipy.integrate.simpson: {simpson_result}, Error: {abs(simpson_result-exact)}")

# Using scipy.integrate.quad (adaptive quadrature)
quad_result, quad_error = integrate.quad(f, a, b)
print(f"scipy.integrate.quad: {quad_result}, Error estimate: {quad_error}")
```
{% endraw %}

### Error Behavior in Newton-Cotes Methods

Each Newton-Cotes formula has a characteristic error behavior:

1. **Trapezoid Rule**: Error ∝ h², where h is the step size
2. **Simpson's 1/3 Rule**: Error ∝ h⁴
3. **Simpson's 3/8 Rule**: Error ∝ h⁴

For a desired accuracy, Simpson's rules typically require far fewer function evaluations than the trapezoid rule.

### Resources

- **Davishahl Numerical Methods Videos**
  - [Trapezoid Rule](#) <!-- Video link to be inserted by you -->
  - [Simpson's Rule](#) <!-- Video link to be inserted by you -->

- **Berkeley Numerical Methods**
  - [Trapezoidal Rule](https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter19.01-Numerical-Integration-Trapezoidal-Rule.html)
  - [Simpson's Rule](https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter19.02-Numerical-Integration-Simpsons-Rule.html)

## 7. Gauss-Legendre Quadrature

### Key Concepts

- **Gaussian Quadrature**: Integration method that optimally selects both evaluation points and weights
- **Gauss-Legendre**: Uses roots of Legendre polynomials as evaluation points
- **Advantages**:
  - Higher accuracy with fewer function evaluations
  - Exact for polynomials of degree up to 2n-1 using n points
- **Limitations**:
  - More complex to implement than Newton-Cotes
  - Requires function evaluations at specific, non-uniform points

### Mathematical Form

The Gauss-Legendre quadrature approximates the integral as:

$$\int_{-1}^{1} f(x) \, dx \approx \sum_{i=1}^{n} w_i f(x_i)$$

Where:
- $x_i$ are the roots of the nth-degree Legendre polynomial
- $w_i$ are the corresponding weights

For an interval [a, b] different from [-1, 1], the formula becomes:

$$\int_{a}^{b} f(x) \, dx \approx \frac{b-a}{2} \sum_{i=1}^{n} w_i f\left(\frac{b-a}{2}x_i + \frac{a+b}{2}\right)$$

### Implementation Using SciPy

SciPy provides Gauss-Legendre quadrature through `scipy.special.legendre` and `scipy.special.roots_legendre`:

{% raw %}
```python
import numpy as np
from scipy import special, integrate
import matplotlib.pyplot as plt

def f(x):
    return np.exp(-x**2)

# Integration limits
a, b = 0, 1

# Exact value (using scipy.integrate.quad)
exact = integrate.quad(f, a, b)[0]

def gauss_legendre(f, a, b, n):
    """
    Gauss-Legendre quadrature with n points.
    """
    # Get Legendre polynomial roots and weights
    x_i, w_i = special.roots_legendre(n)
    
    # Transform from [-1, 1] to [a, b]
    x_transformed = 0.5 * (b - a) * x_i + 0.5 * (a + b)
    
    # Evaluate function at transformed points
    f_values = f(x_transformed)
    
    # Apply quadrature formula
    return 0.5 * (b - a) * np.sum(w_i * f_values)

# Compare with different numbers of points
n_points = [2, 3, 4, 5, 6, 10]
gl_results = []
gl_errors = []

for n in n_points:
    result = gauss_legendre(f, a, b, n)
    error = abs(result - exact)
    gl_results.append(result)
    gl_errors.append(error)
    print(f"Gauss-Legendre with {n} points: {result:.10f}, Error: {error:.10e}")

# Compare with trapezoidal rule (same number of points)
trap_errors = []
simp_errors = []

for n in n_points:
    # For trapezoidal rule with n points
    x = np.linspace(a, b, n)
    y = f(x)
    trap_result = np.trapz(y, x)
    trap_error = abs(trap_result - exact)
    trap_errors.append(trap_error)
    
    # For Simpson's rule (when applicable)
    if n >= 3 and (n - 1) % 2 == 0:  # Need odd number of points for Simpson's rule
        simp_result = integrate.simpson(f(np.linspace(a, b, n)), dx=(b-a)/(n-1))
        simp_errors.append(abs(simp_result - exact))
    else:
        simp_errors.append(np.nan)  # Not applicable

# Plot comparison
plt.figure(figsize=(10, 6))
plt.semilogy(n_points, gl_errors, 'o-', label='Gauss-Legendre')
plt.semilogy(n_points, trap_errors, 's-', label='Trapezoidal')
plt.semilogy([n for n, e in zip(n_points, simp_errors) if not np.isnan(e)], 
             [e for e in simp_errors if not np.isnan(e)], '^-', label='Simpson')
plt.xlabel('Number of Points')
plt.ylabel('Absolute Error (log scale)')
plt.legend()
plt.grid(True)
plt.title('Error Comparison Between Integration Methods')
plt.show()
```
{% endraw %}

### Using SciPy's Fixed-Point Quadrature

SciPy also provides direct implementation of Gaussian quadrature methods:

{% raw %}
```python
from scipy import integrate

def f(x):
    return np.exp(-x**2)

a, b = 0, 1
exact = integrate.quad(f, a, b)[0]

# Using fixed_quad (Gauss-Legendre quadrature)
for n in [5, 10, 15]:
    result, error = integrate.fixed_quad(f, a, b, n=n)
    print(f"fixed_quad with {n} points: {result:.10f}, Error: {abs(result-exact):.10e}")
```
{% endraw %}

### Comparing Efficiency with Newton-Cotes

For a smooth function, Gaussian quadrature typically achieves higher accuracy with fewer function evaluations:

{% raw %}
```python
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import time

def f(x):
    return np.sin(x) * np.exp(-0.1*x)

a, b = 0, 10
exact = integrate.quad(f, a, b)[0]

# Function to test different methods with various numbers of points
def compare_methods(max_points=30):
    n_points = list(range(2, max_points+1))
    
    # Storage for results
    gl_errors = []
    trap_errors = []
    simp_errors = []
    gl_times = []
    trap_times = []
    simp_times = []
    
    for n in n_points:
        # Gauss-Legendre
        start = time.time()
        gl_result, _ = integrate.fixed_quad(f, a, b, n=n)
        gl_time = time.time() - start
        gl_error = abs(gl_result - exact)
        gl_errors.append(gl_error)
        gl_times.append(gl_time)
        
        # Trapezoidal
        start = time.time()
        x = np.linspace(a, b, n)
        y = f(x)
        trap_result = np.trapz(y, x)
        trap_time = time.time() - start
        trap_error = abs(trap_result - exact)
        trap_errors.append(trap_error)
        trap_times.append(trap_time)
        
        # Simpson's (when applicable)
        if n >= 3:
            start = time.time()
            simp_result = integrate.simpson(y, x)
            simp_time = time.time() - start
            simp_error = abs(simp_result - exact)
            simp_errors.append(simp_error)
            simp_times.append(simp_time)
        else:
            simp_errors.append(np.nan)
            simp_times.append(np.nan)
    
    return n_points, gl_errors, trap_errors, simp_errors, gl_times, trap_times, simp_times

# Run comparison
n_points, gl_errors, trap_errors, simp_errors, gl_times, trap_times, simp_times = compare_methods(30)

# Plot error vs. number of points
plt.figure(figsize=(12, 10))
plt.subplot(2, 1, 1)
plt.semilogy(n_points, gl_errors, 'o-', label='Gauss-Legendre')
plt.semilogy(n_points, trap_errors, 's-', label='Trapezoidal')
plt.semilogy([n for n, e in zip(n_points, simp_errors) if not np.isnan(e)], 
             [e for e in simp_errors if not np.isnan(e)], '^-', label='Simpson')
plt.xlabel('Number of Points')
plt.ylabel('Absolute Error (log scale)')
plt.legend()
plt.grid(True)
plt.title('Error Comparison Between Integration Methods')

# Plot computational efficiency (error vs. time)
plt.subplot(2, 1, 2)
plt.loglog(gl_times, gl_errors, 'o-', label='Gauss-Legendre')
plt.loglog(trap_times, trap_errors, 's-', label='Trapezoidal')
plt.loglog([t for t, e in zip(simp_times, simp_errors) if not np.isnan(e)], 
           [e for e in simp_errors if not np.isnan(e)], '^-', label='Simpson')
plt.xlabel('Computation Time (s)')
plt.ylabel('Absolute Error (log scale)')
plt.legend()
plt.grid(True)
plt.title('Error vs. Computation Time')

plt.tight_layout()
plt.show()
```
{% endraw %}

### Resources

- **Davishahl Numerical Methods Videos**
  - [Gauss-Legendre Quadrature](#) <!-- Video link to be inserted by you -->

- **Berkeley Numerical Methods**
  - [Gaussian Quadrature](https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter19.03-Numerical-Integration-Gaussian-Quadrature.html)

## 8. Using NumPy and SciPy for Differentiation and Integration

### Key Functions in NumPy for Differentiation

NumPy provides gradient calculation through `numpy.gradient()`:

{% raw %}
```python
import numpy as np
import matplotlib.pyplot as plt

# Example with evenly spaced data
x = np.linspace(0, 2*np.pi, 100)
y = np.sin(x)

# Calculate gradient 
dy_dx = np.gradient(y, x)  # First argument is values, second is coordinates

# Compare with analytical derivative
plt.figure(figsize=(10, 6))
plt.plot(x, y, label='sin(x)')
plt.plot(x, dy_dx, label='np.gradient(sin(x))')
plt.plot(x, np.cos(x), '--', label='cos(x) (Analytical)')
plt.legend()
plt.title('NumPy Gradient vs. Analytical Derivative')
plt.grid(True)
plt.show()

# Example with unevenly spaced data
x_uneven = np.sort(np.random.uniform(0, 2*np.pi, 100))
y_uneven = np.sin(x_uneven)

# Calculate gradient for uneven spacing
dy_dx_uneven = np.gradient(y_uneven, x_uneven)

# Compare with analytical derivative
plt.figure(figsize=(10, 6))
plt.plot(x_uneven, y_uneven, 'o', markersize=3, label='sin(x) (Uneven)')
plt.plot(x_uneven, dy_dx_uneven, label='np.gradient(sin(x))')
plt.plot(x_uneven, np.cos(x_uneven), '--', label='cos(x) (Analytical)')
plt.legend()
plt.title('NumPy Gradient with Uneven Spacing')
plt.grid(True)
plt.show()

# For 2D data (partial derivatives)
x = np.linspace(-5, 5, 50)
y = np.linspace(-5, 5, 50)
X, Y = np.meshgrid(x, y)
Z = np.sin(X) * np.cos(Y)

# Calculate gradients in both directions
dZ_dx, dZ_dy = np.gradient(Z, x, y)

# Plot the function and its partial derivatives
fig = plt.figure(figsize=(15, 5))

ax1 = fig.add_subplot(131)
c1 = ax1.contourf(X, Y, Z, cmap='viridis')
plt.colorbar(c1, ax=ax1)
ax1.set_title('f(x,y) = sin(x)cos(y)')

ax2 = fig.add_subplot(132)
c2 = ax2.contourf(X, Y, dZ_dx, cmap='viridis')
plt.colorbar(c2, ax=ax2)
ax2.set_title('∂f/∂x')

ax3 = fig.add_subplot(133)
c3 = ax3.contourf(X, Y, dZ_dy, cmap='viridis')
plt.colorbar(c3, ax=ax3)
ax3.set_title('∂f/∂y')

plt.tight_layout()
plt.show()
```
{% endraw %}

### Key Functions in SciPy for Integration

SciPy's `scipy.integrate` module offers various integration methods:

{% raw %}
```python
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

def f(x):
    return np.exp(-x**2) * np.sin(2*x)

# Define integration range
a, b = 0, 3

# 1. General-purpose integration with adaptive algorithm
result, error = integrate.quad(f, a, b)
print(f"quad result: {result:.10f}, error estimate: {error:.10e}")

# 2. Integration of fixed samples using various methods
x = np.linspace(a, b, 100)
y = f(x)

trapz_result = integrate.trapz(y, x)
simpson_result = integrate.simpson(y, x)
romb_result = integrate.romb(f(np.linspace(a, b, 2**4+1)), dx=(b-a)/(2**4))

print(f"trapz result: {trapz_result:.10f}")
print(f"simpson result: {simpson_result:.10f}")
print(f"romb result: {romb_result:.10f}")

# 3. Cumulative integration
cumtrapz_result = integrate.cumtrapz(y, x, initial=0)

plt.figure(figsize=(10, 6))
plt.plot(x, y, label='f(x)')
plt.plot(x, cumtrapz_result, label='Cumulative Integral')
plt.axhline(result, color='r', linestyle='--', label=f'quad result = {result:.4f}')
plt.legend()
plt.title('Function and its Cumulative Integral')
plt.grid(True)
plt.show()

# 4. Double (and multiple) integration
def g(x, y):
    return np.sin(x + y)

def bounds_y(x):
    return [0, 1]  # y varies from 0 to 1 for all x

def bounds_x():
    return [0, 1]  # x varies from 0 to 1

# Compute double integral
result_dblquad = integrate.dblquad(g, *bounds_x(), bounds_y)
print(f"Double integral result: {result_dblquad[0]:.10f}, error estimate: {result_dblquad[1]:.10e}")

# 5. Romberg integration with different orders
for k in range(1, 6):
    n = 2**k + 1
    result_romb = integrate.romb(f(np.linspace(a, b, n)), dx=(b-a)/(n-1))
    print(f"Romberg integration with 2^{k}+1 points: {result_romb:.10f}")
```
{% endraw %}

### SciPy's Integration with Error Control

SciPy's `quad` function provides sophisticated error control:

{% raw %}
```python
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

def f(x):
    return np.sin(100*x) / (x + 0.1)

# Basic usage
result, error = integrate.quad(f, 0, 1)
print(f"Default settings: {result:.10f}, error estimate: {error:.10e}")

# With explicit error control
result_tight, error_tight = integrate.quad(f, 0, 1, epsabs=1e-12, epsrel=1e-12)
print(f"Tight error control: {result_tight:.10f}, error estimate: {error_tight:.10e}")

# With increased maximum evaluations
result_maxeval, error_maxeval = integrate.quad(f, 0, 1, limit=200)
print(f"Increased max evaluations: {result_maxeval:.10f}, error estimate: {error_maxeval:.10e}")

# Compare with fixed-point quadrature
for n in [10, 50, 100, 200]:
    result_fixed = integrate.fixed_quad(f, 0, 1, n=n)[0]
    print(f"fixed_quad with {n} points: {result_fixed:.10f}")
```
{% endraw %}

### Resources

- **SciPy Documentation**
  - [numpy.gradient](https://numpy.org/doc/stable/reference/generated/numpy.gradient.html)
  - [scipy.integrate](https://docs.scipy.org/doc/scipy/reference/integrate.html)

## Tips for Week 7 Assignments

1. **For numerical differentiation**:
   - Always test methods on functions with known analytical derivatives
   - Carefully choose step sizes to balance truncation and round-off errors
   - For noisy data, consider spline differentiation over direct finite differences
   - Use numpy.gradient() for simple cases with tabular data

2. **For numerical integration**:
   - Start with scipy.integrate.quad for most single-variable integrals
   - For tabular data, use scipy.integrate.trapz or scipy.integrate.simpson
   - Compare different methods to understand their efficiency and accuracy characteristics
   - For highly oscillatory or singular functions, use adaptive methods

3. **Error analysis**:
   - Verify numerical results against known analytical solutions when possible
   - Observe how errors scale with step size or number of points
   - Remember that differentiation typically amplifies errors while integration smooths them

4. **Performance considerations**:
   - For repeated integrations with different parameters, consider precomputing Gaussian quadrature points and weights
   - For large datasets, vectorized operations with NumPy are much faster than Python loops

5. **Documentation best practices**:
   - Clearly state the methods used and their expected order of accuracy
   - Document assumptions about the input data
   - Include error estimates when presenting numerical results