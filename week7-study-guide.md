## Overview

This study guide covers numerical methods for differentiation and integration, essential techniques in engineering computation. Week 7 introduces finite difference methods for approximating derivatives, discusses truncation error and accuracy, explores differentiation of splines, and examines various numerical integration techniques including Newton-Cotes formulas and Gaussian quadrature.

## 1. Introduction to Numerical Differentiation

### Key Concepts

- **Differentiation**: The process of finding the rate of change of a function with respect to a variable
- **Analytical vs. Numerical Differentiation**: 
  - Analytical: Exact mathematical expressions for derivatives
  - Numerical: Approximations based on function values at discrete points
- **Applications**:
  - Solving differential equations
  - Analyzing rate of change in experimental data
  - Optimization problems
  - Sensitivity analysis in engineering systems

### Basic Definition of Differentiation

The derivative of a function f(x) at a point x is defined as:

{% raw %}
$$f'(x) = \lim_{h \to 0} \frac{f(x+h) - f(x)}{h}$$
{% endraw %}

This represents the instantaneous rate of change of f(x) with respect to x. Geometrically, it is the slope of the tangent line at point x on the curve of f(x).

In higher dimensions, partial derivatives measure the rate of change with respect to one variable while others are held constant. For a function f(x,y), the partial derivatives are:

{% raw %}
$$\frac{\partial f}{\partial x} = \lim_{h \to 0} \frac{f(x+h,y) - f(x,y)}{h}$$
{% endraw %}

{% raw %}
$$\frac{\partial f}{\partial y} = \lim_{h \to 0} \frac{f(x,y+h) - f(x,y)}{h}$$
{% endraw %}

The gradient of a multivariable function combines these partial derivatives into a vector that points in the direction of maximum increase.

### Challenges in Numerical Differentiation

Numerical differentiation is inherently ill-conditioned compared to numerical integration:

1. **Error Magnification**: Small errors in function values can lead to large errors in derivative approximations. This is because differentiation amplifies high-frequency components (including noise).

2. **Step Size Dilemma**: 
   - Too large: Approximation (truncation) error increases
   - Too small: Round-off error dominates due to subtraction of nearly equal numbers

3. **Noise Sensitivity**: For experimental data with noise, derivatives amplify the noise, often requiring smoothing techniques.

The fundamental mathematical issue can be seen by examining the condition number of differentiation. For a perturbation ε in the function value, the resulting perturbation in the derivative can be approximately {% raw %}$\frac{\varepsilon}{h}${% endraw %}, which grows as h decreases.

### Resources

- **Davishahl Numerical Methods Videos**
  - [Numerical Differentiation](https://youtu.be/XYjWnDJln-Y?si=yluR_buQAzJy7141)
  - [Total Numerical Error](https://youtu.be/cttS8j6VUJk?si=gdQRBp4Adu8u9Z0v)

- **Berkeley Numerical Methods**
  - [Numerical Differentiation](https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter20.00-Numerical-Differentiation.html)

## 2. Finite Difference Formulas

### Key Concepts

- **Finite Difference**: Approximation of derivatives using function values at discrete points
- **Categories**:
  - Forward difference: Uses points ahead of the evaluation point
  - Backward difference: Uses points behind the evaluation point
  - Central difference: Uses points on both sides of the evaluation point
  - Higher-order formulas: Use more points for better accuracy
- **Step Size (h)**: The spacing between points where function values are evaluated
- **Derivation**: Based on Taylor series expansion

### Mathematical Foundation: Taylor Series

The Taylor series expansion of a function f(x) about a point a is:

{% raw %}
$$f(x) = f(a) + f'(a)(x-a) + \frac{f''(a)}{2!}(x-a)^2 + \frac{f'''(a)}{3!}(x-a)^3 + \ldots$$
{% endraw %}

When we evaluate at x = a+h:

{% raw %}
$$f(a+h) = f(a) + f'(a)h + \frac{f''(a)}{2!}h^2 + \frac{f'''(a)}{3!}h^3 + \ldots$$
{% endraw %}

This forms the basis for deriving finite difference formulas by algebraically manipulating these expansions to isolate derivatives.

### Common Finite Difference Formulas

#### First Derivative Approximations

**Forward Difference (First-Order):**
{% raw %}
$$f'(x) \approx \frac{f(x+h) - f(x)}{h} + O(h)$$
{% endraw %}

**Backward Difference (First-Order):**
{% raw %}
$$f'(x) \approx \frac{f(x) - f(x-h)}{h} + O(h)$$
{% endraw %}

**Central Difference (Second-Order):**
{% raw %}
$$f'(x) \approx \frac{f(x+h) - f(x-h)}{2h} + O(h^2)$$
{% endraw %}

**Three-Point Formula (Second-Order):**
{% raw %}
$$f'(x) \approx \frac{-3f(x) + 4f(x+h) - f(x+2h)}{2h} + O(h^2)$$
{% endraw %}

#### Second Derivative Approximations

**Central Difference for Second Derivative (Second-Order):**
{% raw %}
$$f''(x) \approx \frac{f(x+h) - 2f(x) + f(x-h)}{h^2} + O(h^2)$$
{% endraw %}

### Higher-Order Formulas

For increased accuracy, higher-order formulas use more points:

**Fourth-Order Central Difference:**
{% raw %}
$$f'(x) \approx \frac{-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)}{12h} + O(h^4)$$
{% endraw %}

**Fourth-Order Formula for Second Derivative:**
{% raw %}
$$f''(x) \approx \frac{-f(x+2h) + 16f(x+h) - 30f(x) + 16f(x-h) - f(x-2h)}{12h^2} + O(h^4)$$
{% endraw %}

These higher-order formulas reduce truncation error but require more function evaluations and may be more sensitive to noise.

### Basic Implementation in Python

{% raw %}
```python
import numpy as np

def forward_difference(f, x, h):
    return (f(x+h) - f(x)) / h

def central_difference(f, x, h):
    return (f(x+h) - f(x-h)) / (2*h)

def second_derivative(f, x, h):
    return (f(x+h) - 2*f(x) + f(x-h)) / h**2

# Example
def f(x): return np.sin(x)
x0 = np.pi/4
h = 0.01

print(f\"Forward difference: {forward_difference(f, x0, h)}\")
print(f\"Central difference: {central_difference(f, x0, h)}\")
print(f\"Second derivative: {second_derivative(f, x0, h)}\")
print(f\"Exact first derivative: {np.cos(x0)}\")
print(f\"Exact second derivative: {-np.sin(x0)}\")
```
{% endraw %}

### Applying Finite Differences to Data Arrays

For discrete data points (x₁, y₁), (x₂, y₂), ..., (xₙ, yₙ), where the spacing may not be uniform:

**First Derivative at Interior Points:**
{% raw %}
$$f'(x_i) \approx \frac{y_{i+1} - y_{i-1}}{x_{i+1} - x_{i-1}}$$
{% endraw %}

**At Endpoints (Forward and Backward):**
{% raw %}
$$f'(x_1) \approx \frac{y_2 - y_1}{x_2 - x_1}$$
{% endraw %}

{% raw %}
$$f'(x_n) \approx \frac{y_n - y_{n-1}}{x_n - x_{n-1}}$$
{% endraw %}

### Using NumPy's gradient Function

For practical applications, NumPy's `gradient()` function efficiently computes numerical derivatives:

{% raw %}
```python
import numpy as np

# Generate data
x = np.linspace(0, 2*np.pi, 100)
y = np.sin(x)

# Calculate gradient
dy_dx = np.gradient(y, x)  # First argument is values, second is coordinates
```
{% endraw %}

## Resources
 - **NumPy Documentation**
   - [numpy.gradient](https://numpy.org/doc/stable/reference/generated/numpy.gradient.html)  

## 3. Truncation Error and Order of Accuracy

### Key Concepts

- **Truncation Error**: Error resulting from approximating continuous operations by discrete operations
- **Taylor Series Expansion**: The mathematical foundation for analyzing truncation error
- **Order of Accuracy**: The rate at which the error decreases as the step size decreases
- **Big O Notation**: Used to express the leading term in the error equation

### Mathematical Analysis of Truncation Error

The truncation error arises from neglecting higher-order terms in the Taylor series. For the forward difference formula:

Starting with the Taylor expansion:
{% raw %}
$$f(x+h) = f(x) + hf'(x) + \frac{h^2}{2}f''(x) + \frac{h^3}{6}f'''(x) + O(h^4)$$
{% endraw %}

Rearranging to isolate $f'(x)$:
{% raw %}
$$f'(x) = \frac{f(x+h) - f(x)}{h} - \frac{h}{2}f''(x) - \frac{h^2}{6}f'''(x) + O(h^3)$$
{% endraw %}

So the approximation {% raw %}$f'(x) \approx \frac{f(x+h) - f(x)}{h}${% endraw %} has a truncation error of {% raw %}$-\frac{h}{2}f''(x) + O(h^2)${% endraw %}.

The leading term containing $h$ determines the order of accuracy. Here, it's $O(h)$, making it a first-order method.

### Error Analysis for Common Formulas

**Forward Difference (First-Order):**
{% raw %}
$$\frac{f(x+h) - f(x)}{h} = f'(x) + \frac{h}{2}f''(x) + O(h^2)$$
{% endraw %}
Error: $O(h)$

**Central Difference (Second-Order):**
{% raw %}
$$\frac{f(x+h) - f(x-h)}{2h} = f'(x) + \frac{h^2}{6}f'''(x) + O(h^4)$$
{% endraw %}
Error: $O(h^2)$

Notice how the central difference formula eliminates the $h$ term, making it second-order accurate. This is why central differences are generally preferred when possible.

### Total Error: Balancing Truncation and Round-off Error

The total error in numerical differentiation combines:

1. **Truncation Error**: {% raw %}$E_t \approx Ch^p${% endraw %} (decreases with smaller h)
2. **Round-off Error**: {% raw %}$E_r \approx \frac{M\varepsilon}{h}${% endraw %} (increases with smaller h)

Where:
- $p$ is the order of the method
- $C$ is a constant depending on the derivatives of $f$
- $\varepsilon$ is the machine epsilon (floating-point precision)
- $M$ is a constant

The total error is:
{% raw %}
$$E_{total}(h) \approx Ch^p + \frac{M\varepsilon}{h}$$
{% endraw %}

### Optimal Step Size

To minimize the total error, we differentiate with respect to $h$ and set to zero:
{% raw %}
$$\frac{dE_{total}}{dh} = pCh^{p-1} - \frac{M\varepsilon}{h^2} = 0$$
{% endraw %}

Solving for $h$:
{% raw %}
$$h_{optimal} \approx \left( \frac{M\varepsilon}{pC} \right)^{1/(p+1)}$$
{% endraw %}

This shows that:
1. Higher-order methods (larger $p$) allow for larger optimal step sizes
2. Higher machine precision (smaller $\varepsilon$) allows for smaller step sizes

## 4. Differentiating Splines in SciPy

### Key Concepts

- **Spline Interpolation**: Creating a smooth curve through discrete data points
- **Spline Differentiation**: Computing derivatives of the spline function
- **SciPy Implementation**: Using `scipy.interpolate` for both interpolation and differentiation
- **Advantages**: Smoother derivatives compared to direct finite differences on noisy data

### Mathematical Background of Splines

A spline of degree $k$ is a piecewise polynomial function where pieces are connected with continuity in derivatives up to order $k-1$. The most common are cubic splines (degree 3), which maintain continuity in the function and its first and second derivatives.

For $n+1$ data points $(x_i, y_i)$ where $i = 0, 1, ..., n$, a cubic spline $S(x)$ consists of $n$ cubic polynomials $S_i(x)$ defined on intervals $[x_i, x_{i+1}]$:

{% raw %}
$$S_i(x) = a_i + b_i(x-x_i) + c_i(x-x_i)^2 + d_i(x-x_i)^3$$
{% endraw %}

The coefficients $a_i, b_i, c_i, d_i$ are determined by:
1. Interpolation: $S_i(x_i) = y_i$ and $S_i(x_{i+1}) = y_{i+1}$
2. Continuity: $S_i'(x_{i+1}) = S_{i+1}'(x_i)$ and $S_i''(x_{i+1}) = S_{i+1}''(x_i)$
3. Boundary conditions (commonly, natural: $S''(x_0) = S''(x_n) = 0$)

### Differentiation of Cubic Splines

Once a cubic spline is constructed, its derivative at any point is:

{% raw %}
$$S_i'(x) = b_i + 2c_i(x-x_i) + 3d_i(x-x_i)^2$$
{% endraw %}

And its second derivative is:

{% raw %}
$$S_i''(x) = 2c_i + 6d_i(x-x_i)$$
{% endraw %}

This provides a smooth, continuous approximation of the first and second derivatives.

### Spline Differentiation in SciPy

SciPy's `CubicSpline` class allows easy differentiation:

{% raw %}
```python
import numpy as np
from scipy.interpolate import CubicSpline

# Generate sample data
x = np.linspace(0, 10, 15)
y = np.sin(x) + 0.1 * np.random.randn(len(x))  # Add some noise

# Create cubic spline
cs = CubicSpline(x, y)

# Evaluate points and derivatives
x_fine = np.linspace(0, 10, 100)
y_spline = cs(x_fine)                # Function values
y_prime = cs(x_fine, 1)              # First derivative (nu=1)
y_second = cs(x_fine, 2)             # Second derivative (nu=2)
```
{% endraw %}

### Advantages for Noisy Data

Spline differentiation provides several benefits for noisy data:

1. **Smoothing**: Reduces noise amplification by fitting a smooth function
2. **Global Information**: Uses information from multiple points, not just neighbors
3. **Analytic Differentiation**: Differentiates the function analytically rather than numerically
4. **Consistency**: Derivatives are continuous across the entire domain

### Resources

- **SciPy Documentation**
  - [scipy.interpolate](https://docs.scipy.org/doc/scipy/reference/interpolate.html)
  - [CubicSpline](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CubicSpline.html)

## 5. Introduction to Numerical Integration

### Key Concepts

- **Numerical Integration**: Approximating definite integrals using function values at discrete points
- **Quadrature**: The numerical computation of an integral
- **Applications**:
  - Volume and area calculations
  - Work and energy computations
  - Probability and statistics
  - Signal processing
- **Advantages over Analytical Integration**:
  - Can handle functions without closed-form antiderivatives
  - Works with tabular data
  - Often simple to implement and computationally efficient

### Mathematical Foundation

A definite integral is defined as:

{% raw %}
$$\int_{a}^{b} f(x) \, dx = \lim_{n \to \infty} \sum_{i=1}^{n} f(x_i^*) \Delta x$$
{% endraw %}

Where {% raw %}$\Delta x = \frac{b-a}{n}${% endraw %} and $x_i^*$ is any point in the interval $[x_{i-1}, x_i]$.

Numerical integration methods approximate this sum using finite values of $n$ and specific choices of evaluation points.

### Riemann Sums

The simplest form of numerical integration uses Riemann sums:

1. **Left Riemann Sum**: {% raw %}$\sum_{i=1}^{n} f(x_{i-1}) \Delta x${% endraw %}
2. **Right Riemann Sum**: {% raw %}$\sum_{i=1}^{n} f(x_i) \Delta x${% endraw %}
3. **Midpoint Riemann Sum**: {% raw %}$\sum_{i=1}^{n} f\left(\frac{x_{i-1}+x_i}{2}\right) \Delta x${% endraw %}

These are first-order methods with error proportional to the step size.

### Classification of Integration Methods

Numerical integration methods can be classified by:

1. **Newton-Cotes Formulas**:
   - Based on polynomial interpolation of equally spaced points
   - Open formulas: Don't include endpoints
   - Closed formulas: Include endpoints
   - Examples: Trapezoid rule, Simpson's rules

2. **Gaussian Quadrature**:
   - Based on orthogonal polynomials
   - Uses optimally placed evaluation points and weights
   - Higher accuracy with fewer function evaluations
   - Examples: Gauss-Legendre, Gauss-Hermite

3. **Adaptive Methods**:
   - Dynamically adjust the step size based on local error estimates
   - Concentrate evaluations where needed
   - Examples: Adaptive Simpson's, Gauss-Kronrod

### Resources

- **Berkeley Numerical Methods**
  - [Numerical Integration](https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter21.00-Numerical-Integration.html)