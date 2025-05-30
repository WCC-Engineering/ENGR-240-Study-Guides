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

print(f"Forward difference: {forward_difference(f, x0, h)}")
print(f"Central difference: {central_difference(f, x0, h)}")
print(f"Second derivative: {second_derivative(f, x0, h)}")
print(f"Exact first derivative: {np.cos(x0)}")
print(f"Exact second derivative: {-np.sin(x0)}")
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

## 6. Newton-Cotes Formulas

### Key Concepts

- **Newton-Cotes Formulas**: Family of integration methods based on polynomial interpolation
- **Common Methods**:
  - Trapezoid rule (linear interpolation)
  - Simpson's 1/3 rule (quadratic interpolation)
  - Simpson's 3/8 rule (cubic interpolation)
- **Composite Rules**: Apply basic formulas over multiple subintervals for better accuracy

### Mathematical Derivation

Newton-Cotes formulas approximate {% raw %}$\int_a^b f(x) dx${% endraw %} by:
1. Dividing $[a, b]$ into $n$ subintervals
2. Fitting a polynomial through points in each subinterval
3. Integrating the polynomial analytically

The general form is:
{% raw %}
$$\int_a^b f(x) dx \approx (b-a) \sum_{i=0}^n w_i f(x_i)$$
{% endraw %}

Where $w_i$ are weights derived from the integral of the Lagrange basis polynomials.

### Trapezoid Rule

The trapezoid rule approximates the integral by connecting function values with straight lines:

**Single Interval Formula:**
{% raw %}
$$\int_{a}^{b} f(x) \, dx \approx \frac{b-a}{2} [f(a) + f(b)]$$
{% endraw %}

**Composite Formula (n subintervals):**
{% raw %}
$$\int_{a}^{b} f(x) \, dx \approx \frac{h}{2} [f(x_0) + 2f(x_1) + 2f(x_2) + \ldots + 2f(x_{n-1}) + f(x_n)]$$
{% endraw %}

where {% raw %}$h = \frac{b-a}{n}${% endraw %} and $x_i = a + ih$.

**Error Term:**
{% raw %}
$$E_{trap} = -\frac{(b-a)^3}{12n^2} f''(\xi)$$
{% endraw %}

for some $\xi \in [a,b]$, making it a second-order method $O(h^2)$.

### Simpson's 1/3 Rule

Simpson's 1/3 rule fits parabolas through three points:

**Single Interval Formula (over two subintervals):**
{% raw %}
$$\int_{a}^{b} f(x) \, dx \approx \frac{b-a}{6} [f(a) + 4f\left(\frac{a+b}{2}\right) + f(b)]$$
{% endraw %}

**Composite Formula (n must be even):**
{% raw %}
$$\int_{a}^{b} f(x) \, dx \approx \frac{h}{3} [f(x_0) + 4f(x_1) + 2f(x_2) + 4f(x_3) + \ldots + 2f(x_{n-2}) + 4f(x_{n-1}) + f(x_n)]$$
{% endraw %}

**Error Term:**
{% raw %}
$$E_{simp} = -\frac{(b-a)^5}{180n^4} f^{(4)}(\xi)$$
{% endraw %}

making it a fourth-order method $O(h^4)$.

### Simpson's 3/8 Rule

Simpson's 3/8 rule uses cubic polynomials through four points:

**Single Interval Formula (over three subintervals):**
{% raw %}
$$\int_{a}^{b} f(x) \, dx \approx \frac{b-a}{8} [f(a) + 3f\left(\frac{2a+b}{3}\right) + 3f\left(\frac{a+2b}{3}\right) + f(b)]$$
{% endraw %}

**Composite Formula (n must be multiple of 3):**
{% raw %}
$$\int_{a}^{b} f(x) \, dx \approx \frac{3h}{8} [f(x_0) + 3f(x_1) + 3f(x_2) + 2f(x_3) + 3f(x_4) + \ldots + f(x_n)]$$
{% endraw %}

**Error Term:**
{% raw %}
$$E_{3/8} = -\frac{(b-a)^5}{80n^4} f^{(4)}(\xi)$$
{% endraw %}

making it also a fourth-order method $O(h^4)$.

### Basic Implementation Example

{% raw %}
```python
import numpy as np
from scipy import integrate

def f(x):
    return np.sin(x)**2

# Define integration limits
a, b = 0, np.pi
n = 10  # Number of subintervals

# Exact integral value for comparison
exact = integrate.quad(f, a, b)[0]

# Calculate using different methods
x = np.linspace(a, b, n+1)
y = f(x)

# Using NumPy/SciPy built-in functions
trapz_result = np.trapz(y, x)
simpson_result = integrate.simpson(y, x)
quad_result = integrate.quad(f, a, b)[0]

print(f"Trapezoid Rule: {trapz_result}, Error: {abs(trapz_result-exact)}")
print(f"Simpson's Rule: {simpson_result}, Error: {abs(simpson_result-exact)}")
print(f"Adaptive Quadrature: {quad_result}, Error: {abs(quad_result-exact)}")
```
{% endraw %}

### Error Analysis and Convergence

The error of Newton-Cotes formulas depends on:
1. The step size {% raw %}$h = \frac{b-a}{n}${% endraw %}
2. The smoothness of the integrand (higher derivatives)
3. The order of the method

For smooth functions:
- Trapezoid rule error scales as $O(h^2)$ or $O(n^{-2})$
- Simpson's rules error scales as $O(h^4)$ or $O(n^{-4})$

This rapid convergence for Simpson's rules makes them very efficient for smooth functions.

### Resources

- **Davishahl Numerical Methods Videos**
  - [Newton Cotes Integration](https://youtu.be/IH4wZ-cahXM?si=pzoZ5nHw3d2nZLqd)

- **Berkeley Numerical Methods**
  - [Trapezoidal Rule](https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter21.03-Trapezoid-Rule.html)
  - [Simpson's Rule](https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter21.04-Simpsons-Rule.html)

## 7. Gauss Quadrature

Gauss quadrature is an advanced numerical integration technique that achieves remarkable accuracy by strategically choosing both the evaluation points and their weights. Unlike other methods that use predetermined points (like endpoints or equally spaced intervals), Gauss quadrature optimizes both the locations where we evaluate the function and how much weight we give each evaluation.

### 7.1 Conceptual Foundation

The key insight behind Gauss quadrature is **error balancing**. Instead of simply connecting endpoints with a straight line (which creates systematic error), Gauss quadrature chooses evaluation points that balance positive and negative errors, leading to much higher accuracy.

**Visual Understanding:**
- **Trapezoid method**: Uses endpoints, creating systematic over- or under-estimation
- **2-point Gauss**: Chooses two interior points that balance positive and negative errors
- **3-point Gauss**: Uses three strategically placed points for even higher accuracy

This error balancing allows an n-point Gauss quadrature formula to integrate polynomials of degree 2n-1 exactly, which is remarkable efficiency.

### 7.2 Gauss-Legendre Formulas

Gauss quadrature approximates integrals using the form:

$$I \approx c_0 f(x_0) + c_1 f(x_1) + \cdots + c_{n-1} f(x_{n-1}) = \sum_{i=0}^{n-1} c_i f(x_i)$$

where both the weights $c_i$ and the evaluation points $x_i$ are optimally determined.

**Standard Form**: Gauss-Legendre quadrature is derived for the interval [-1, 1]. For other intervals [a, b], we need coordinate transformation.

### 7.3 Gauss-Legendre Constants

The following table provides the weights and evaluation points for common Gauss-Legendre formulas:

| **Points** | **Weighting Factors, $c_i$** | **Function Arguments, $x_i$** |
|------------|------------------------------|--------------------------------|
| 2 | $c_0 = 1.0$ | $x_0 = -1/\sqrt{3} \approx -0.5774$ |
|   | $c_1 = 1.0$ | $x_1 = 1/\sqrt{3} \approx 0.5774$ |
| 3 | $c_0 = 5/9 \approx 0.5556$ | $x_0 = -\sqrt{3/5} \approx -0.7746$ |
|   | $c_1 = 8/9 \approx 0.8889$ | $x_1 = 0.0$ |
|   | $c_2 = 5/9 \approx 0.5556$ | $x_2 = \sqrt{3/5} \approx 0.7746$ |
| 4 | $c_0 = (18-\sqrt{30})/36 \approx 0.3479$ | $x_0 = -\sqrt{525+70\sqrt{30}}/35 \approx -0.8612$ |
|   | $c_1 = (18+\sqrt{30})/36 \approx 0.6521$ | $x_1 = -\sqrt{525-70\sqrt{30}}/35 \approx -0.3400$ |
|   | $c_2 = (18+\sqrt{30})/36 \approx 0.6521$ | $x_2 = \sqrt{525-70\sqrt{30}}/35 \approx 0.3400$ |
|   | $c_3 = (18-\sqrt{30})/36 \approx 0.3479$ | $x_3 = \sqrt{525+70\sqrt{30}}/35 \approx 0.8612$ |

### 7.4 Coordinate Transformation

Since Gauss-Legendre formulas are defined for [-1, 1], we must transform integrals over [a, b]:

$$\int_a^b f(x) dx \Rightarrow \int_{-1}^1 f(x_d) dx_d$$

**Linear Transformation**: If $x_d$ is the variable on [-1, 1] and $x$ is the variable on [a, b], then:

$$x = a_1 + a_2 x_d$$

**Finding the coefficients**:
- When $x_d = -1$: $x = a$ → $a_1 + a_2(-1) = a$
- When $x_d = 1$: $x = b$ → $a_1 + a_2(1) = b$

Solving these equations:
- $a_1 = \frac{a + b}{2}$
- $a_2 = \frac{b - a}{2}$

Therefore: $x = \frac{a + b}{2} + \frac{b - a}{2} x_d$

And: $dx = \frac{b - a}{2} dx_d$

### 7.5 Complete Example

**Problem**: Use three-point Gauss-Legendre quadrature to estimate:
$$I = \int_1^3 \left(2x + \frac{3}{x}\right)^2 dx$$

**Step 1**: Transform to [-1, 1]
- $a = 1$, $b = 3$
- $a_1 = \frac{1 + 3}{2} = 2$
- $a_2 = \frac{3 - 1}{2} = 1$
- $x = 2 + x_d$, so $dx = dx_d$

**Step 2**: Apply three-point Gauss-Legendre formula
From the table:
- $x_0 = -\sqrt{3/5} \approx -0.7746$, $c_0 = 5/9$
- $x_1 = 0.0$, $c_1 = 8/9$  
- $x_2 = \sqrt{3/5} \approx 0.7746$, $c_2 = 5/9$

**Step 3**: Calculate function values
Transform the x-coordinates:
- At $x_d = -\sqrt{3/5}$: $x = 2 + (-\sqrt{3/5}) = 2 - \sqrt{3/5} \approx 1.2254$
- At $x_d = 0$: $x = 2 + 0 = 2.0$
- At $x_d = \sqrt{3/5}$: $x = 2 + \sqrt{3/5} \approx 2.7746$

Evaluate $f(x) = \left(2x + \frac{3}{x}\right)^2$:
- $f(1.2254) = (2(1.2254) + 3/1.2254)^2 = (2.451 + 2.449)^2 = 24.01$
- $f(2.0) = (2(2.0) + 3/2.0)^2 = (4.0 + 1.5)^2 = 30.25$
- $f(2.7746) = (2(2.7746) + 3/2.7746)^2 = (5.549 + 1.081)^2 = 43.96$

**Step 4**: Apply the formula
$$I = \frac{5}{9} f(1.2254) + \frac{8}{9} f(2.0) + \frac{5}{9} f(2.7746)$$
$$I = \frac{5}{9}(24.01) + \frac{8}{9}(30.25) + \frac{5}{9}(43.96)$$
$$I = 13.34 + 26.89 + 24.42 = 64.65$$

**Verification**: The exact value is $I_{true} = 64.667$, so our estimate has excellent accuracy!

### 7.6 Python Implementation

```python
import numpy as np

def gauss_legendre_2point(f, a, b):
    """Two-point Gauss-Legendre quadrature"""
    # Gauss-Legendre constants for 2 points
    x_vals = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
    weights = np.array([1.0, 1.0])
    
    # Transform to [a, b]
    x_transformed = (a + b)/2 + (b - a)/2 * x_vals
    
    # Calculate integral
    integral = (b - a)/2 * np.sum(weights * f(x_transformed))
    return integral

def gauss_legendre_3point(f, a, b):
    """Three-point Gauss-Legendre quadrature"""
    # Gauss-Legendre constants for 3 points
    x_vals = np.array([-np.sqrt(3/5), 0.0, np.sqrt(3/5)])
    weights = np.array([5/9, 8/9, 5/9])
    
    # Transform to [a, b]
    x_transformed = (a + b)/2 + (b - a)/2 * x_vals
    
    # Calculate integral
    integral = (b - a)/2 * np.sum(weights * f(x_transformed))
    return integral

# Example usage
def f(x):
    return (2*x + 3/x)**2

# Compare methods
result_2pt = gauss_legendre_2point(f, 1, 3)
result_3pt = gauss_legendre_3point(f, 1, 3)

print(f"2-point Gauss-Legendre: {result_2pt:.4f}")
print(f"3-point Gauss-Legendre: {result_3pt:.4f}")
print(f"Exact value: 64.6667")
```

### 7.7 Key Advantages

1. **High Accuracy**: n-point Gauss quadrature integrates polynomials of degree 2n-1 exactly
2. **Efficiency**: Fewer function evaluations needed compared to other methods for similar accuracy
3. **Optimal Point Placement**: Points are strategically chosen to minimize error
4. **Error Balancing**: Positive and negative errors cancel out effectively

### 7.8 When to Use Gauss Quadrature

**Best for**:
- Smooth functions over finite intervals
- When high accuracy is needed with minimal function evaluations
- Integrands that are well-approximated by polynomials

**Consider alternatives for**:
- Functions with discontinuities or sharp peaks
- Infinite intervals (use other specialized methods)
- When you need to evaluate the function at specific predetermined points

### Resources

- **Davishahl Numerical Methods Videos**
  - [Gauss Quadrature](https://youtu.be/yF-JlQxaZQs?si=Kxrc6qlQvAw1xD8s)

## 8. Using NumPy and SciPy for Differentiation and Integration

### Key Functions in NumPy for Differentiation

NumPy provides gradient calculation through `numpy.gradient()`:

{% raw %}
```python
import numpy as np

# Example with evenly spaced data
x = np.linspace(0, 2*np.pi, 100)
y = np.sin(x)

# Calculate gradient (first derivative)
dy_dx = np.gradient(y, x)  # First argument is values, second is coordinates

# For 2D data (partial derivatives)
X, Y = np.meshgrid(x, x)
Z = np.sin(X) * np.cos(Y)
dZ_dx, dZ_dy = np.gradient(Z, x, x)
```
{% endraw %}

Key features of `numpy.gradient()`:
- Handles 1D, 2D, and higher-dimensional arrays
- Supports both uniform and non-uniform spacing
- Uses second-order central differences for interior points
- Uses forward/backward differences at boundaries
- Returns gradients along each axis

The mathematical implementation of `numpy.gradient()` for non-uniform spacing is:

{% raw %}
$$f'(x_i) \approx \frac{f(x_{i+1}) - f(x_i)}{x_{i+1} - x_i} \cdot \frac{x_{i+1} - x_i}{x_{i+1} - x_{i-1}} + \frac{f(x_i) - f(x_{i-1})}{x_i - x_{i-1}} \cdot \frac{x_i - x_{i-1}}{x_{i+1} - x_{i-1}}$$
{% endraw %}

which is a weighted average of forward and backward differences.

### Key Functions in SciPy for Integration

SciPy's `scipy.integrate` module offers various integration methods:

{% raw %}
```python
import numpy as np
from scipy import integrate

def f(x):
    return np.exp(-x**2) * np.sin(2*x)

# General-purpose adaptive integration
result, error = integrate.quad(f, 0, 3)

# Integration of fixed samples
x = np.linspace(0, 3, 100)
y = f(x)
trapz_result = integrate.trapz(y, x)
simpson_result = integrate.simpson(y, x)

# Cumulative integration
cumulative = integrate.cumtrapz(y, x, initial=0)

# For double integrals
def g(x, y):
    return np.sin(x + y)

result_dblquad = integrate.dblquad(g, 0, 1, lambda x: 0, lambda x: 1)
```
{% endraw %}

Key integration functions in SciPy:

1. **integrate.quad**: General-purpose adaptive integration
   - Based on QUADPACK Fortran library
   - Automatically handles singularities and oscillatory integrands
   - Provides error estimates
   - Parameters for controlling accuracy and evaluation limits

2. **integrate.trapz, integrate.simpson**: Integration of tabular data
   - Implement composite trapezoid and Simpson's rules
   - Work directly with arrays of values
   - Suitable for experimental data or pre-computed function values

3. **integrate.cumtrapz**: Cumulative integration (indefinite integral)
   - Returns the running integral at each point
   - Useful for analyzing accumulated effects

4. **integrate.dblquad, integrate.tplquad**: Multi-dimensional integration
   - Handle double and triple integrals
   - Support variable bounds
   - Based on adaptive quadrature

5. **integrate.fixed_quad**: Non-adaptive Gaussian quadrature
   - Uses Gauss-Legendre quadrature with specified order
   - Very efficient for smooth functions

6. **integrate.romberg**: Romberg integration
   - Uses Richardson extrapolation on the trapezoid rule
   - Efficient for smooth functions

### Mathematical Comparison of Integration Methods

For a function with continuous derivatives up to order $p$, the error behavior of different methods can be summarized as follows:

| Method | Error Behavior | Function Evaluations | Best For |
|--------|---------------|----------------------|----------|
| Trapezoid Rule | $O(h^2)$ | $n+1$ | Simple implementation, low accuracy |
| Simpson's 1/3 | $O(h^4)$ | $n+1$ (n even) | Good balance of accuracy and simplicity |
| Simpson's 3/8 | $O(h^4)$ | $n+1$ (n multiple of 3) | Similar to Simpson's 1/3 |
| Gauss-Legendre (n points) | $O(h^{2n})$ | $n$ | Highest accuracy per function evaluation |
| Adaptive Quadrature | Varies | Varies | Unknown behavior, singularities |

### Error Control in SciPy's quad

SciPy's `quad` function provides sophisticated error control through:

1. **Absolute Tolerance (epsabs)**: Maximum absolute error allowed
2. **Relative Tolerance (epsrel)**: Maximum relative error allowed
3. **Limit (limit)**: Maximum number of subintervals
4. **Points (points)**: Explicitly specified evaluation points

The algorithm subdivides the interval adaptively until error estimates meet the specified tolerances.

The mathematical principles behind adaptive quadrature include:
- Estimate local error by comparing results from different-order rules
- Subdivide regions with highest estimated error
- Allocate more points where the function changes rapidly
- Combine results using appropriate rules to preserve accuracy

### Example: When to Use Different Methods

1. **Use quad for general-purpose integration**:
   - Unknown function behavior
   - Need for error estimates
   - Highest reliability

2. **Use Simpson's rule for smooth, tabular data**:
   - Pre-computed function values
   - Known smooth behavior
   - Moderate accuracy requirements

3. **Use Gauss-Legendre for high-precision needs**:
   - Very smooth functions
   - High accuracy requirements
   - Minimal function evaluations

4. **Use adaptive methods for challenging integrands**:
   - Highly oscillatory functions
   - Functions with sharp peaks
   - Integrals over large or infinite domains

### The Mathematical Relationship Between Differentiation and Integration

The Fundamental Theorem of Calculus establishes the relationship between differentiation and integration:

For a continuous function $f$ on $[a,b]$ and the antiderivative {% raw %}$F(x) = \int_a^x f(t) dt${% endraw %}:

{% raw %}
$$\frac{d}{dx}F(x) = \frac{d}{dx}\int_a^x f(t) dt = f(x)$$
{% endraw %}

This relationship has important implications for numerical methods:
- Integration is generally well-conditioned (smooths errors)
- Differentiation is ill-conditioned (amplifies errors)
- Numerical differentiation of an integral can reduce error
- Numerical integration of a derivative often improves accuracy

### Resources
- **Berkeley Numerical Methods**
  - [Computing Integrals in Python](https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter21.05-Computing-Integrals-in-Python.html)

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

4. **Mathematical insights**:
   - Differentiation is inherently ill-conditioned: small errors in function values lead to large errors in derivatives
   - Integration is well-conditioned: small errors in function values lead to small errors in integrals
   - Higher-order methods give better accuracy but might be more sensitive to noise
   - For practical applications, built-in NumPy and SciPy functions are usually preferred over manual implementations

5. **Documentation best practices**:
   - Clearly state the methods used and their expected order of accuracy
   - Document assumptions about the input data
   - Include error estimates when presenting numerical results
