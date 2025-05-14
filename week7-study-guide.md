}{12n^2} f''(\xi)$$

for some $\xi \in [a,b]$, making it a second-order method $O(h^2)$.

### Simpson's 1/3 Rule

Simpson's 1/3 rule fits parabolas through three points:

**Single Interval Formula (over two subintervals):**
$$\int_{a}^{b} f(x) \, dx \approx \frac{b-a}{6} [f(a) + 4f\left(\frac{a+b}{2}\right) + f(b)]$$

**Composite Formula (n must be even):**
$$\int_{a}^{b} f(x) \, dx \approx \frac{h}{3} [f(x_0) + 4f(x_1) + 2f(x_2) + 4f(x_3) + \ldots + 2f(x_{n-2}) + 4f(x_{n-1}) + f(x_n)]$$

**Error Term:**
$$E_{simp} = -\frac{(b-a)^5}{180n^4} f^{(4)}(\xi)$$

making it a fourth-order method $O(h^4)$.

### Simpson's 3/8 Rule

Simpson's 3/8 rule uses cubic polynomials through four points:

**Single Interval Formula (over three subintervals):**
$$\int_{a}^{b} f(x) \, dx \approx \frac{b-a}{8} [f(a) + 3f\left(\frac{2a+b}{3}\right) + 3f\left(\frac{a+2b}{3}\right) + f(b)]$$

**Composite Formula (n must be multiple of 3):**
$$\int_{a}^{b} f(x) \, dx \approx \frac{3h}{8} [f(x_0) + 3f(x_1) + 3f(x_2) + 2f(x_3) + 3f(x_4) + \ldots + f(x_n)]$$

**Error Term:**
$$E_{3/8} = -\frac{(b-a)^5}{80n^4} f^{(4)}(\xi)$$

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
1. The step size $h = \frac{b-a}{n}$
2. The smoothness of the integrand (higher derivatives)
3. The order of the method

For smooth functions:
- Trapezoid rule error scales as $O(h^2)$ or $O(n^{-2})$
- Simpson's rules error scales as $O(h^4)$ or $O(n^{-4})$

This rapid convergence for Simpson's rules makes them very efficient for smooth functions.

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

### Mathematical Theory

Gaussian quadrature is based on the idea that we can choose both the evaluation points $x_i$ and weights $w_i$ to maximize accuracy. For an $n$-point quadrature:

$$\int_{-1}^{1} f(x) \, dx \approx \sum_{i=1}^{n} w_i f(x_i)$$

For Gauss-Legendre quadrature:
- The points $x_i$ are the roots of the $n$th-degree Legendre polynomial $P_n(x)$
- The weights are given by $w_i = \frac{2}{(1-x_i^2)[P_n'(x_i)]^2}$

This scheme exactly integrates polynomials of degree $2n-1$ or less, making it extremely efficient.

### Legendre Polynomials

Legendre polynomials $P_n(x)$ are orthogonal polynomials on $[-1,1]$ with the property:

$$\int_{-1}^{1} P_m(x)P_n(x) \, dx = \frac{2}{2n+1}\delta_{mn}$$

They can be generated using the recurrence relation:
$$(n+1)P_{n+1}(x) = (2n+1)xP_n(x) - nP_{n-1}(x)$$

with $P_0(x) = 1$ and $P_1(x) = x$.

The first few Legendre polynomials are:
- $P_0(x) = 1$
- $P_1(x) = x$
- $P_2(x) = \frac{1}{2}(3x^2 - 1)$
- $P_3(x) = \frac{1}{2}(5x^3 - 3x)$
- $P_4(x) = \frac{1}{8}(35x^4 - 30x^2 + 3)$

### Integration on General Intervals

For an integral over $[a,b]$ instead of $[-1,1]$, we use the transformation:

$$\int_{a}^{b} f(x) \, dx = \frac{b-a}{2} \int_{-1}^{1} f\left(\frac{b-a}{2}t + \frac{a+b}{2}\right) \, dt$$

Applying Gauss-Legendre quadrature:

$$\int_{a}^{b} f(x) \, dx \approx \frac{b-a}{2} \sum_{i=1}^{n} w_i f\left(\frac{b-a}{2}x_i + \frac{a+b}{2}\right)$$

### Implementation with SciPy

SciPy provides Gauss-Legendre quadrature through `scipy.special`:

{% raw %}
```python
import numpy as np
from scipy import special, integrate

def f(x):
    return np.exp(-x**2)

# Integration limits
a, b = 0, 1

# Using SciPy's fixed_quad (Gauss-Legendre)
result, error = integrate.fixed_quad(f, a, b, n=5)
print(f"5-point Gauss-Legendre: {result}")

# For comparison
exact = integrate.quad(f, a, b)[0]
print(f"Exact: {exact}, Error: {abs(result-exact)}")
```
{% endraw %}

### Effectiveness of Gauss-Legendre Quadrature

Gauss-Legendre quadrature is remarkably efficient for smooth functions:

1. An $n$-point formula exactly integrates polynomials up to degree $2n-1$
2. Error for analytic functions decreases exponentially with $n$
3. For highly oscillatory or singular functions, specialized Gaussian rules may be more appropriate (e.g., Gauss-Chebyshev, Gauss-Laguerre)

The theoretical error term for Gauss-Legendre quadrature is:

$$E_n = \frac{2^{2n+1}(n!)^4}{(2n+1)[(2n)!]^3} f^{(2n)}(\xi)$$

for some $\xi \in [a,b]$. This displays super-polynomial convergence for smooth functions.

### Other Gaussian Quadrature Methods

Different weight functions lead to different Gaussian quadrature methods:

1. **Gauss-Chebyshev**: For integrals of the form $\int_{-1}^{1} \frac{f(x)}{\sqrt{1-x^2}} \, dx$
2. **Gauss-Laguerre**: For integrals of the form $\int_{0}^{\infty} e^{-x}f(x) \, dx$
3. **Gauss-Hermite**: For integrals of the form $\int_{-\infty}^{\infty} e^{-x^2}f(x) \, dx$
4. **Gauss-Jacobi**: For integrals of the form $\int_{-1}^{1} (1-x)^{\alpha}(1+x)^{\beta}f(x) \, dx$

Each is optimal for specific types of integrals and decay behaviors.

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

$$f'(x_i) \approx \frac{f(x_{i+1}) - f(x_i)}{x_{i+1} - x_i} \cdot \frac{x_{i+1} - x_i}{x_{i+1} - x_{i-1}} + \frac{f(x_i) - f(x_{i-1})}{x_i - x_{i-1}} \cdot \frac{x_i - x_{i-1}}{x_{i+1} - x_{i-1}}$$

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

For a continuous function $f$ on $[a,b]$ and the antiderivative $F(x) = \int_a^x f(t) dt$:

$$\frac{d}{dx}F(x) = \frac{d}{dx}\int_a^x f(t) dt = f(x)$$

This relationship has important implications for numerical methods:
- Integration is generally well-conditioned (smooths errors)
- Differentiation is ill-conditioned (amplifies errors)
- Numerical differentiation of an integral can reduce error
- Numerical integration of a derivative often improves accuracy

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

4. **Mathematical insights**:
   - Differentiation is inherently ill-conditioned: small errors in function values lead to large errors in derivatives
   - Integration is well-conditioned: small errors in function values lead to small errors in integrals
   - Higher-order methods give better accuracy but might be more sensitive to noise
   - For practical applications, built-in NumPy and SciPy functions are usually preferred over manual implementations

5. **Documentation best practices**:
   - Clearly state the methods used and their expected order of accuracy
   - Document assumptions about the input data
   - Include error estimates when presenting numerical results