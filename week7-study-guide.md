### Other Gaussian Quadrature Methods

Different weight functions lead to different Gaussian quadrature methods:

1. **Gauss-Chebyshev**: For integrals of the form {% raw %}$\int_{-1}^{1} \frac{f(x)}{\sqrt{1-x^2}} \, dx${% endraw %}
2. **Gauss-Laguerre**: For integrals of the form {% raw %}$\int_{0}^{\infty} e^{-x}f(x) \, dx${% endraw %}
3. **Gauss-Hermite**: For integrals of the form {% raw %}$\int_{-\infty}^{\infty} e^{-x^2}f(x) \, dx${% endraw %}
4. **Gauss-Jacobi**: For integrals of the form {% raw %}$\int_{-1}^{1} (1-x)^{\alpha}(1+x)^{\beta}f(x) \, dx${% endraw %}

Each is optimal for specific types of integrals and decay behaviors.

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