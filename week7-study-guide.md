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