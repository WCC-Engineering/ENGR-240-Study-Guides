# Week 3 Study Guide: Root Finding Algorithms

## Overview

This study guide covers numerical methods for finding roots of nonlinear equations, a fundamental problem in engineering and scientific computing. The activities, worksheets, and programming assignments for Week 3 introduce several algorithms for root finding and methods for measuring error in numerical approximations.

## 1. Understanding Root Finding Problems

### Key Concepts

- **Definition of a Root**: A value of x where f(x) = 0
- **Root Finding Problem**: Given a function f(x), find the value(s) of x where f(x) = 0
- **Graphical Interpretation**: Roots are the x-intercepts of the function's graph
- **Common Engineering Applications**:
  - Solving nonlinear equations in thermodynamics
  - Finding equilibrium points in mechanical systems
  - Determining critical values in electrical circuits
  - Computing flow rates in fluid mechanics

### Formal Definitions and Measures of Accuracy

#### Root Finding Problem Statement

For a continuous function $f: \mathbb{R} \rightarrow \mathbb{R}$, find values $x_T$ such that:

$$f(x_T) = 0$$

A function may have:
- No roots
- A single root
- Multiple distinct roots
- Roots with multiplicity (where both $f(x)=0$ and $f'(x)=0$ at the root)

#### Bracketing a Root

If $f(a)$ and $f(b)$ have opposite signs (i.e., $f(a) \cdot f(b) < 0$), and $f$ is continuous on the interval $[a,b]$, then by the Intermediate Value Theorem, there exists at least one value $c \in (a,b)$ such that $f(c) = 0$.

#### Measures of Accuracy

When approximating a root $x_T$ with a numerical value $x_a$, we can measure the accuracy in several ways:

1. **True Error**: The absolute difference between the exact value and approximation
   $$E_t = |x_T - x_a|$$

2. **True Relative Error**: The true error normalized by the exact value
   $$\varepsilon_t = \left|\frac{x_T - x_a}{x_T}\right| \times 100\%$$

3. **Approximate Relative Error**: Since $x_T$ is usually unknown, we estimate the error using consecutive approximations
   $$\varepsilon_a = \left|\frac{x_{a,new} - x_{a,old}}{x_{a,new}}\right| \times 100\%$$

4. **Residual Error**: The absolute value of the function at the approximate root
   $$\varepsilon_r = |f(x_a)|$$

#### Stopping Criteria for Iterative Methods

Numerical root-finding algorithms typically use one or more of these stopping criteria:

1. **Relative Error Criterion**: Stop when $\varepsilon_a < \varepsilon_s$, where $\varepsilon_s$ is a specified tolerance

2. **Residual Criterion**: Stop when $|f(x_a)| < \delta$, where $\delta$ is a specified tolerance

3. **Maximum Iterations**: Stop after a specified number of iterations to prevent infinite loops

4. **Stagnation Criterion**: Stop when the change in successive iterations becomes too small
   
### Resources

- **Sullivan Numerical Methods Textbook**:
  - [Introduction to Root Finding](https://numericalmethodssullivan.github.io/ch-algebra.html#intro-to-numerical-root-finding)

## 2. Bisection Method

### Key Concepts

- **Principle**: If f(a) and f(b) have opposite signs, and f is continuous, there must be a root between a and b
- **Algorithm**:
  1. Start with an interval [a, b] where f(a) and f(b) have opposite signs
  2. Calculate the midpoint c = (a + b) / 2
  3. Evaluate f(c)
  4. If f(c) = 0 (or close enough), c is the root
  5. If f(c) has the same sign as f(a), update a = c; otherwise, update b = c
  6. Repeat steps 2-5 until convergence criteria are met
- **Convergence**: Always converges if initial conditions are valid
- **Rate of Convergence**: Linear (relatively slow)
- **Maximum Error Bound**: Error ≤ (b - a) / 2^n where n is the number of iterations

### Bisection Method Algorithm

1. **Initialize**:
   - Start with an interval [a, b] where f(a) and f(b) have opposite signs
   - Set error tolerance and maximum iterations

2. **Verify bracket validity**:
   - Check that f(a) × f(b) < 0 to ensure a root exists in the interval

3. **Iterate until convergence**:
   - Calculate midpoint c = (a + b) / 2
   - Evaluate f(c)
   - Calculate approximate relative error between current and previous midpoints
   - Check if error is less than tolerance or |f(c)| is sufficiently small
   - Update the interval:
     - If f(a) × f(c) < 0, the root is in [a, c], so set b = c
     - Otherwise, the root is in [c, b], so set a = c
   - Continue until error tolerance or maximum iterations reached

4. **Return results**:
   - Final approximation of the root
   - Number of iterations performed
   - Error convergence history

The bisection method is one of the most reliable root-finding techniques because it always converges for continuous functions if a valid initial bracket is provided. However, it typically requires more iterations than other methods due to its linear convergence rate.

### Resources

- **Davishahl Numerical Methods**
  - [Bisection](https://youtu.be/PdkCEp5M4J0?si=lyW0W1fFXKPdzLth)

- **Berkeley Python Numerical Methods**:
  - [Bracketing Methods](https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter19.03-Bisection-Method.html)

## 3. Newton-Raphson Method

### Key Concepts

- **Principle**: Uses function value and derivative to approximate roots via tangent lines
- **Formula**: $x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}$
- **Algorithm**:
  1. Start with an initial guess x₀
  2. Calculate the next approximation using the formula above
  3. Repeat until convergence criteria are met
- **Advantages**:
  - Converges quadratically when close to the root
  - Generally requires fewer iterations than bisection
- **Disadvantages**:
  - Requires the derivative of the function
  - May not converge if the initial guess is poor
  - Can fail if f'(x) is zero or very small

### Newton-Raphson Method Algorithm

1. **Initialize**:
   - Start with an initial guess x₀
   - Set error tolerance and maximum iterations

2. **Iterate until convergence**:
   - Calculate function value f(xₙ)
   - Calculate derivative value f'(xₙ)
   - Verify that f'(xₙ) is not too close to zero to avoid division problems
   - Calculate next approximation: xₙ₊₁ = xₙ - f(xₙ)/f'(xₙ)
   - Calculate approximate relative error: |xₙ₊₁ - xₙ|/|xₙ₊₁|
   - Check if error is less than tolerance or |f(xₙ₊₁)| is sufficiently small
   - Continue until error tolerance or maximum iterations reached

3. **Return results**:
   - Final approximation of the root
   - Number of iterations performed
   - Error convergence history

The Newton-Raphson method offers quadratic convergence, making it much faster than bisection when close to a root. However, it requires the derivative of the function and can fail if the derivative becomes zero or if the initial guess is poor. The method can also diverge for functions with complex behavior.

### Resources

- **Davishahl Numerical Methods**
  - [Newton-Raphson](https://youtu.be/ZSTSsG3I1z8?si=InK8iC9F7Uj54rn_)

- **Berkeley Python Numerical Methods**:]
  - [Newton-Raphson Method](https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter19.04-Newton-Raphson-Method.html)

## 4. Secant Method

### Key Concepts

- **Principle**: Approximates the derivative using two previous points
- **Formula**: $x_{i+1} = x_i - f(x_i) \frac{x_i - x_{i-1}}{f(x_i) - f(x_{i-1})}$
- **Algorithm**:
  1. Start with two initial guesses x₀ and x₁
  2. Calculate the next approximation using the formula above
  3. Repeat until convergence criteria are met
- **Advantages**:
  - Doesn't require the derivative of the function
  - Converges superlinearly (order of 1.618)
  - Generally faster than bisection
- **Disadvantages**:
  - Requires two initial points
  - May not converge if initial guesses are poor
  - Not as fast as Newton-Raphson when derivatives are available

### Secant Method Algorithm

{% raw %}
1. **Initialize**:
   - Start with two initial approximations x₀ and x₁
   - Evaluate f(x₀) and f(x₁)
   - Set error tolerance and maximum iterations

2. **Iterate until convergence**:
   - Check if f(x₁) - f(x₀) is close to zero to avoid division problems
   - Calculate next approximation: x₂ = x₁ - f(x₁)·(x₁ - x₀)/(f(x₁) - f(x₀))
   - Calculate approximate relative error: |x₂ - x₁|/|x₂|
   - Check if error is less than tolerance or |f(x₂)| is sufficiently small
   - Update values for next iteration: x₀ = x₁, x₁ = x₂, f(x₀) = f(x₁), f(x₁) = f(x₂)
   - Continue until error tolerance or maximum iterations reached

3. **Return results**:
   - Final approximation of the root
   - Number of iterations performed
   - Error convergence history
{% endraw %}

The secant method offers a good compromise between bisection and Newton-Raphson. It doesn't require the derivative, making it suitable for functions where derivatives are difficult to compute or unavailable. Its convergence rate is superlinear (order of approximately 1.618), which is slower than Newton-Raphson but faster than bisection.

### Resources

- **Davishahl Numerical Methods**
  - [Secant Method](https://youtu.be/IGB9V7mNOv8?si=pp5TJUmv9EcVoIiI)

## 5. Using SciPy for Root Finding

### Key Concepts

- **SciPy** provides optimized implementations of several root-finding algorithms
- **scipy.optimize.root_scalar**: Function for finding a root of a scalar function
- **scipy.optimize.fsolve**: General-purpose root-finding function
- **Available Methods**:
  - Bisection
  - Newton-Raphson ('newton')
  - Secant method
  - Brentq (combination of bisection, secant, and inverse quadratic interpolation)
  - Many others
- **Advantages**: Pre-implemented, well-tested, efficient algorithms

### SciPy Root Finding Overview

SciPy provides several optimized methods for root finding through the `scipy.optimize` module:

#### 1. Using `root_scalar` for Single-Variable Problems

`scipy.optimize.root_scalar` provides a unified interface to various root-finding algorithms:

```python
from scipy import optimize

# Basic syntax:
result = optimize.root_scalar(function, method='method_name', **method_specific_parameters)
```

**Common Methods and Their Parameters:**

- **Bisection Method**: 
  ```python
  result = optimize.root_scalar(f, method='bisect', bracket=[a, b])
  ```

- **Newton-Raphson Method**:
  ```python
  result = optimize.root_scalar(f, method='newton', x0=initial_guess, fprime=derivative_function)
  ```

- **Secant Method**:
  ```python
  result = optimize.root_scalar(f, method='secant', x0=guess1, x1=guess2)
  ```

- **Brentq Method** (Recommended general-purpose method):
  ```python
  result = optimize.root_scalar(f, method='brentq', bracket=[a, b])
  ```

**Result Properties:**
- `result.root`: The solution (root)
- `result.iterations`: Number of iterations needed
- `result.function_calls`: Number of function evaluations
- `result.converged`: Boolean indicating success

#### 2. Using `fsolve` for Multi-Variable Problems

For systems of nonlinear equations:

```python
solution = optimize.fsolve(system_function, initial_guess)
```

SciPy's implementations are generally more robust, efficient, and well-tested than custom implementations, making them preferable for most practical applications.

### Resources

- **SciPy Documentation**:
  - [scipy.optimize.root_scalar](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root_scalar.html)
  - [scipy.optimize.fsolve](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fsolve.html#scipy.optimize.fsolve)

## 6. Comparing Root Finding Methods

### Key Concepts

- **Selection Criteria**:
  - Convergence rate
  - Reliability
  - Implementation complexity
  - Derivative availability
  - Initial guess availability
- **Performance Comparison**:
  - Bisection: Slow but reliable, needs bracket
  - Newton-Raphson: Fast (quadratic) convergence, requires derivative
  - Secant: Good compromise, no derivative needed, needs two points
  - Brentq: Robust general purpose, needs bracket
- **Method Selection Guidelines**:
  - Use Bisection when a reliable method is needed and speed is not critical
  - Use Newton-Raphson when the derivative is available and a good initial guess is possible
  - Use Secant when the derivative is difficult to compute
  - Use Brentq when a robust general-purpose method is needed

### Comparing Root Finding Methods

When selecting a root-finding method for a specific problem, consider these key aspects:

#### Convergence Rate
- **Bisection**: Linear convergence, approximately halves the error each iteration
- **Newton-Raphson**: Quadratic convergence, roughly doubles the number of correct digits each iteration
- **Secant**: Superlinear convergence (order ≈ 1.618), faster than bisection but slower than Newton-Raphson

#### Requirements
- **Bisection**: Requires a bracket where the function changes sign
- **Newton-Raphson**: Requires the derivative of the function
- **Secant**: Requires two initial points but no derivative
- **Brentq**: Requires a bracket where the function changes sign

#### Reliability
- **Bisection**: Most reliable, guaranteed to converge if initial conditions are valid
- **Newton-Raphson**: Can fail if derivative approaches zero or initial guess is poor
- **Secant**: More reliable than Newton-Raphson but less than bisection
- **Brentq**: Very reliable, combines advantages of bisection and inverse quadratic interpolation

#### Efficiency
- Newton-Raphson typically requires the fewest iterations but has higher computational cost per iteration
- Brentq offers an excellent balance of reliability and efficiency
- Secant method is a good compromise when derivatives are unavailable

When comparing methods visually, error plots on semi-log scales show the different convergence rates clearly, with Newton-Raphson showing the steepest descent.

## 7. Best Practices for Root Finding

### Key Considerations

- **Problem Analysis**:
  - Understand the function's behavior before attempting to find roots
  - Identify possible regions where roots might exist
  - Determine if the function has multiple roots

- **Method Selection**:
  - Choose bracketing methods (bisection, Brent) for reliability
  - Use open methods (Newton-Raphson, secant) for speed when good initial guesses are available
  - Consider SciPy implementations for efficiency and robustness

- **Initial Guess Selection**:
  - For bracketing methods, ensure f(a) and f(b) have opposite signs
  - For open methods, choose initial points as close as possible to the expected root
  - Plot the function to visually identify good starting points

- **Convergence Criteria**:
  - Use both relative error and residual for stopping conditions
  - Set appropriate tolerance levels for your application
  - Implement maximum iteration limits as a safety measure

- **Error Handling**:
  - Check for divergence conditions
  - Handle special cases (e.g., derivative = 0 in Newton-Raphson)
  - Implement fallback methods when primary methods fail

### Best Practices for Robust Root Finding

When implementing root-finding algorithms in engineering applications, consider these best practices:

#### 1. Strategic Method Selection

Create a hierarchical approach to method selection:
- If a valid bracketing interval is available, use Brentq method (combines reliability of bisection with efficiency of inverse quadratic interpolation)
- If derivative is available and a reasonable initial guess exists, use Newton-Raphson
- If two initial guesses are available but no derivative, use the secant method
- Fall back to more reliable methods if the preferred method fails

#### 2. Comprehensive Error Handling

- Check for invalid inputs (e.g., non-bracketing intervals)
- Handle special cases like zero derivatives in Newton-Raphson
- Implement fallback strategies when primary methods fail
- Provide informative error messages about what was tried and why it failed

#### 3. Multiple Starting Points

For complex functions or when roots might be missed:
- Try multiple initial guesses in different regions
- Combine global search methods with local refinement
- Use bracketing methods to locate regions containing roots
- Implement grid searches if necessary

#### 4. Convergence Criteria

Use multiple stopping criteria:
- Approximate relative error below tolerance
- Function value (residual) below tolerance
- Maximum iterations reached

#### 5. Validation

Always validate the found root by:
- Checking that |f(root)| is sufficiently small
- Verifying the solution makes sense physically
- Testing with different methods to confirm consistency

## Tips for Week 3 Assignments

1. **Visualize the function** before attempting to find roots to understand its behavior.

2. **Test your code with simple examples** where you know the roots analytically (e.g., x² - 4 = 0 has roots at x = ±2).

3. **Implement error checking** for edge cases like division by zero or invalid intervals.

4. **Monitor convergence** by tracking both relative error and function value residual.

5. **Consider the tradeoffs** between different methods based on the specific problem requirements.

6. **Use SciPy's implementations** when appropriate rather than reimplementing algorithms from scratch.

7. **Document your functions** thoroughly with docstrings explaining parameters, return values, and behavior.

8. **Test with multiple initial guesses** to ensure your solution is robust.
