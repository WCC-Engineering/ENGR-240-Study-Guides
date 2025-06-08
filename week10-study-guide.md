## Overview

This study guide covers numerical methods for solving boundary value problems (BVPs). Week 10 introduces the fundamental differences between initial value problems and boundary value problems, and explores two primary numerical approaches: the shooting method and the finite difference method.

## 1. Introduction to Boundary Value Problems

### Key Concepts

- **Boundary Value Problem (BVP)**: An ODE where conditions are specified at two or more different points
- **Initial Value Problem (IVP)**: An ODE where all conditions are specified at a single starting point
- **Two-Point BVP**: Conditions given at two points: $y(a) = \alpha$, $y(b) = \beta$

### What Defines a Boundary Value Problem?

A second-order BVP has the general form:
$$\frac{d^2y}{dx^2} = f\left(x, y, \frac{dy}{dx}\right), \quad a \leq x \leq b$$

with boundary conditions at both ends.

### Types of Boundary Conditions

1. **Dirichlet**: $y(a) = \alpha$, $y(b) = \beta$
2. **Neumann**: $y'(a) = \alpha$, $y'(b) = \beta$  
3. **Mixed**: $y(a) = \alpha$, $y'(b) = \beta$

### Why Can't We Just Use IVP Solvers?

1. **Missing Initial Conditions**: We don't know all conditions at one point
2. **Global Constraints**: Solution must satisfy conditions at multiple points
3. **Iterative Nature**: Most BVP methods require iteration

### Engineering Applications

- **Heat Transfer**: Steady-state temperature: $\frac{d^2T}{dx^2} = 0$
- **Beam Deflection**: $\frac{d^2y}{dx^2} = \frac{M(x)}{EI}$
- **Fluid Flow**: Pressure distribution in pipes

### Resources

- **Davishahl Numerical Methods**
  - [Introduction to Boundary Value Problems](placeholder-for-video-link)

## 2. The Shooting Method

### Key Concepts

- **Convert BVP to IVPs**: "Shoot" from one boundary toward the other
- **Root-Finding**: Find initial condition that satisfies far boundary
- **Linear vs Nonlinear**: Different strategies for each case
- **Recommended Solver**: Use `solve_ivp` with `method='RK45'` for most cases

### Basic Algorithm

For BVP: $y'' = f(x,y,y')$, $y(a) = \alpha$, $y(b) = \beta$

1. **Guess** initial slope: $s = y'(a)$
2. **Solve IVP** from $a$ to $b$ using RK45
3. **Check** if $y(b) = \beta$  
4. **Adjust** and repeat

### Linear Case: Direct Solution

For linear BVPs: $y'' + p(x)y' + q(x)y = r(x)$

Use **linear interpolation**:
1. Solve with two different slopes $s_1, s_2$
2. Get end values $y_1(b), y_2(b)$
3. Correct slope: $s = s_1 + (s_2 - s_1) \frac{\beta - y_1(b)}{y_2(b) - y_1(b)}$

### Nonlinear Case: Root Finding

Define: $G(s) = y(b; s) - \beta$

Solve $G(s) = 0$ using:
- **fsolve** (recommended)
- Bisection method
- Secant method

{% raw %}
```python
def shooting_method_nonlinear(ode_func, a, b, alpha, beta, s_guess=0):
    """
    Solve nonlinear BVP using shooting method with root finding.
    
    Parameters:
    ode_func: function that returns [y', y''] given (x, y, y')
    a, b: domain boundaries
    alpha, beta: boundary values y(a)=alpha, y(b)=beta
    s_guess: initial guess for the slope y'(a)
    """
    
    def ode_system(x, y):
        """
        Convert 2nd order ODE to system of 1st order ODEs
        y[0] = y (function value)
        y[1] = y' (derivative value)
        """
        y1, y2 = y  # Extract function and derivative values
        derivatives = ode_func(x, y1, y2)  # Get [y', y''] from user function
        return derivatives
    
    def boundary_residual(s):
        """
        Function whose root gives the correct initial slope.
        
        Parameters:
        s: trial value for initial slope y'(a)
        
        Returns:
        Residual: y(b) - beta (should be zero for correct slope)
        """
        # Solve IVP with current slope guess
        sol = solve_ivp(ode_system, [a, b], [alpha, s], 
                       method='RK45', rtol=1e-8)
        
        # Return difference between computed y(b) and target beta
        return sol.y[0, -1] - beta  # y(b) - beta
    
    # Find the root (correct slope) using scipy's robust solver
    s_correct = fsolve(boundary_residual, s_guess)[0]
    
    # Solve one final time with the correct slope for full solution
    sol_final = solve_ivp(ode_system, [a, b], [alpha, s_correct], 
                         method='RK45', dense_output=True, rtol=1e-8)
    
    return sol_final, s_correct
```
{% endraw %}

### Strategies for Initial Guesses

1. **Linear interpolation**: $s = \frac{\beta - \alpha}{b - a}$
2. **Physical intuition**: Engineering judgment
3. **Bracketing**: Two guesses with opposite residual signs
4. **Multiple attempts**: Try different starting values

### Advantages/Disadvantages

**Advantages:**
- Uses robust IVP solvers (RK45)
- High accuracy possible
- Works for nonlinear problems

**Disadvantages:**
- Sensitive to initial guesses
- Instability for some problems
- Multiple IVP solutions required

### Resources

- **Davishahl Numerical Methods**
  - [Shooting Method for Linear BVPs](placeholder-for-video-link)
  - [Shooting Method for Nonlinear BVPs](placeholder-for-video-link)

## 3. Finite Difference Method (Linear BVPs)

### Key Concepts

- **Discretize domain**: Replace derivatives with finite differences
- **Matrix formulation**: Convert to linear system $\mathbf{Ay} = \mathbf{b}$
- **Direct solution**: Solve system for all unknowns simultaneously

### Mesh Generation

Create uniform mesh: $x_i = a + ih$, $h = \frac{b-a}{n+1}$
- Interior points: $x_1, x_2, \ldots, x_n$
- Boundary points: $x_0 = a$, $x_{n+1} = b$

### Finite Difference Approximations

**Second derivative (central)**: $\frac{d^2y}{dx^2}\bigg|_{x_i} \approx \frac{y_{i+1} - 2y_i + y_{i-1}}{h^2}$

**First derivative (central)**: $\frac{dy}{dx}\bigg|_{x_i} \approx \frac{y_{i+1} - y_{i-1}}{2h}$

### General Linear BVP

For: $y'' + p(x)y' + q(x)y = r(x)$

At each interior point:
$$\frac{y_{i+1} - 2y_i + y_{i-1}}{h^2} + p_i\frac{y_{i+1} - y_{i-1}}{2h} + q_i y_i = r_i$$

Rearranging:
$$\left(1 + \frac{hp_i}{2}\right)y_{i+1} + \left(-2 + h^2q_i\right)y_i + \left(1 - \frac{hp_i}{2}\right)y_{i-1} = h^2r_i$$

### Special Case: No First Derivative Term

For the common case $y'' + q(x)y = r(x)$ (no $y'$ term), the finite difference equation simplifies to:

$$\frac{y_{i+1} - 2y_i + y_{i-1}}{h^2} + q_i y_i = r_i$$

Rearranging: $y_{i-1} + (-2 + h^2 q_i)y_i + y_{i+1} = h^2 r_i$

{% raw %}
```python
def setup_tridiagonal_matrix_simple(q_func, r_func, a, b, alpha, beta, n):
    """
    Set up tridiagonal matrix for y'' + q(x)*y = r(x)
    with Dirichlet BCs: y(a) = alpha, y(b) = beta
    """
    # Create mesh
    h = (b - a) / (n + 1)
    x_interior = np.linspace(a + h, b - h, n)  # Interior points only
    
    # Evaluate functions at interior points
    q_vals = np.array([q_func(xi) for xi in x_interior])
    r_vals = np.array([r_func(xi) for xi in x_interior])
    
    # Set up tridiagonal matrix coefficients
    # Lower diagonal: all ones (coefficient of y_{i-1})
    lower_diag = np.ones(n-1)
    
    # Main diagonal: -2 + h^2 * q(x_i)
    main_diag = -2 + h**2 * q_vals
    
    # Upper diagonal: all ones (coefficient of y_{i+1})
    upper_diag = np.ones(n-1)
    
    # Create sparse tridiagonal matrix
    A = diags([lower_diag, main_diag, upper_diag], [-1, 0, 1], 
              shape=(n, n), format='csr')
    
    # Set up right-hand side vector
    b_vec = h**2 * r_vals
    
    # Apply Dirichlet boundary conditions
    # Modify first equation: subtract alpha from RHS
    b_vec[0] -= alpha
    
    # Modify last equation: subtract beta from RHS  
    b_vec[-1] -= beta
    
    return A, b_vec, x_interior
```
{% endraw %}

### Matrix System

Results in tridiagonal system:

{% raw %}
```python
# Main diagonal
main_diag = -2 + h**2 * q_vals

# Off-diagonals  
upper_diag = 1 + h * p_vals[:-1] / 2
lower_diag = 1 - h * p_vals[1:] / 2

# Create matrix
A = diags([lower_diag, main_diag, upper_diag], [-1, 0, 1])
```
{% endraw %}

### Implementing Boundary Conditions

For Dirichlet BCs: $y(a) = \alpha$, $y(b) = \beta$

Modify first/last equations:
- Move known boundary values to RHS
- First equation: subtract $\left(1 - \frac{hp_1}{2}\right)\alpha$
- Last equation: subtract $\left(1 + \frac{hp_n}{2}\right)\beta$

### Convergence Behavior

**Second-order convergence**: Error = $O(h^2)$
- Halving mesh size reduces error by factor of 4
- Very predictable convergence

### Advantages/Disadvantages

**Advantages:**
- Systematic approach
- Sparse matrices
- Predictable convergence
- Direct solution

**Disadvantages:**
- Linear problems only
- Uniform mesh typically required
- Special handling for different BC types

### Resources

- **Davishahl Numerical Methods**
  - [Finite Difference Method for BVPs](placeholder-for-video-link)

## 4. Finite Difference Method for Nonlinear BVPs (Conceptual)

### Key Concepts

For nonlinear BVP: $y'' = f(x, y, y')$

Discretization yields **system of nonlinear algebraic equations**:
$$F_i(y_1, y_2, \ldots, y_n) = 0, \quad i = 1, 2, \ldots, n$$

### Solution Strategy (High-Level)

1. **Newton's Method** for systems:
   $$\mathbf{y}^{(k+1)} = \mathbf{y}^{(k)} - \mathbf{J}^{-1}(\mathbf{y}^{(k)}) \mathbf{F}(\mathbf{y}^{(k)})$$

2. **Jacobian Matrix**: $J_{ij} = \frac{\partial F_i}{\partial y_j}$ (tridiagonal structure)

3. **Iteration Process**:
   - Start with initial guess
   - Solve linear system at each iteration
   - Update solution until convergence

### Why Not Covered in Detail

Requires understanding of:
- Systems of nonlinear algebraic equations
- Advanced linear algebra techniques
- Jacobian computation methods
- Convergence theory for Newton's method

### Alternative for Nonlinear Problems

**Shooting method** often more straightforward for nonlinear BVPs than finite difference approach.

## Tips for Week 10 Assignments

1. **Shooting Method**:
   - Start with simple initial guesses
   - Use physical intuition for estimates
   - Try multiple starting values for nonlinear problems
   - Check convergence of root-finding algorithm
   - Use `method='RK45'` in `solve_ivp` for reliability

2. **Finite Difference Method**:
   - Start with coarse mesh, refine gradually
   - Verify second-order convergence
   - Use sparse matrix routines for efficiency
   - Check boundary condition implementation

3. **Method Selection**:
   - **Linear problems**: Either method works well
   - **Nonlinear problems**: Shooting often easier
   - **Complex geometry**: Finite difference more flexible
   - **High accuracy needed**: Both can achieve machine precision

4. **Debugging**:
   - Plot solutions to check reasonableness
   - Compare with analytical solutions when available
   - Check error convergence rates
   - Verify boundary conditions are satisfied

5. **Common Pitfalls**:
   - Poor initial guesses in shooting method
   - Incorrect boundary condition implementation
   - Not checking convergence
   - Using too coarse mesh in finite difference

