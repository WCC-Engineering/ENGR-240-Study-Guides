## Overview

This study guide covers numerical methods for solving ordinary differential equations (ODEs), a critical topic in engineering computations. Week 8 introduces the fundamental concepts of differential equations, various time-stepping methods for solving initial value problems, analysis of solution accuracy and stability, and techniques for solving systems of ODEs.

## 1. Introduction to Differential Equations

### Key Concepts

- **Differential Equation**: An equation containing derivatives of an unknown function
- **Ordinary Differential Equation (ODE)**: Contains derivatives with respect to a single independent variable
- **Order**: The highest derivative that appears in the equation
- **Initial Value Problem (IVP)**: An ODE with conditions specified at the starting point
- **Boundary Value Problem (BVP)**: An ODE with conditions specified at multiple points
- **Analytical vs. Numerical Solutions**:
  - Analytical: Exact mathematical expressions
  - Numerical: Approximations using computational methods

### Classification of Differential Equations

#### By Type
- **Ordinary Differential Equations (ODEs)**: Involve derivatives with respect to a single variable
  
  Example: $\frac{dy}{dt} = -ky$ (First-order ODE describing exponential decay)
  
- **Partial Differential Equations (PDEs)**: Involve partial derivatives with respect to multiple variables
  
  Example: $\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} = 0$ (Laplace's equation)

#### By Order
- **First-Order**: Contains only first derivatives
  
  Example: $\frac{dy}{dt} = f(t, y)$
  
- **Second-Order**: Contains second derivatives
  
  Example: $\frac{d^2y}{dt^2} + \omega^2 y = 0$ (Simple harmonic oscillator)
  
- **Higher-Order**: Contains derivatives of order three or higher

#### By Linearity
- **Linear**: The dependent variable and its derivatives appear only to the first power and are not multiplied together
  
  Example: $\frac{d^2y}{dt^2} + a\frac{dy}{dt} + by = f(t)$
  
- **Nonlinear**: Contains products of the dependent variable with itself or its derivatives, or nonlinear functions of these
  
  Example: $\frac{dy}{dt} = y^2 - t$

### Initial Value Problems (IVPs)

An Initial Value Problem consists of:
1. A differential equation: $\frac{dy}{dt} = f(t, y)$
2. Initial condition(s): $y(t_0) = y_0$

For an nth-order ODE, n initial conditions are required, typically the values of the function and its first n-1 derivatives at the starting point.

**Example**: First-order IVP
- Differential equation: $\frac{dy}{dt} = -ky$
- Initial condition: $y(0) = y_0$
- Analytical solution: $y(t) = y_0 e^{-kt}$

**Example**: Second-order IVP
- Differential equation: $\frac{d^2y}{dt^2} + \omega^2 y = 0$
- Initial conditions: $y(0) = y_0$, $\frac{dy}{dt}(0) = v_0$
- Analytical solution: $y(t) = y_0 \cos(\omega t) + \frac{v_0}{\omega} \sin(\omega t)$

### Boundary Value Problems (BVPs)

A Boundary Value Problem consists of:
1. A differential equation: $\frac{d^2y}{dx^2} = f(x, y, \frac{dy}{dx})$
2. Boundary conditions specified at different points: $y(a) = \alpha$, $y(b) = \beta$

BVPs typically arise in steady-state problems and require different solution techniques than IVPs.

**Example**: Second-order BVP
- Differential equation: $\frac{d^2y}{dx^2} + q(x)y = f(x)$
- Boundary conditions: $y(0) = 0$, $y(L) = 0$
- This could represent the deflection of a beam with fixed ends

### Common Engineering Applications

Differential equations are ubiquitous in engineering:

1. **Mechanical Systems**:
   - Spring-mass-damper: $m\frac{d^2x}{dt^2} + c\frac{dx}{dt} + kx = F(t)$
   - Structural vibrations and dynamics

2. **Electrical Systems**:
   - RLC circuits: $L\frac{d^2I}{dt^2} + R\frac{dI}{dt} + \frac{1}{C}I = V(t)$
   - Signal processing and control systems

3. **Thermal Systems**:
   - Heat conduction: $\frac{\partial T}{\partial t} = \alpha \nabla^2 T$
   - Cooling and heating processes

4. **Fluid Dynamics**:
   - Navier-Stokes equations
   - Pipe flow and hydraulic systems

5. **Chemical Reactions**:
   - Reaction kinetics: $\frac{dC}{dt} = -kC^n$
   - Batch and continuous reactors

### Resources

- **Berkeley Numerical Methods**
  - [Ordinary Differential Equations](https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter22.00-Ordinary-Differential-Equations.html)

## 2. Numerical Methods for Initial Value Problems

### Key Concepts

- **Time-Stepping Methods**: Algorithms that advance a solution from time t to t+h
- **One-Step Methods**: Use only the current state to calculate the next state
- **Multi-Step Methods**: Use information from multiple previous steps
- **Explicit Methods**: Calculate the next state directly from current values
- **Implicit Methods**: Require solving an equation involving the next state
- **Local Truncation Error**: Error in a single step
- **Global Truncation Error**: Accumulated error over multiple steps
- **Stability**: Behavior of numerical solutions as step size and integration time vary

### Common Terminology

- **Step Size (h)**: The interval between successive solution points
- **Approximation**: The calculated value $y_n$ at time $t_n$
- **Exact Solution**: The true value $y(t_n)$ at time $t_n$
- **Converged Solution**: Solution obtained as step size approaches zero
- **Order of Accuracy**: Rate at which error decreases as step size decreases

### Mathematical Framework

For the first-order ODE:
$$\frac{dy}{dt} = f(t, y), \quad y(t_0) = y_0$$

We seek to find approximate values $y_1, y_2, ..., y_n$ corresponding to times $t_1, t_2, ..., t_n$.

The general form of a one-step method can be written as:
$$y_{n+1} = y_n + h \cdot \phi(t_n, y_n, h, f)$$

where $\phi$ is the increment function that depends on the specific method.

### Resources

- **Davishahl Numerical Methods Videos**
  - [Intro to ODEs](INSERT VIDEO LINK HERE)
  - [Euler's Method](INSERT VIDEO LINK HERE)
  - [Improved Euler's Method](INSERT VIDEO LINK HERE)
