# Week 4 Study Guide: Linear Systems of Equations

## Overview

This study guide covers numerical methods for solving systems of linear equations, a fundamental problem in engineering computation. The activities, worksheets, and programming assignments for Week 4 introduce matrix notation, Gauss elimination, LU factorization, built-in NumPy/SciPy functions for solving linear systems, and the concept of condition number.

## 1. Matrix Notation for Linear Systems

### Key Concepts

- **Systems of Linear Equations**: A set of equations where each variable appears at most to the first power
- **Matrix Notation**: A compact way to represent systems using matrices and vectors
- **Standard Form**: Ax = b, where A is the coefficient matrix, x is the vector of unknowns, and b is the right-hand side vector
- **Solution Types**:
  - Unique solution: Exactly one value for each unknown
  - Infinite solutions: Underdetermined system
  - No solution: Overdetermined or inconsistent system

### Matrix Representation of Linear Systems

A system of linear equations:

$$
\begin{align}
a_{11}x_1 + a_{12}x_2 + \ldots + a_{1n}x_n &= b_1\\
a_{21}x_1 + a_{22}x_2 + \ldots + a_{2n}x_n &= b_2\\
\vdots\\
a_{m1}x_1 + a_{m2}x_2 + \ldots + a_{mn}x_n &= b_m
\end{align}
$$

Can be represented in matrix form as:

$$
\begin{bmatrix}
a_{11} & a_{12} & \cdots & a_{1n} \\
a_{21} & a_{22} & \cdots & a_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
a_{m1} & a_{m2} & \cdots & a_{mn}
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
\vdots \\
x_n
\end{bmatrix}
=
\begin{bmatrix}
b_1 \\
b_2 \\
\vdots \\
b_m
\end{bmatrix}
$$

Or simply:

$$\textbf{A}\textbf{x} = \textbf{b}$$

### Properties of Linear Systems

- **Square System**: Number of equations equals number of unknowns (m = n)
- **Overdetermined System**: More equations than unknowns (m > n)
- **Underdetermined System**: Fewer equations than unknowns (m < n)
- **Homogeneous System**: Right-hand side is zero (b = 0)
- **Non-homogeneous System**: Right-hand side is non-zero (b ≠ 0)

### Existence and Uniqueness of Solutions

- A square system (n equations, n unknowns) has a unique solution if and only if the determinant of the coefficient matrix is non-zero (det(A) ≠ 0)
- An equivalent condition is that the matrix A must be invertible (non-singular)
- When det(A) = 0, the matrix is singular and either:
  - The system has no solution (inconsistent)
  - The system has infinitely many solutions (dependent equations)

### Engineering Applications of Linear Systems

- Structural analysis (force distributions in trusses)
- Circuit analysis (Kirchhoff's laws)
- Heat transfer calculations
- Material balance equations in chemical engineering
- Traffic flow modeling
- Mechanical equilibrium problems

### Resources

- **Sullivan Numerical Methods Textbook**:
  - [Introduction to Linear Systems](https://numericalmethodssullivan.github.io/ch-linear.html)

## 2. Gauss Elimination

### Key Concepts

- **Principle**: Systematic elimination of variables to convert the system to an upper triangular form
- **Forward Elimination**: The process of obtaining an upper triangular matrix
- **Back Substitution**: The process of finding solutions starting from the last equation
- **Pivoting**: Techniques to improve numerical stability
  - Partial pivoting: Row exchanges
  - Complete pivoting: Row and column exchanges
- **Elementary Row Operations**:
  - Multiply a row by a non-zero constant
  - Add a multiple of one row to another
  - Interchange two rows

### Gauss Elimination Algorithm

1. **Forward Elimination Phase**:
   - For each pivot position (diagonal element):
     - Select the pivot (optionally implement pivoting strategies)
     - Eliminate all elements below the pivot by subtracting multiples of the pivot row

2. **Back Substitution Phase**:
   - Start with the last equation (now with only one unknown)
   - Solve for that unknown
   - Substitute the value into the previous equation and solve
   - Continue until all unknowns are found

### Mathematical Representation

For a 3×3 system:

$$
\begin{bmatrix}
a_{11} & a_{12} & a_{13} \\
a_{21} & a_{22} & a_{23} \\
a_{31} & a_{32} & a_{33}
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix}
=
\begin{bmatrix}
b_1 \\
b_2 \\
b_3
\end{bmatrix}
$$

**Forward Elimination**:

First, eliminate $x_1$ from the second and third equations:
- Compute multipliers: $m_{21} = a_{21}/a_{11}$ and $m_{31} = a_{31}/a_{11}$
- Update second row: $[a_{21}~a_{22}~a_{23}~b_2] = [a_{21}~a_{22}~a_{23}~b_2] - m_{21}[a_{11}~a_{12}~a_{13}~b_1]$
- Update third row: $[a_{31}~a_{32}~a_{33}~b_3] = [a_{31}~a_{32}~a_{33}~b_3] - m_{31}[a_{11}~a_{12}~a_{13}~b_1]$

Then, eliminate $x_2$ from the third equation:
- Compute multiplier: $m_{32} = a_{32}/a_{22}$
- Update third row: $[a_{32}~a_{33}~b_3] = [a_{32}~a_{33}~b_3] - m_{32}[a_{22}~a_{23}~b_2]$

Resulting in upper triangular form:

$$
\begin{bmatrix}
a_{11} & a_{12} & a_{13} \\
0 & a_{22}' & a_{23}' \\
0 & 0 & a_{33}'
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix}
=
\begin{bmatrix}
b_1 \\
b_2' \\
b_3'
\end{bmatrix}
$$

**Back Substitution**:

Solve equations in reverse order:
- $x_3 = b_3'/a_{33}'$
- $x_2 = (b_2' - a_{23}'x_3)/a_{22}'$
- $x_1 = (b_1 - a_{12}x_2 - a_{13}x_3)/a_{11}$

### Operation Counting

For an n×n system, Gauss elimination requires:
- Forward elimination: approximately $\frac{n^3}{3} + \frac{n^2}{2} - \frac{5n}{6}$ operations
- Back substitution: approximately $n^2$ operations
- Total operations: $O(n^3)$ (cubic complexity)

This operation count demonstrates why solving large linear systems is computationally intensive.

### Pivoting for Numerical Stability

When the pivot element is zero or very small, numerical errors can be magnified. Pivoting strategies include:

**Partial Pivoting**:
- Search for the largest absolute value in the current column below the diagonal
- Swap rows to place this element in the pivot position
- Reduces numerical errors and prevents division by zero

**Complete Pivoting**:
- Search for the largest absolute value in the remaining submatrix
- Swap both rows and columns to place this element in the pivot position
- Maximum numerical stability but more computationally expensive

### Resources

- **Video Lectures**:
  - [Gauss Elimination Introduction](#) <!-- Placeholder for video link -->
  - [Gauss Elimination with Pivoting](#) <!-- Placeholder for video link -->

## 3. LU Factorization

### Key Concepts

- **Principle**: Decompose the coefficient matrix A into a product of lower (L) and upper (U) triangular matrices
- **Form**: A = LU
- **Advantages**:
  - Efficient for solving multiple systems with the same coefficient matrix but different right-hand sides
  - Simpler operations for computing determinants and inverses
  - Reduced computational cost for repeated solutions
- **Variants**:
  - LU decomposition without pivoting
  - LU decomposition with partial pivoting (LUP factorization)
  - Cholesky decomposition for symmetric positive-definite matrices

### LU Factorization Process

The LU factorization decomposes a matrix A into:

$$A = LU$$

Where:
- L is a lower triangular matrix with ones on the diagonal
- U is an upper triangular matrix

For a 3×3 system:

$$
\begin{bmatrix}
a_{11} & a_{12} & a_{13} \\
a_{21} & a_{22} & a_{23} \\
a_{31} & a_{32} & a_{33}
\end{bmatrix}
=
\begin{bmatrix}
1 & 0 & 0 \\
l_{21} & 1 & 0 \\
l_{31} & l_{32} & 1
\end{bmatrix}
\begin{bmatrix}
u_{11} & u_{12} & u_{13} \\
0 & u_{22} & u_{23} \\
0 & 0 & u_{33}
\end{bmatrix}
$$

### Computing LU Factorization

The LU factorization can be computed using several methods:

1. **Doolittle Algorithm**:
   - U is the same as the upper triangular matrix obtained from Gauss elimination
   - L contains the multipliers used during elimination

2. **Crout Algorithm**:
   - Different algorithm leading to the same result
   - More numerically stable in some cases

3. **LU with Partial Pivoting**:
   - Adds row interchanges with a permutation matrix P
   - Results in PA = LU

### Solving Systems with LU Factorization

Once A = LU factorization is available, solving Ax = b becomes:
1. LUx = b
2. Let Ux = y, then Ly = b

The solution process:
1. Forward substitution: Solve Ly = b for y
2. Back substitution: Solve Ux = y for x

This is particularly efficient when solving multiple systems with the same A but different b:
- Factorize A = LU once: $O(n^3)$ operations
- Solve Ly = b and Ux = y for each new b: $O(n^2)$ operations per system

### Advantages of LU Factorization

1. **Computational Efficiency**: For multiple right-hand sides (b vectors)
   - Single factorization: $O(n^3)$ operations
   - Each additional solution: $O(n^2)$ operations

2. **Matrix Properties**:
   - Determinant calculation: det(A) = det(L) × det(U) = product of diagonal elements of U
   - Inverse calculation: More efficient through LU factorization

3. **Stability Analysis**:
   - Error analysis is more straightforward
   - Pivoting strategies can be incorporated

### Resources

- **Video Lectures**:
  - [LU Factorization Explained](#) <!-- Placeholder for video link -->
  - [Using LU Factorization to Solve Systems](#) <!-- Placeholder for video link -->

## 4. NumPy/SciPy Functions for Linear Systems

### Key Concepts

- **NumPy Functions**: np.linalg module provides core linear algebra operations
- **SciPy Functions**: scipy.linalg module offers additional specialized functions
- **Key Operations**:
  - Direct solving of linear systems
  - Computing matrix inverse
  - LU decomposition
  - Analysis of matrix properties (determinant, rank, etc.)
- **Advantages**: Optimized implementations, vectorized operations, numerical stability

### NumPy Functions for Linear Systems

NumPy's `numpy.linalg` module provides essential functions for linear algebra operations:

```python
import numpy as np

# Define coefficient matrix and right-hand side
A = np.array([[4, 2, 3], [3, 1, -2], [2, -3, 1]])
b = np.array([7, -1, 2])

# Solve the system using numpy.linalg.solve
x = np.linalg.solve(A, b)
print("Solution:", x)

# Check the solution
print("Verification:", np.allclose(np.dot(A, x), b))

# Compute the inverse of A
A_inv = np.linalg.inv(A)
print("Matrix inverse:\n", A_inv)

# Compute the determinant
det_A = np.linalg.det(A)
print("Determinant:", det_A)

# Compute the rank
rank_A = np.linalg.matrix_rank(A)
print("Rank:", rank_A)
```

### SciPy Functions for Linear Systems

SciPy's `scipy.linalg` module extends NumPy's capabilities with more specialized functions:

```python
import numpy as np
import scipy.linalg as spla

# Define coefficient matrix and right-hand side
A = np.array([[4, 2, 3], [3, 1, -2], [2, -3, 1]])
b = np.array([7, -1, 2])

# LU decomposition with SciPy
P, L, U = spla.lu(A)
print("L matrix:\n", L)
print("U matrix:\n", U)
print("P matrix:\n", P)

# Verify A = P.T @ L @ U
print("Verification:", np.allclose(A, P.T @ L @ U))

# Solve using LU factorization
# Forward substitution
y = spla.solve_triangular(L, P @ b, lower=True)
# Back substitution
x = spla.solve_triangular(U, y, lower=False)
print("Solution via LU:", x)

# More efficient diagonal solving for triangular systems
y = spla.solve_triangular(L, P @ b, lower=True, unit_diagonal=True)
```