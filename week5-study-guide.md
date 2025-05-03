# Week 5 Study Guide: Curve Fitting with Regression Analysis

## Overview

This study guide covers curve fitting and regression analysis techniques used in engineering and scientific computing. The activities, worksheets, and programming assignments for Week 5 introduce several approaches to fitting mathematical models to experimental data, including linear regression, nonlinear regression, and general linear least squares.

## 1. Introduction to Curve Fitting

### Key Concepts

- **Definition**: Curve fitting is the process of finding a mathematical function that best represents a set of data points
- **Purpose**:
  - Interpolation: Estimating values between known data points
  - Extrapolation: Predicting values beyond the range of known data
  - Parameter estimation: Determining coefficients for theoretical models
  - Data smoothing: Reducing noise in experimental measurements
- **Common Engineering Applications**:
  - Material property characterization
  - Reaction kinetics modeling
  - Calibration of sensors and instruments
  - Analysis of experimental data
  - System identification and modeling

### Mathematical Representation

For a set of data points $(x_i, y_i)$ for $i = 1, 2, ..., n$, we seek a function $f(x)$ with parameters $p_1, p_2, ..., p_m$ such that:

$$f(x_i; p_1, p_2, ..., p_m) \approx y_i \text{ for all } i$$

The goal is to determine the parameter values that minimize some measure of the difference between the observed $y_i$ values and the predicted $f(x_i)$ values.

### Types of Curve Fitting

1. **Linear Regression**:
   - Fits data to a linear function $f(x) = ax + b$
   - Parameters have a direct interpretation (slope and intercept)
   - Can be solved analytically using the method of least squares

2. **Polynomial Regression**:
   - Fits data to a polynomial function $f(x) = a_0 + a_1x + a_2x^2 + ... + a_mx^m$
   - More flexible than linear regression but may overfit with high-order terms
   - Still linear in parameters, so can use linear least squares methods

3. **Nonlinear Regression**:
   - Fits data to nonlinear functions like exponential, logarithmic, or power law models
   - Often requires iterative numerical optimization techniques
   - May have parameters with direct physical interpretation

4. **Linearized Nonlinear Models**:
   - Transforms nonlinear models to create linear relationships
   - Allows use of linear least squares methods for nonlinear forms
   - May introduce bias in the error distribution

### Resources

- **Davishahl Numerical Methods Videos**:
  - [Linear Regression](https://youtu.be/qo6-Dp_smIg?si=DBpo4WbImf5VdX2L)

- **Berkeley Numerical Methods**:
  - [Least-Squares Regression](https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter16.00-Least-Squares-Regression.html)

## 2. Method of Least Squares

### Key Concepts

- **Principle**: Minimizes the sum of squared differences between observed and predicted values
- **Objective Function**: Minimize $\sum_{i=1}^{n} [y_i - f(x_i)]^2$
- **Linear Least Squares**: Closed-form solution exists for linear models
- **Nonlinear Least Squares**: Requires iterative methods to minimize the objective function
- **Weighted Least Squares**: Applies different weights to different data points based on their reliability

### Linear Least Squares Derivation

For a linear model $f(x) = ax + b$, we minimize:

$$S = \sum_{i=1}^{n} [y_i - (ax_i + b)]^2$$

Taking partial derivatives with respect to $a$ and $b$ and setting them to zero:

$$\frac{\partial S}{\partial a} = -2\sum_{i=1}^{n} x_i[y_i - (ax_i + b)] = 0$$

$$\frac{\partial S}{\partial b} = -2\sum_{i=1}^{n} [y_i - (ax_i + b)] = 0$$

Rearranging:

$$a\sum_{i=1}^{n} x_i^2 + b\sum_{i=1}^{n} x_i = \sum_{i=1}^{n} x_iy_i$$

$$a\sum_{i=1}^{n} x_i + bn = \sum_{i=1}^{n} y_i$$

Solving these "normal equations" yields the optimal values for $a$ and $b$.

### Matrix Formulation for Linear Least Squares

For a general linear model $f(x) = \beta_1 Z_1(x) + \beta_2 Z_2(x) + ... + \beta_m Z_m(x)$ where $Z_j(x)$ are basis functions:

1. Form the design matrix $\mathbf{Z}$ where $Z_{ij} = Z_j(x_i)$
2. Form the vector of observations $\mathbf{y} = [y_1, y_2, ..., y_n]^T$
3. The least squares solution is given by $\boldsymbol{\beta} = (\mathbf{Z}^T\mathbf{Z})^{-1}\mathbf{Z}^T\mathbf{y}$

In NumPy, this can be computed using `np.linalg.lstsq(Z, y)`.

### Resources

- **Davishahl Numerical Methods Videos**:

  - [Linear Regression](https://youtu.be/qo6-Dp_smIg?si=DBpo4WbImf5VdX2L)


## 3. Linear Regression

### Key Concepts

- **Simple Linear Regression**: Fits data to a straight line $y = mx + b$
- **Multiple Linear Regression**: Fits data to a function of multiple variables
- **Polynomial Regression**: Fits data to a polynomial function
- **Advantages**:
  - Analytically solvable (no iterative methods required)
  - Computationally efficient
  - Easy to interpret results
- **Limitations**:
  - Only applicable to linear (or linearized) relationships
  - May perform poorly on inherently nonlinear data

### Implementation in Python

{% raw %}
```python
import numpy as np
import matplotlib.pyplot as plt

# Sample data
x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
y = np.array([1.1, 1.9, 3.2, 3.8, 5.1])

# Method 1: Using numpy.polyfit
m, b = np.polyfit(x, y, 1)  # 1 indicates linear (first-order polynomial)
print(f"Slope: {m}, Intercept: {b}")

# Method 2: Using numpy.linalg.lstsq
A = np.vstack([x, np.ones(len(x))]).T  # Design matrix
m_lstsq, b_lstsq = np.linalg.lstsq(A, y, rcond=None)[0]
print(f"Slope: {m_lstsq}, Intercept: {b_lstsq}")

# Plot the data and fitted line
plt.figure(figsize=(8, 6))
plt.scatter(x, y, color='blue', label='Data points')
plt.plot(x, m*x + b, color='red', label=f'Fitted line: y = {m:.2f}x + {b:.2f}')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.show()
```
{% endraw %}

### Applications and Extensions

1. **Polynomial Regression**:
   - Uses higher-order polynomials: $y = a_0 + a_1x + a_2x^2 + ... + a_nx^n$
   - Implemented with `np.polyfit(x, y, n)` where `n` is the polynomial degree
   - Higher flexibility but risk of overfitting with high-order terms

2. **Multiple Linear Regression**:
   - Handles multiple independent variables: $y = \beta_0 + \beta_1x_1 + \beta_2x_2 + ... + \beta_mx_m$
   - Design matrix has a column for each variable plus a column of ones for the intercept
   - Useful for modeling complex relationships with multiple factors

## 4. Linearizing Nonlinear Models

### Key Concepts

- **Principle**: Transform nonlinear equations into linear form through algebraic manipulation
- **Common Transformations**:
  - Logarithmic: For exponential, power law, and multiplicative models
  - Reciprocal: For certain rational functions
  - Polynomial: For polynomial-like functions
- **Advantages**:
  - Uses simple linear regression techniques for nonlinear relationships
  - Often computationally less intensive than direct nonlinear regression
- **Limitations**:
  - Transformation may distort error distribution
  - May not be as accurate as direct nonlinear regression
  - Not all nonlinear models can be linearized

### Common Linearizable Models

1. **Exponential Model**:
   - Original form: $y = ae^{bx}$
   - Taking natural logarithm: $\ln(y) = \ln(a) + bx$
   - Linear form: $Y = A + bX$ where $Y = \ln(y)$ and $A = \ln(a)$

2. **Power Law Model**:
   - Original form: $y = ax^b$
   - Taking natural logarithm: $\ln(y) = \ln(a) + b\ln(x)$
   - Linear form: $Y = A + bX$ where $Y = \ln(y)$, $X = \ln(x)$, and $A = \ln(a)$

3. **Saturation Model (Michaelis-Menten)**:
   - Original form: $y = \frac{ax}{b+x}$
   - Taking reciprocal: $\frac{1}{y} = \frac{b+x}{ax} = \frac{b}{ax} + \frac{1}{a}$
   - Linear form: $Y = mX + c$ where $Y = \frac{1}{y}$, $X = \frac{1}{x}$, $m = \frac{b}{a}$, and $c = \frac{1}{a}$

### Example: Bacterial Growth Model

{% raw %}
```python
import numpy as np
import matplotlib.pyplot as plt

# Bacterial growth data: substrate concentration (c) and growth rate (k)
c = np.array([0.5, 0.8, 1.5, 2.5, 4.0])  # Substrate concentration (mg/L)
k = np.array([1.0, 2.5, 5.1, 7.3, 9.1])  # Growth rate (number/day)

# Model: k = kmax * c^2 / (cs + c^2)
# Linearized form: 1/k = (cs/kmax) * (1/c^2) + (1/kmax)

# Transform data
x_transformed = 1 / c**2  # 1/c^2
y_transformed = 1 / k     # 1/k

# Linear regression on transformed data
m, b = np.polyfit(x_transformed, y_transformed, 1)

# Calculate original model parameters
kmax = 1 / b
cs = m * kmax

print(f"Linear regression results: y = {m:.4f}x + {b:.4f}")
print(f"Model parameters: kmax = {kmax:.4f}, cs = {cs:.4f}")

# Define the model function
def kmodel(c, kmax, cs):
    """Growth rate model: k = kmax * c^2 / (cs + c^2)"""
    return kmax * c**2 / (cs + c**2)

# Generate smooth curve
c_model = np.linspace(0, 5, 100)
k_model = kmodel(c_model, kmax, cs)

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(c, k, 'ko', label='Experimental data')
plt.plot(c_model, k_model, 'r-', label='Fitted model')
plt.xlabel('Substrate Concentration (mg/L)')
plt.ylabel('Growth Rate (number/day)')
plt.legend()
plt.grid(True)
plt.show()
```
{% endraw %}

### Resources

- **Davishahl Numerical Methods Videos**:
  - [Linear Regression with Nonlinear Models](https://youtu.be/FdhCnHXmoi8?si=itNgU6zfabQzRj0O)

## 5. Direct Nonlinear Regression

### Key Concepts

- **Principle**: Directly fit nonlinear models without transformation
- **Methods**:
  - Iterative numerical optimization techniques
  - Minimizes sum of squared residuals (or other objective functions)
- **Common Algorithms**:
  - Levenberg-Marquardt
  - Gauss-Newton
  - Nelder-Mead simplex method
  - Gradient descent
- **Advantages**:
  - Can fit any nonlinear model
  - Preserves original error structure
  - Often more accurate than linearization
- **Limitations**:
  - Computationally more intensive
  - May converge to local rather than global optima
  - Requires good initial parameter estimates

### Implementation using SciPy

The primary SciPy functions for nonlinear regression are:

1. **`scipy.optimize.curve_fit`**:
   - Purpose-built for curve fitting
   - Based on Levenberg-Marquardt algorithm
   - Returns optimal parameters and covariance matrix

2. **`scipy.optimize.minimize`**:
   - General-purpose optimization function
   - Multiple algorithms available
   - More flexible but requires manual setup

### Example using curve_fit

{% raw %}
```python
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

# Sample data
x = np.array([0, 1, 2, 3, 4, 5])
y = np.array([1.0, 1.8, 3.5, 7.4, 14.6, 29.5])

# Define nonlinear model function
def exponential_model(x, a, b):
    """Exponential model: y = a * exp(b * x)"""
    return a * np.exp(b * x)

# Fit using curve_fit
popt, pcov = optimize.curve_fit(exponential_model, x, y, p0=[1.0, 0.5])
a_fit, b_fit = popt
perr = np.sqrt(np.diag(pcov))

print(f"Fitted parameters: a = {a_fit:.4f} ± {perr[0]:.4f}, b = {b_fit:.4f} ± {perr[1]:.4f}")

# Generate curve with fitted parameters
x_smooth = np.linspace(0, 5, 100)
y_fit = exponential_model(x_smooth, a_fit, b_fit)

# Plot results
plt.figure(figsize=(10, 6))
plt.scatter(x, y, color='blue', label='Data points')
plt.plot(x_smooth, y_fit, color='red', label=f'Fitted: y = {a_fit:.2f} * exp({b_fit:.2f} * x)')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.show()
```
{% endraw %}

### Resources

- **SciPy Documentation**:
  - [scipy.optimize.curve_fit](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html)
  - [scipy.optimize.minimize](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html)

## 6. General Linear Least Squares (GLLS)

### Key Concepts

- **Definition**: A method for fitting models that are linear in their parameters but may include nonlinear basis functions
- **General form**: $y = \beta_1 Z_1(x) + \beta_2 Z_2(x) + ... + \beta_m Z_m(x)$
- **Key components**:
  - Basis functions $Z_j(x)$ (can be nonlinear)
  - Coefficients $\beta_j$ (linear parameters to be determined)
  - Design matrix construction
  - Linear least squares solution
- **Applications**:
  - Multiple linear regression
  - Polynomial regression
  - Fourier series fitting
  - Custom basis function models

### GLLS Implementation Steps

1. **Define the basis functions** $Z_j(x)$ appropriate for your problem
2. **Construct the design matrix** Z where each column corresponds to a basis function evaluated at all data points
3. **Solve the linear system** using `np.linalg.lstsq(Z, y)`
4. **Calculate fit statistics** to evaluate the quality of the fit

### Example: Thermistor Calibration with Steinhart-Hart Equation

{% raw %}
```python
import numpy as np
import matplotlib.pyplot as plt

# Thermistor data: resistance (ohms) and temperature (Celsius)
R_data = np.array([6519, 4212, 2815, 1921, 1340, 944, 682, 496])
T_data = np.array([0, 10, 20, 30, 40, 50, 60, 70])

# Steinhart-Hart equation: 1/T = A + B*ln(R) + C*[ln(R)]^3
# This is already in linear form with respect to parameters A, B, and C

# Create design matrix Z with columns for each basis function
lnR = np.log(R_data)
Z = np.column_stack((np.ones_like(R_data), lnR, lnR**3))

# Create the dependent variable vector (1/T)
y = 1 / (T_data + 273.15)  # Convert to Kelvin for the model

# Solve using np.linalg.lstsq
coefficients, residuals, rank, s = np.linalg.lstsq(Z, y, rcond=None)
A, B, C = coefficients

print(f"Coefficients (A, B, C): {A:.8f}, {B:.8f}, {C:.8f}")

# Create prediction function
def predict_temperature(R, coefficients):
    """
    Predict temperature from resistance using the Steinhart-Hart equation.
    Returns temperature in Celsius.
    """
    A, B, C = coefficients
    lnR = np.log(R)
    inv_T = A + B * lnR + C * lnR**3  # 1/T in Kelvin
    T_kelvin = 1 / inv_T
    return T_kelvin - 273.15  # Convert back to Celsius

# Generate predicted temperatures for original data points
T_predicted = predict_temperature(R_data, coefficients)
```
{% endraw %}

### Resources

- **Davishahl Numerical Methods Videos**:
  - [General Linear Least Squares](https://youtu.be/593YmHgIoMI?si=Ey3AWmC-gijkjbqc)

## 7. Fit Quality Statistics

### Key Concepts

- **Purpose**: Quantify how well a model fits the data
- **Common Statistics**:
  - Coefficient of Determination (R²)
  - Standard Error of the Estimate (Syx)
  - Mean Squared Error (MSE)
  - Root Mean Squared Error (RMSE)
  - Mean Absolute Error (MAE)
- **Importance**:
  - Model selection and comparison
  - Assessment of prediction accuracy
  - Detection of overfitting or underfitting
  - Reporting confidence in results

### Key Statistics Definitions

1. **Coefficient of Determination (R²)**:
   - Represents the proportion of variance in the dependent variable explained by the model
   - Formula: $R^2 = 1 - \frac{SS_{res}}{SS_{tot}}$
   - Where:
     - $SS_{res} = \sum_i (y_i - f_i)^2$ is the sum of squares of residuals
     - $SS_{tot} = \sum_i (y_i - \bar{y})^2$ is the total sum of squares
   - Ranges from 0 to 1, with 1 being a perfect fit
   - May be negative for nonlinear models or models without an intercept

2. **Standard Error of the Estimate (Syx)**:
   - Measures the average distance data points fall from the regression line
   - Formula: $S_{yx} = \sqrt{\frac{SS_{res}}{n-p}}$
   - Where:
     - $n$ is the number of observations
     - $p$ is the number of parameters in the model
   - Units match the dependent variable for easy interpretation

3. **Root Mean Squared Error (RMSE)**:
   - Measures the square root of the average squared differences between predicted and actual values
   - Formula: $RMSE = \sqrt{\frac{1}{n}\sum_i (y_i - f_i)^2}$
   - Penalizes large errors more than small ones

### Implementation in Python

{% raw %}
```python
def fit_quality_stats(y_observed, y_predicted, num_params):
    """
    Calculate fit quality statistics.
    
    Parameters:
    -----------
    y_observed : array-like
        Observed data values
    y_predicted : array-like
        Model predicted values
    num_params : int
        Number of parameters in the model
        
    Returns:
    --------
    dict
        Dictionary containing various fit statistics
    """
    n = len(y_observed)
    
    # Mean of observed data
    y_mean = np.mean(y_observed)
    
    # Calculate residuals
    residuals = y_observed - y_predicted
    
    # Sum of squares
    ss_residuals = np.sum(residuals**2)
    ss_total = np.sum((y_observed - y_mean)**2)
    
    # R-squared
    r_squared = 1 - (ss_residuals / ss_total)
    
    # Adjusted R-squared (accounts for model complexity)
    adjusted_r_squared = 1 - (ss_residuals / ss_total) * ((n - 1) / (n - num_params - 1))
    
    # Standard error of the estimate
    if n > num_params:  # Ensure we don't divide by zero
        std_error = np.sqrt(ss_residuals / (n - num_params))
    else:
        std_error = np.nan
    
    # Root mean squared error
    rmse = np.sqrt(np.mean(residuals**2))
    
    return {
        'r_squared': r_squared,
        'adjusted_r_squared': adjusted_r_squared,
        'std_error': std_error,
        'root_mean_squared_error': rmse
    }
```
{% endraw %}

## 8. Visualization Conventions for Curve Fitting

### Key Concepts

- **Purpose**: Effective visualization is essential for evaluating fit quality and communicating results
- **Common Elements**:
  - Raw data representation
  - Fitted curve visualization
  - Error bars when applicable
  - Residuals analysis

### Best Practices

1. **Data Point Representation**:
   - Use distinctive markers for data points
   - Include error bars if measurement uncertainty is known

2. **Fitted Curve Display**:
   - Use a smooth line with sufficient points
   - Use different line styles for multiple models
   - Extend the curve slightly beyond the data range

3. **Axis Scaling**:
   - Use logarithmic scales for data spanning multiple orders of magnitude
   - Include zero in the axis range when appropriate

4. **Annotation**:
   - Include fit equation on the plot with parameter values
   - Display fit quality statistics (R², RMSE)

### Essential Visualization Types

1. **Data and Curve Plot**: Shows both data points and fitted model
2. **Residuals Plot**: Helps identify patterns and systematic biases
3. **Multiple Model Comparison**: Compares different fitting approaches

## 9. Tips for Week 5 Assignments

1. **Understand the Physical Meaning** of model parameters before fitting. This helps with setting reasonable initial guesses for nonlinear regression.

2. **Choose the Right Fitting Method** based on the model structure:
   - Use linear regression for linear relationships
   - Use linearization when appropriate for nonlinear models
   - Use direct nonlinear regression for complex models

3. **Compare Multiple Approaches** when fitting nonlinear models:
   - Linearization followed by linear regression
   - Direct nonlinear regression with curve_fit
   - Direct nonlinear regression with minimize

4. **Always Visualize Results** including:
   - Data points with the fitted curve
   - Residuals to check for patterns
   - Multiple models for comparison

5. **Calculate and Report Fit Quality Statistics**:
   - R² to measure goodness of fit
   - Standard error to quantify prediction accuracy
   - Parameter confidence intervals when available

6. **Check for Overfitting** particularly with high-order polynomial models:
   - Compare models with different numbers of parameters
   - Use adjusted R² which penalizes model complexity

7. **Be Mindful of Transformations** and their effects on error distribution:
   - Logarithmic transformations compress large values and expand small values
   - Error minimization in transformed space may not match error minimization in original space

8. **Document Your Process** including:
   - Mathematical model and its physical significance
   - Transformation approach if used
   - Fitting method and implementation details
   - Quality assessment of the fit

These tips will help you successfully complete the Week 5 programming assignments and develop a solid understanding of curve fitting techniques for engineering applications.
