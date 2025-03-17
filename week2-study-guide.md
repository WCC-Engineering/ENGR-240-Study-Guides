# Week 2 Study Guide: Functions, Loops, and Conditionals

## Overview

This study guide covers the fundamental Python concepts related to program flow control, which are crucial for implementing computational algorithms in engineering. The activities, worksheets, and programming assignments for Week 2 are designed to provide practice with these concepts and techniques.

## 1. Functions in Python

### Key Concepts

- **Function Definition**: A reusable block of code that performs a specific task
- **Function Components**:
  - `def` keyword to define a function
  - Function name
  - Parameters (inputs)
  - Function body
  - Return statement (optional)
- **Function Documentation**: Using docstrings to document function behavior
- **Function Scope**: Variables defined inside a function are local to that function
- **Return Values**: Functions can return results to the caller

### Example Code

```python
# A simple function to calculate the area of a circle
def calculate_circle_area(radius):
    """
    Calculate the area of a circle given its radius.
    
    Parameters:
    -----------
    radius : float
        The radius of the circle
        
    Returns:
    --------
    area : float
        The area of the circle
    """
    import math
    area = math.pi * radius**2
    return area

# Call the function
radius = 5
area = calculate_circle_area(radius)
print(f"The area of a circle with radius {radius} is {area:.2f}")
```

### Lambda Functions

Lambda functions are small, anonymous functions defined with the `lambda` keyword.

```python
# Traditional function
def square(x):
    return x**2

# Equivalent lambda function
square_lambda = lambda x: x**2

# Using the lambda function
print(square_lambda(5))  # Output: 25

# Lambda functions are useful for short operations
# Example: Using a lambda function as an argument to another function
numbers = [1, 5, 3, 9, 2, 6]
sorted_numbers = sorted(numbers, key=lambda x: x**2)
print(sorted_numbers)  # Sort based on square of each number
```

### Resources

- **Python Documentation**:
  - [Defining Functions](https://docs.python.org/3/tutorial/controlflow.html#defining-functions)
  - [Lambda Expressions](https://docs.python.org/3/tutorial/controlflow.html#lambda-expressions)
 
- **Berkeley Python Numerical Methods**:
  - [Functions](https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter03.00-Functions.html)
 
- **Sullivan Numerical Methods YouTube Playlist**:
  - [Functions and Lambda Functions](https://youtu.be/ItNJaUiieog?si=UBOAjk4AsrtLEc6U)
  - [More Complicated Python Functions](https://youtu.be/SnceGI35xaw?si=3htda-le280sWAvC)
 
- **Python for Everybody Videos**
  - [Functions part 1](https://youtu.be/5Kzw-0-DQAk?si=c___UF6gJPB4yFQ2)
  - [Functions part 2](https://youtu.be/AJVNYRqn8kM?si=5wHG5ABxYbtDznt5)

## 2. Loops for Iteration

### Key Concepts

- **For Loops**: Iterate over a sequence (lists, tuples, strings, etc.)
- **Range Function**: Generate a sequence of numbers for iteration
  - `range(stop)`: 0 to stop-1
  - `range(start, stop)`: start to stop-1
  - `range(start, stop, step)`: start to stop-1 with step increment
- **Loop Control Statements**:
  - `break`: Exit the loop completely
  - `continue`: Skip the current iteration and continue with the next
- **Nested Loops**: Loops inside loops for multidimensional iteration
- **Loop Variables**: Keep track of the current element in each iteration

### Example Code

```python
# Example 1: Basic for loop with range
for i in range(5):
    print(i)  # Prints 0, 1, 2, 3, 4

# Example 2: Loop through a list
fruits = ['apple', 'banana', 'cherry']
for fruit in fruits:
    print(f"I like {fruit}s")

# Example 3: Using break and continue
for i in range(10):
    if i == 3:
        continue  # Skip 3
    if i == 7:
        break  # Stop at 7
    print(i)  # Prints 0, 1, 2, 4, 5, 6

# Example 4: Nested loops
for i in range(3):
    for j in range(2):
        print(f"({i}, {j})")
```

### Series Calculations with Loops

Loops are particularly useful for calculating series, which are common in engineering analysis:

```python
# Calculate the first N terms of a geometric series with ratio 1/2
def geometric_series(N):
    """
    Calculate terms and sum of geometric series 1/2 + 1/4 + 1/8 + ...
    
    Parameters:
    -----------
    N : int
        Number of terms to calculate
        
    Returns:
    --------
    terms : list
        Individual terms of the series
    series_sum : float
        Sum of the first N terms
    """
    terms = []
    series_sum = 0
    
    for k in range(1, N+1):
        term = 1 / (2**k)
        terms.append(term)
        series_sum += term
    
    return terms, series_sum

# Calculate first 10 terms
terms, total = geometric_series(10)
print(f"Terms: {terms}")
print(f"Sum: {total}")
```

### Resources

- **Python Documentation**:
  - [For Statements](https://docs.python.org/3/tutorial/controlflow.html#for-statements)
  - [The range() Function](https://docs.python.org/3/tutorial/controlflow.html#the-range-function)
 
- **Berkeley Python Numerical Methods**:
  - [Iteration](https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter05.00-Iteration.html)

- **Sullivan Numerical Methods YouTube Playlist**:
  - [For Loops](https://youtu.be/plOu_VAFUSg?si=WmvDvCpJ7VlPuv6o)
  - [Unpacking FOR Loops](https://youtu.be/askBlkCnB5s?si=xIHe_ykmLDXN4x2-)
 
- **Python for Everybody Videos**
  - [Loops and Iteration part 1](https://youtu.be/FzpurxjwmsM?si=bq_2iJw7wa4EW7Jq)
  - [Loops and Iteration part 2](https://youtu.be/5QDrj5ogPYc?si=D4VoBhraGi603fwZ)
  - [Loops and Iteration part 3](https://youtu.be/xsavQp8hd78?si=Znk48wDCNn_RVHhr)
  - [Loops and Iteration part 4](https://youtu.be/yjlMMwf9Y5I?si=hvS1W2_M6sd5SmPH)

## 3. Conditional Statements

### Key Concepts

- **Simple Conditionals**: `if`, `else` statements
- **Multiple Conditions**: `if`, `elif`, `else` statements
- **Comparison Operators**: `==`, `!=`, `<`, `>`, `<=`, `>=`
- **Logical Operators**:
  - `and`: Both conditions must be true
  - `or`: At least one condition must be true
  - `not`: Negates a condition
- **Nested Conditionals**: Conditionals inside other conditionals
- **Conditional Expressions**: Compact one-line conditionals (ternary operators)

### Example Code

```python
# Example 1: Simple if-else
temperature = 28
if temperature > 25:
    print("It's warm outside")
else:
    print("It's not that warm")

# Example 2: Multiple conditions
score = 85
if score >= 90:
    grade = 'A'
elif score >= 80:
    grade = 'B'
elif score >= 70:
    grade = 'C'
elif score >= 60:
    grade = 'D'
else:
    grade = 'F'
print(f"Grade: {grade}")

# Example 3: Logical operators
temp = 22
humidity = 80
if temp > 25 and humidity > 60:
    print("Hot and humid")
elif temp > 25 or humidity > 60:
    print("Either hot or humid")
else:
    print("Neither hot nor humid")

# Example 4: Conditional expression (ternary operator)
age = 20
status = "Adult" if age >= 18 else "Minor"
print(status)
```

### Flow Regime Classification Example

This example demonstrates how to use conditionals to classify flow regimes in fluid mechanics:

```python
def calculate_reynolds(velocity, diameter, kinematic_viscosity):
    """Calculate the Reynolds number for pipe flow."""
    reynolds = velocity * diameter / kinematic_viscosity
    return reynolds

def classify_flow_regime(reynolds):
    """
    Classify the flow regime based on the Reynolds number.
    
    Parameters:
    ----------
    reynolds : float
        Reynolds number
        
    Returns:
    -------
    regime : str
        Flow regime ('laminar', 'transitional', or 'turbulent')
    """
    if reynolds < 2300:
        regime = 'laminar'
    elif reynolds <= 4000:
        regime = 'transitional'
    else:
        regime = 'turbulent'
    
    return regime

# Example usage
water_viscosity = 1.004e-6  # m²/s
pipe_diameter = 0.05  # m
flow_velocity = 1.0  # m/s

re = calculate_reynolds(flow_velocity, pipe_diameter, water_viscosity)
flow_type = classify_flow_regime(re)
print(f"Reynolds number: {re:.0f}")
print(f"Flow regime: {flow_type}")
```

### Resources

- **Python Documentation**:
  - [If Statements](https://docs.python.org/3/tutorial/controlflow.html#if-statements)
 
- **Berkeley Python Numerical Methods**:
  - [Conditional Statements](https://pythonnumericalmethods.berkeley.edu/notebooks/chapter03.02-Conditional-Statements.html)

## 4. Combining Functions, Loops, and Conditionals

### Key Concepts

- **Modular Programming**: Breaking complex problems into smaller, manageable parts
- **Function Composition**: Using the output of one function as input to another
- **Control Flow Patterns**:
  - Looping with conditional exits
  - Processing collections with filters
  - Implementing mathematical algorithms

### Example Code: Friction Factor Calculation

This example shows how to integrate functions, loops, and conditionals to solve an engineering problem:

```python
def calculate_reynolds(velocity, diameter, kinematic_viscosity):
    """Calculate the Reynolds number for pipe flow."""
    reynolds = velocity * diameter / kinematic_viscosity
    return reynolds

def calculate_friction_factor(reynolds, relative_roughness=0):
    """
    Calculate the Darcy friction factor based on Reynolds number and relative roughness.
    
    Parameters:
    ----------
    reynolds : float
        Reynolds number
    relative_roughness : float, optional
        Relative pipe roughness (ε/D), dimensionless
        Default is 0 (smooth pipe)
        
    Returns:
    -------
    friction_factor : float
        Darcy friction factor
    """
    # Define the friction factor equations as lambda functions
    laminar_f = lambda re: 64 / re
    
    # Blasius equation for turbulent flow in smooth pipes
    blasius_f = lambda re: 0.316 / (re ** 0.25)
    
    # Haaland equation for turbulent flow in rough pipes
    haaland_f = lambda re, roughness: (-1.8 * np.log10((6.9/re) + (roughness/3.7)**1.11))**(-2)
    
    # Implement conditional logic based on flow regime
    if reynolds < 2300:  # Laminar flow
        friction_factor = laminar_f(reynolds)
    elif reynolds <= 4000:  # Transitional flow
        # Linear interpolation between laminar and turbulent values
        f_laminar = laminar_f(2300)
        f_turbulent = blasius_f(4000) if relative_roughness == 0 else haaland_f(4000, relative_roughness)
        
        # Linear interpolation factor
        t = (reynolds - 2300) / (4000 - 2300)
        friction_factor = f_laminar + t * (f_turbulent - f_laminar)
    else:  # Turbulent flow
        if relative_roughness == 0:
            # Smooth pipe - use Blasius equation
            friction_factor = blasius_f(reynolds)
        else:
            # Rough pipe - use Haaland equation
            friction_factor = haaland_f(reynolds, relative_roughness)
    
    return friction_factor

# Example of using these functions together
import numpy as np
import matplotlib.pyplot as plt

# Define parameters for water in a steel pipe
water_viscosity = 1.004e-6  # m²/s
pipe_diameter = 0.05  # m
relative_roughness = 0.0001  # A typical value for commercial steel pipes

# Create a range of velocities to analyze
velocities = np.linspace(0.01, 3.0, 100)

# Calculate Reynolds numbers and friction factors
reynolds_numbers = []
friction_factors = []

for velocity in velocities:
    re = calculate_reynolds(velocity, pipe_diameter, water_viscosity)
    reynolds_numbers.append(re)
    f = calculate_friction_factor(re, relative_roughness)
    friction_factors.append(f)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(reynolds_numbers, friction_factors)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Reynolds Number')
plt.ylabel('Friction Factor')
plt.title('Friction Factor vs. Reynolds Number')
plt.grid(True, which='both', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()
```

## 5. Best Practices

### Function Design

- **Single Responsibility**: Functions should do one thing well
- **Meaningful Names**: Function names should indicate what they do, not how they do it
- **Documentation**: Use docstrings to document parameters, return values, and behavior
- **Error Handling**: Consider how your function will handle invalid inputs
- **Default Parameters**: Use default parameter values for optional arguments

### Control Flow Design

- **Early Returns**: Exit functions early for special cases
- **Avoid Deep Nesting**: Too many nested conditionals or loops make code hard to read
- **Keep It Simple**: Complicated control flow is harder to debug and maintain
- **Be Consistent**: Use the same style for similar operations

### Code Example: Early Returns vs Deep Nesting

```python
# Deep nesting approach - harder to read
def process_data(data):
    if data is not None:
        if len(data) > 0:
            if isinstance(data[0], (int, float)):
                # Process the data
                return sum(data)
            else:
                return None
        else:
            return None
    else:
        return None

# Early returns approach - cleaner
def process_data_better(data):
    if data is None:
        return None
    if len(data) == 0:
        return None
    if not isinstance(data[0], (int, float)):
        return None
        
    # Process the data
    return sum(data)
```

## Tips for Week 2 Assignments

1. **Start with function definitions**: Clearly define the inputs and outputs of each function before implementing the details.

2. **Test incrementally**: Test each function separately before integrating them together.

3. **Debug using print statements**: Insert print statements to check values at different stages of your program.

4. **Use helper functions**: Break complex operations into smaller, reusable helper functions.

5. **Validate inputs**: Check that inputs to your functions are valid and handle errors gracefully.

6. **Comment your code**: Add comments to explain your logic, especially for complex conditionals or loop structures.

7. **Use docstrings**: Document your functions with docstrings to clarify their purpose, parameters, and return values.
