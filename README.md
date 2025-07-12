
# Fractional Stefan Problem: Semi-Analytical Solution in 1D

This repository provides a Python implementation for studying the one-dimensional fractional Stefan problem with a convective boundary condition. The solution is constructed using a semi-analytical approach based on series expansions and numerical methods to compute the moving phase boundary.

## Files

- `MAIN.py`: Main script. Runs the full simulation, computes the moving boundary location using the bisection method, and generates plots for the phase front and temperature solution.
- `Auxiliar_functions.py`: Contains auxiliary functions used in constructing the solution (series integrals tailored to the model).
- `Methods.py`: Implements the bisection method to determine the scaling parameter `δ` (delta), and defines functions to compute the interface position and temperature.
- `Ploters.py`: Visualization module. Includes functions to plot the phase front, 3D temperature evolution, and multiple phase fronts for various α values.

## Requirements

- Python 3.x
- Libraries:
  - `numpy`
  - `matplotlib`

Install dependencies via:

```bash
pip install numpy matplotlib
```

## Usage

Run `MAIN.py` to compute and visualize the results using the parameters defined in the script. It performs the following steps:

1. Computes `δ` for a given α using the bisection method.
2. Plots the phase front ξ(τ).
3. Displays a 3D surface of the temperature solution `u(y, τ)`.
4. Compares phase fronts for multiple fractional orders α.

## Physical Parameters

The model depends on several physical and numerical parameters:

- `α`: fractional order of the time derivative.
- `Ste`: Stefan number.
- `Bi`: Biot number.
- `U_0`: characteristic temperature.
- `U_m`: melting temperature.
- `U_inf`: ambient temperature. 
- `M`: number of terms in the series expansion.
- `N`: number of iterations in the bisection method.

These can be modified directly within `MAIN.py`.

## References

This project is based on the fractional formulation of the classical Stefan problem, commonly used in phase-change and heat transfer modeling. The approach considers a convective-type boundary condition and solves the problem via a semi-analytical series method.
