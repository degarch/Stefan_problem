# MAIN.py
import numpy as np
from Methods import bisection
from Auxiliar_functions import f
from Ploters import plot_phase_front, plot_solution_3D, plot_multiple_phase_fronts


Ste = 1.0              # Stefan number
Bi = 1.0               # Biot number
U_0 = 1.0              # Characteristic temperature
U_m = 0.0              # Melting temperature 
U_inf = 1.0            # Ambient temperature
alpha = 0.5            # Fractional order
L = 2                  # Length of the domain
# Parameters
N = 50                 # Number of iterations for bisection method
M = 30                 # Number of terms in the series expansion

# Time parameter discretization 
tau_min, tau_max = 0.001, 5
tau_range = np.linspace(tau_min, tau_max, 100)

# Compute delta using bisection method
delta = bisection(f, alpha, Ste, Bi, U_inf, U_0, U_m, 0, L, M, N)
print("Computed delta:", delta)

# Plot phase front for fixed alpha
#plot_phase_front(alpha, tau_range, delta)

# Plot solution surface
#plot_solution_3D(alpha, 0, L, tau_min, tau_max, M, delta, Bi, U_0, U_inf, U_m)


# ----------------------------------Experiment 1-----------------------------------
# Plot the evolution of the phase front for different values of the fractional order α

plot_multiple_phase_fronts([0.5, 0.7, 0.8, 1.0], [Ste], [Bi], U_inf, U_0, U_m, tau_range, M, N, L)

# ----------------------------------Experiment 2-----------------------------------
# Plot the evolution of the phase front for different values of the Stefan number

plot_multiple_phase_fronts([0.5], [0.001,0.1,0.5,1,2,3,3.5,4,4.5,5], [Bi], U_inf, U_0, U_m, tau_range, M, N, L)

# ----------------------------------Experiment 1-----------------------------------
# Plot the evolution of the phase front for different values of the Biot number

plot_multiple_phase_fronts([0.5], [Ste], [0.001, 0.1, 1, 10], U_inf, U_0, U_m, tau_range, M, N, L)