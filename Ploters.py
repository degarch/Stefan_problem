import numpy as np
import matplotlib.pyplot as plt
from Methods import interface_xi, solution, bisection
from Auxiliar_functions import f

def plot_phase_front(alpha, tau_range, delta):
    """
    Plot the phase front xi(tau) for fixed alpha.
    
    :param alpha: fractional order
    :param tau_range: list of tau values
    :param delta: scaling factor
    """
    xi_values = [interface_xi(tau, alpha, delta) for tau in tau_range]
    
    plt.figure(figsize=(10, 6))
    plt.plot(xi_values, tau_range, label=f'alpha = {alpha}')
    plt.title('Phase Front xi(tau)')
    plt.xlabel('xi(tau)')
    plt.ylabel('tau')
    plt.legend()
    plt.grid()
    plt.show()


def plot_solution_3D(alpha, y_min, y_max, tau_min, tau_max, M, delta, Bi, U_0, U_inf, U_m, y_points=100, tau_points=100):
    """
    Plot 3D surface of the temperature solution.
    
    :param alpha: fractional order
    :param y_min, y_max: spatial domain limits
    :param tau_min, tau_max: temporal domain limits
    :param M: number of series terms
    :param delta: scaling parameter for phase front
    :param Um: exterior temperature
    :param y_points, tau_points: resolution of the mesh
    """
    y = np.linspace(y_min, y_max, y_points)
    tau = np.linspace(tau_min, tau_max, tau_points)
    Y, T = np.meshgrid(y, tau)

    Z = np.array([[solution(alpha, yi, taui, M, delta, Bi, U_0, U_inf, U_m) for yi, taui in zip(y_row, t_row)]
                  for y_row, t_row in zip(Y, T)])

    # Mask values beyond the moving boundary
    for i in range(len(y)):
        for j in range(len(tau)):
            if y[i] > interface_xi(tau[j], alpha, delta):
                Z[j, i] = U_m

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(Y, T, Z, cmap='viridis', edgecolor='none')
    ax.set_title('3D Temperature Solution')
    ax.set_xlabel('y')
    ax.set_ylabel('tau')
    ax.set_zlabel('Temperature')
    fig.colorbar(surf, shrink=0.5, aspect=10)
    plt.show()


def plot_multiple_phase_fronts(alphas, Stefans, Biots, U_inf, U_0, U_m, tau_range, M, N, L):
    """
    Plot multiple phase fronts for different alpha values.
    
    :param alphas: list of alpha values
    :param Stefans: list of Stefan number values
    :param Biots: list of Biot number values
    :param tau_range: list of tau values
    :param M: number of series terms
    :param N: number of bisection iterations
    :param L: spatial domain length
    """
    plt.figure(figsize=(10, 6))

    for Ste in Stefans:
        for Bi in Biots:
            for alpha in alphas:
                delta = bisection(f, alpha, Ste, Bi, U_inf, U_0, U_m, 0, L, M, N)
                if delta==None:
                    break
                print(f"Alpha = {alpha}, Delta = {delta}")
                xi_values = [interface_xi(tau, alpha, delta) for tau in tau_range]
                if len(alphas)!=1:
                    plt.plot(xi_values, tau_range, label=f'α = {alpha}')
                if len(Stefans)!=1:
                    plt.plot(xi_values, tau_range, label=f'Ste = {Ste}')
                if len(Biots)!=1:
                    plt.plot(xi_values, tau_range, label=f'Bi = {Bi}')
                    
    if len(alphas)!=1:
        plt.title('Phase Fronts for Multiple α')
    if len(Stefans)!=1:
        plt.title('Phase Fronts for Multiple Ste numbers')
    if len(Biots)!=1:
        plt.title('Phase Fronts for Multiple Bi numbers')
    plt.xlabel('xi(tau)')
    plt.ylabel('tau')
    plt.legend()
    plt.grid()
    plt.show()
