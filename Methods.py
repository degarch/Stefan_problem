from Auxiliar_functions import f, int_sigma
import math

def bisection(f, alpha, Ste, Bi, U_inf, U_0, U_m, a, b, M, N):
    """
    Bisection method to find the root of f(y) = 0.
    
    :param f: function for which the root is sought
    :param alpha: fractional order
    :param a, b: interval endpoints
    :param M: number of series terms
    :param N: number of iterations
    :return: approximate root
    """
    if f(alpha, a, M, Ste, Bi, U_inf, U_0, U_m) * f(alpha, b, M, Ste, Bi, U_inf, U_0, U_m) >= 0:
        print("Bisection method failed: invalid interval for delta.")
        return None
    a_n = a
    b_n = b
    for _ in range(N):
        m_n = (a_n + b_n) / 2
        f_m_n = f(alpha, m_n, M, Ste, Bi, U_inf, U_0, U_m)
        if f(alpha, a_n, M, Ste, Bi, U_inf, U_0, U_m) * f_m_n < 0:
            b_n = m_n
        elif f(alpha, b_n, M, Ste, Bi, U_inf, U_0, U_m) * f_m_n < 0:
            a_n = m_n
        else:
            return m_n
    return (a_n + b_n) / 2


def interface_xi(tau, alpha, delta):
    """
    Compute the moving boundary position interface_xi(tau).
    
    :param tau: time variable
    :param alpha: fractional order
    :param delta: scaling constant
    :return: interface_xi(tau)
    """
    return delta * math.pow(tau, 1 / (1 + alpha))


def solution(alpha, y, tau, M, delta, Bi, U_0, U_inf, U_m):
    """
    Compute the explicit temperature solution.
    
    :param alpha: fractional order
    :param y: spatial coordinate
    :param tau: time
    :param M: number of terms in the series
    :param delta: front scaling parameter
    :param U0: temperature at fixed boundary
    :param Um: temperature at moving boundary
    :param h0: heat transfer parameter
    :return: temperature at (y, tau)
    """
    B = (U_m-U_inf)/U_0 *1/(int_sigma(alpha, delta, M) +  math.gamma(alpha) * (1 / Bi))
    A = U_inf/U_0 + B*  math.gamma(alpha)/ Bi
    transformed_y = y / math.pow(tau, 1 / (1 + alpha))
    return A + B * int_sigma(alpha, transformed_y, M)
