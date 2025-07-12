import math

def int_sigma(alpha, y, M):
    """
    Compute the series integral Σ for the Stefan problem.
    
    :param alpha: fractional order
    :param y: spatial variable
    :param M: number of series terms
    :return: approximate value of integral
    """
    result = math.pow(y, alpha) / alpha
    cn = 1.0
    for i in range(1, M + 1):
        cn *= math.gamma(i * (alpha + 1)) / math.gamma(i * (alpha + 1) + alpha)
        term = cn * math.pow(-1, i) / math.pow(1 + alpha, i) \
               * (1 / ((1 + alpha) * (i + 1) - 1)) \
               * math.pow(y, (i + 1) * (alpha + 1) - 1)
        result += term
    return result


def int_p_sigma(alpha, y, M):
    """
    Compute the series integral ∫σp for the Stefan problem.
    
    :param alpha: fractional order
    :param y: spatial variable
    :param M: number of series terms
    :return: approximate value of integral
    """
    result = math.pow(y, alpha + 1) / (alpha + 1)
    cn = 1.0
    for i in range(1, M + 1):
        cn *= math.gamma(i * (alpha + 1)) / math.gamma(i * (alpha + 1) + alpha)
        term = cn * math.pow(-1, i) / math.pow(1 + alpha, i) \
               * (1 / ((1 + alpha) * (i + 1))) \
               * math.pow(y, (i + 1) * (alpha + 1))
        result += term
    return result


def f(alpha, y, M, Ste=1.0, Bi=1.0, U_inf=1.0, U_0=1.0, Um=0.0):
    """
    Auxiliary function for finding delta using bisection.
    
    :param alpha: fractional order
    :param y: guess for phase boundary parameter
    :param M: number of terms
    :param U_0: caracteristic temperature  
    :param U_inf: ambient temperature 
    :param Um: boundary temperature at moving front
    :param Bi: Biot number
    :param Ste: Stefan number
    :return: f(y) = H(alpha, y) - y
    """
    integral_p = int_p_sigma(alpha, y, M)
    integral = int_sigma(alpha, y, M)
    H =Ste*(U_inf - Um)/U_0 * ((1 + alpha) * math.gamma(alpha) - integral_p) \
        / (integral + math.gamma(alpha) * (1 / Bi))
    return H - y

