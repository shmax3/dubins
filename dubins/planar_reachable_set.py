from numpy import linspace, array, pi, sin, cos, tan, arccos, sqrt, nan_to_num
from numpy import arcsin, where
from scipy.optimize import brenth


N = 100


def D_boundary(N=N):
    tau = linspace(0, 2*pi, N)
    return 1 + cos(tau), sin(tau), -1 + cos(tau), sin(tau)

def D(x, y):
    return (1 - abs(x))**2 + y**2 <= 1

def A_boundary(N=N):
    tau = linspace(0, arcsin(2*sqrt(2)/3), N)
    alpha = linspace(0, pi, N)
    return (-1 + 3*cos(tau), 3*sin(tau), 1 - cos(alpha), sin(alpha),
            1 - 3*cos(tau), 3*sin(tau), -1 - cos(alpha), sin(alpha))

def A(x, y):
    return (y >= 0) & ((1+abs(x))**2 + y**2 <= 9) & ((1-abs(x))**2 + y**2 >= 1)

def x_RS(tau, t):
    return (t-tau)*sin(tau) - (cos(tau) - 1)

def y_RS(tau, t):
    return (t-tau)*cos(tau) + sin(tau)

def x_LS(tau, t):
    return -x_RS(tau, t) 

def y_LS(tau, t):
    return y_RS(tau, t)

def x_LR(tau, t):
    return (2-cos(t-tau))*cos(tau) - sin(t-tau)*sin(tau) - 1

def y_LR(tau, t):
    return (2-cos(t-tau))*sin(tau) + sin(t-tau)*cos(tau)

def x_RL(tau, t):
    return -x_LR(tau, t) 

def y_RL(tau, t):
    return y_LR(tau, t)

def B(t, N=N, M=N):
    return (*B_CS(t, N), *B_CC(t, M))

def ranges_for_B_CS(t, N):
    if t <= (3*pi/2 + 1):
        tau1 = linspace(0, t, N)
        tau2 = array([])
    elif t <= 2*pi:
        tau1_max = brenth(lambda tau: x_RS(tau, t), pi, 3*pi/2)
        tau2_min = brenth(lambda tau: x_RS(tau, t), 3*pi/2, t)
        tau1 = linspace(0, tau1_max, int(N*tau1_max/(t-tau2_min+tau1_max)))
        tau2 = linspace(tau2_min, t, int(N*(t-tau2_min)/(t-tau2_min+tau1_max)))
    else:
        tau1_max = brenth(lambda tau: x_RS(tau, t), pi, 3*pi/2)
        tau1 = linspace(0, tau1_max, N)
        tau2 = array([])
    return tau1, tau2

def B_RS(t, N=N):
    tau1, tau2 = ranges_for_B_CS(t, N)
    return x_RS(tau1, t), y_RS(tau1, t), x_RS(tau2, t), y_RS(tau2, t)

def B_LS(t, N=N):
    tau1, tau2 = ranges_for_B_CS(t, N)
    return x_LS(tau1, t), y_LS(tau1, t), x_LS(tau2, t), y_LS(tau2, t)

def B_CS(t, N=N):
    return (*B_RS(t, int(N/2)), *B_LS(t, int(N/2)))

def range_for_B_CC(t, N):
    """For the tau_min and tau_max analytical solution could be found,
       since the equation could be simplified to quartic."""
    if t >= 2*pi + arccos(23/27):
        tau = array([])
    else:
        if t <= 2*pi:
            tau_min = 0
        else:
            tau_min = brenth(lambda tau: x_LR(tau, t), t-2*pi, max(0, (t-pi)/3))
        tau_max = brenth(lambda tau: x_LR(tau, t),
                            max(0, (t-pi)/3), min(t, (t+pi)/3))
        tau = linspace(tau_min, tau_max, N)
    return tau

def B_LR(t, N=N):
    tau = range_for_B_CC(t, N)
    return x_LR(tau, t), y_LR(tau, t)

def B_RL(t, N=N):
    tau = range_for_B_CC(t, N)
    return x_RL(tau, t), y_RL(tau, t)

def B_CC(t, N=N):
    return (*B_RL(t, int(N/2)), *B_LR(t, int(N/2))) 

def alpha_C(x, y):
    return ((1 - abs(x) + y*sqrt((1-abs(x))**2 + y**2 - 1))
            / ((1-abs(x))**2 + y**2))

def alpha_S(x, y):
    return (where(y > 0, x != 0, 1)
            * (y - (1 - abs(x))*sqrt((1-abs(x))**2 + y**2 - 1))
            / ((1-abs(x))**2 + y**2))

def theta_CS(x, y):
    s = alpha_S(x, y) >= 0
    return (2*s-1)*arccos(alpha_C(x, y)) + 2*(1-s)*pi

def V_CS(x, y):
    return theta_CS(x, y) + sqrt((1-abs(x))**2 + y**2 - 1)

def alpha(x, y):
    return (5 - (1 + abs(x))**2 - y**2) / 4

def theta_CC_plus(x, y):
    return arccos(((1+abs(x))*(2-alpha(x, y)) + y*sqrt(1-alpha(x, y)**2))
                  / ((1+abs(x))**2 + y**2))

def V_CC_plus(x, y):
    A = ((y >= 0) & ((abs(x)- 1)**2 + y**2 >= 1)
         & ((1+abs(x))**2 + y**2 <= 9))
    return theta_CC_plus(x, y)*sqrt(2*A-1) + arccos(alpha(x, y))

def theta_CC_minus(x, y):
    return arccos(((1+abs(x))*(2-alpha(x, y)) - y*sqrt(1-alpha(x, y)**2))
                  / ((1+abs(x))**2 + y**2))

def V_CC_minus(x, y):
    return theta_CC_minus(x, y) + 2*pi - arccos(alpha(x, y))

def R(t, x, y):
    D = (abs(x)- 1)**2 + y**2 <= 1
    A = ((y >= 0) & ((abs(x)- 1)**2 + y**2 >= 1)
         & ((1+abs(x))**2 + y**2 <= 9))
    abnormal = (A & nan_to_num(((t >= V_CC_minus(x, y))
                                | (t <= V_CC_plus(x, y)))
                               & (t >= V_CS(x, y)), nan=False))
    normal = (~A & ((t >= V_CS(x, y))
                    | ((t >= V_CC_minus(x, y)) & D)))
    return normal | abnormal