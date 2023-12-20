##############################################################################################################
##                                             MILESTONE 2                                                  ##
##                             Prototypes to integrate orbits with functions                                ##
##############################################################################################################
##                                                                                                          ##
##                                       Rodrigo Torrico Gijon                                              ##
##############################################################################################################

from re import A
from numpy import array, linspace
from ODEs.Temporal_Schemes import Euler, Crank_Nicolson, Inverse_Euler, RK4,ERK, GBS_Scheme
from ODEs.Cauchy_Problem import Cauchy_Problem
from Physics.Orbits import Kepler, Arenstorf
import matplotlib.pyplot as plt



temporal_schemes = {
    "Euler": Euler,
    "Crank_Nicolson": Crank_Nicolson,
    "Inverse_Euler": Inverse_Euler,
    "RK4": RK4,
    "ERK":ERK,
    "GBS":GBS_Scheme
}



def Simulation(U0, tf, N):
    method = input("Introduzca metodo: ")
    
    if method in temporal_schemes:
        temporal_scheme = temporal_schemes[method]
    else:
        raise ValueError("Metodo no valido")

    time_domain = linspace(0, tf, N)
   # U = Cauchy_Problem(time_domain, temporal_scheme, Kepler, U0)
    U = Cauchy_Problem(time_domain, temporal_scheme, Arenstorf, U0,0)
    plt.axis('equal')
    plt.plot(U[:, 0], U[:, 1])
    title = 'Metodo de '
    plt.title(title + method)
    plt.xlabel('Eje x')
    plt.ylabel('Eje y')
    plt.show()

    return U

Simulation([0.994, 0, 0, -2.0015851063798025664053786222], 17.06521656015796, 10000)

