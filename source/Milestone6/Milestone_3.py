##############################################################################################################
##                                             MILESTONE 3                                                  ##
##                                Error estimation of numerical solutions                                   ##
##############################################################################################################
##                                                                                                          ##
##                                       Rodrigo Torrico Gijon                                              ##
##############################################################################################################

# Inicializacion de variables
from numpy import linspace, array, log10
from Physics.Orbits import Kepler
from ODEs.Temporal_Schemes import Euler, Inverse_Euler, RK4, Crank_Nicolson
from ODEs.Temporal_Error import Error_Richardson, convergence_rate
import matplotlib.pyplot as plt


temporal_schemes = {
    "Euler": Euler, 
    "Crank_Nicolson": Crank_Nicolson, 
    "Inverse_Euler": Inverse_Euler, 
    "RK4": RK4
    }


def Milestone3(U0, tf, N):
    
    method = input("Introduzca metodo: ")
    F = Kepler
    q = int(input("Introduzca orden: "))
    n = 5
    
    if method in temporal_schemes:
        temporal_scheme = temporal_schemes[method]
    else:
        raise ValueError("Metodo no valido")
    
    time_domain = linspace(0, tf, N)

    E, S = Error_Richardson(time_domain,temporal_scheme, F, U0, q)
    
   
    print(f'Numerical integration error for the Cauchy problem, using {method}')
    cauchy_error, U = Error_Richardson(time_domain, temporal_scheme, Kepler, U0, q)
    
    plt.figure (1)
    plt.plot(time_domain, cauchy_error[:,0], color = 'red' )
    plt.xlabel('t')
    plt.ylabel(f'Error Kepler orbit with {method}')
    plt.title(f'{method} time scheme')  
    plt.grid()

    
    
    print("Convergence rate of the inverse Euler time scheme")
    order, log_E, log_E2, log_n = convergence_rate(time_domain, temporal_scheme, Kepler, U0, n)
    
    print( "order =", order)
    
    plt.figure (2)
    plt.plot(log_n, log_E, color = 'red' )
    plt.xlabel('log(N)')
    plt.ylabel('log(E)')
    plt.title(f'{method} time scheme convergence rate')  
    plt.grid()
    plt.show()
    


Milestone3([1,0,0,1],10,2000)