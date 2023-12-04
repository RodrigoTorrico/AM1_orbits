
from numpy import array,linspace
from Cauchy_Problem import Cauchy_Problem
from Temporal_Schemes import Euler, Crank_Nicolson, Inverse_Euler, RK4
from Orbits import Kepler
import matplotlib.pyplot as plt


temporal_schemes = {
    "Euler": Euler, 
    "Crank_Nicolson": Crank_Nicolson, 
    "Inverse_Euler": Inverse_Euler, 
    "RK4": RK4
    }


def Simulation(U0, tf, N):
    
    method = input("Introduzca metodo: ")
    
    if method in temporal_schemes:
        temporal_scheme = temporal_schemes[method]
    else:
        raise ValueError("Metodo no valido")
    
    time_domain = linspace(0, tf, N)
    
    U = Cauchy_Problem(time_domain, temporal_scheme, Kepler, U0)
    
    plt.axis('equal')  
    plt.plot( U[:,0], U[:,1],".")
    titulo ='Metodo de '
    plt.title (titulo + method)
    plt.xlabel('Eje x')
    plt.ylabel('Eje y')
    plt.show()
   
Simulation(array([1,0,0,1]),10,10000)  
