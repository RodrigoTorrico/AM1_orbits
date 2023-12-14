# Inicializacion de variables
from numpy import linspace, array, log10
from Orbits import Kepler
from Cauchy_Problem import Cauchy_Problem
from Temporal_Schemes import Euler, Inverse_Euler, RK4, Crank_Nicolson
from Temporal_Error import Error_Richardson, convergence_rate
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
    
   
    plt.axis('equal')  
    plt.plot( time_domain,E[:,0])
    titulo ='Error con metodo de '
    plt.title (titulo + method)
    plt.xlabel('timepo')
    plt.ylabel('Error')
    plt.show()
    
    print ("Order Euler ")
    order, log_E, log_E2, log_n = convergence_rate(time_domain,temporal_scheme, F, U0, n)
    
    print("Order = ", order)
    
    fig1 = plt.figure(figsize=(10, 4))
    
    ax1 = fig1.add_subplot(121)
    ax1.axis('equal')  
    ax1.plot(log_n,log_E)
    titulo ='Convergence rate '
    ax1.set_title (titulo + method)
    ax1.set_label('log(n)')
    ax1.set_ylabel('log(|E|)')
    
    ax2 = fig1.add_subplot(122)
    ax2.axis('equal')  
    ax2.plot(log_n,log_E2)
    titulo ='Convergence rate '
    ax2.set_title (titulo + method)
    ax2.set_xlabel('log(n)')
    ax2.set_ylabel('log(|E|)')
    
    plt.show()
    
Milestone3([1,0,0,1],10,2000)
