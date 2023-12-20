##############################################################################################################
##                                             MILESTONE 5                                                  ##
##                                           N body problem                                                 ##
##############################################################################################################
##                                                                                                          ##
##                                       Rodrigo Torrico Gijon                                              ##
##############################################################################################################

from numpy import  zeros, reshape, linspace
from ODEs.Temporal_Schemes import RK4
from ODEs.Cauchy_Problem import Cauchy_Problem
from Physics.NBody import N_Body
import matplotlib.pyplot as plt



def Initial_Conditions(Nc, Nb): 
 
    U0 = zeros(2*Nc*Nb)
    U1 = reshape(U0, (Nb, Nc, 2))  
    r0 = reshape(U1[:, :, 0], (Nb, Nc))     
    v0 = reshape(U1[:, :, 1], (Nb, Nc))

    # body 1 
    r0[0,:] = [ 1, 0, 0]
    v0[0,:] = [ 0, 0.4, 0]

    # body 2 
    v0[1,:] = [ 0, -0.4, 0] 
    r0[1,:] = [ -1, 0, 0]

    # body 3 
    r0[2, :] = [ 0, 1, 0 ] 
    v0[2, :] = [ -0.4, 0., 0. ] 
         
    # body 4 
    r0[3, :] = [ 0, -1, 0 ] 
    v0[3, :] = [ 0.4, 0., 0. ]  

    return U0 

def Integration_N_Body():
    
    def F(U,t):
        return N_Body(U, t, Nb, Nc)
    
    N = 1000
    Nb = 4
    Nc = 3
    T = (N + 1)*2*Nc*Nb
     
    time_domain = linspace(0, 4*3.14, N + 1)
    U0 = Initial_Conditions (Nc, Nb)
    U = Cauchy_Problem(time_domain, RK4, F, U0)
     
    U_p  = reshape( U, (N+1, Nb, Nc, 2) ) 
    r   = reshape( U_p[:, :, :, 0], (N+1, Nb, Nc) ) 
   
    for i in range(Nb):
        plt.plot(r[:, i, 0], r[:, i, 1], label=f'Cuerpo {i + 1}')
    plt.axis('equal')
    plt.title(f'Problema de los N cuerpos (N = {Nb}) - 2D')
    plt.xlabel('Coordenada X')
    plt.ylabel('Coordenada Y')
    plt.legend()
    plt.grid()
    plt.show()
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i in range(Nb):
        ax.plot(r[:, i, 0], r[:, i, 1], time_domain, label=f'Cuerpo {i + 1}')

    ax.set_xlabel('Coordenada X')
    ax.set_ylabel('Coordenada Y')
    ax.set_zlabel('Coordenada Z')
    ax.set_title(f'Problema de los N cuerpos (N = {Nb}) - 3D')
    ax.legend()
    plt.show()
    
Integration_N_Body()