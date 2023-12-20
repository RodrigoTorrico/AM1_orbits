from numpy import  zeros, reshape
from numpy.linalg import norm

def N_Body(U, t, Nb, Nc): 
     
     #Nb = len(mass)
     # Se crea un puntero donde Nb es el numero de cuerpos y Nc el numero de coordenadas    
     U_p  = reshape( U, (Nb, Nc, 2) )  
     F =  zeros(len(U))   
     dU_p = reshape( F, (Nb, Nc, 2) )  
     
     # Se crean punteros posicion y velocidad. Esto hace el codigo mas facil de seguir
     r = reshape( U_p[:, :, 0], (Nb, Nc) )     
     v = reshape( U_p[:, :, 1], (Nb, Nc) )
     
     # Se crean punteros para las derivadas
     drdt = reshape(dU_p[:, :, 0], (Nb, Nc) )
     dvdt = reshape(dU_p[:, :, 1], (Nb, Nc) )
    
     dvdt[:,:] = 0  

     # Se introduce el problema de los N cuerpos
     for i in range(Nb):   
        drdt[i,:] = v[i,:]
        for j in range(Nb): 
            if j != i:  
                d = r[j,:] - r[i,:]
                dvdt[i,:] = dvdt[i,:] + d[:] / norm(d)**3     
     return F