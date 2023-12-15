from numpy import linspace,zeros, transpose
from Temporal_Schemes import Euler, RK4, Inverse_Euler, Crank_Nicolson
import matplotlib.pyplot as plt


def Stability_Region(Temporal_Scheme, x0, y0, xf, yf, N):
    x, y = linspace(x0, xf, N), linspace (y0, yf, N)
    rho = zeros((N, N))
    for i in range(N): 
      for j in range(N):

          w = complex(x[i], y[j])
          r = Temporal_Scheme( 1., 1., 0., lambda u, t: w*u )
          rho[i, j] = abs(r) 

    return rho, x, y 

def Test():
   Temporal_Schemes = [Euler, Inverse_Euler, RK4, Crank_Nicolson]
   for Temporal_Scheme in Temporal_Schemes:
        rho, x, y = Stability_Region(Temporal_Scheme, -4, -4, 4, 4, 100)
        plt.contour( x, y, transpose(rho), linspace(0, 1, 11) )
        plt.axis('equal')
        plt.grid()
        plt.title ('Region de estabilidad de ' + Temporal_Scheme.__name__)
        plt.xlabel('Eje Re')
        plt.ylabel('Eje Im')
        plt.show()
        
if __name__ == '__main__':  
    Test()
