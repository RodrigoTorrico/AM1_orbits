##############################################################################################################
##                                             MILESTONE 4                                                  ##
##                           Linear problems. Regions of absolute stability                                 ##
##############################################################################################################
##                                                                                                          ##
##                                       Rodrigo Torrico Gijon                                              ##
##############################################################################################################

from numpy import linspace,zeros, transpose
from ODEs.Temporal_Schemes import Euler, RK4, Inverse_Euler, Crank_Nicolson
from ODEs.Stability_Region import Stability_Region
import matplotlib.pyplot as plt


 

def Milestone4():
   Temporal_Schemes = [Euler, Inverse_Euler, RK4, Crank_Nicolson]
   for Temporal_Scheme in Temporal_Schemes:
        rho, x, y = Stability_Region(Temporal_Scheme, -4, -4, 4, 4, 100)
        plt.contour( x, y, transpose(rho), linspace(0, 1, 11) )
        plt.axis('equal')
        plt.grid()
        plt.title ('Stability Region of ' + Temporal_Scheme.__name__)
        plt.xlabel('Re')
        plt.ylabel('Im')
        plt.show()
        
Milestone4()
        
