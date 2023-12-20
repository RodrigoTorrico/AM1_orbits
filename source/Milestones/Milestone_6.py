##############################################################################################################
##                                             MILESTONE 6                                                  ##
##                                  Lagrange points and their stability                                     ##
##############################################################################################################
##                                                                                                          ##
##                                       Rodrigo Torrico Gijon                                              ##
##############################################################################################################

from numpy import array, zeros, linspace, random, around, size
from ODEs.Cauchy_Problem import Cauchy_Problem
from Physics.CR3B import CR3B
from Physics.Lagrange import Lagrange_points, Stability_LP, FL
from ODEs.Temporal_Schemes import Euler, RK4, Crank_Nicolson, Inverse_Euler, ERK
from Physics.ERK_functions import butcher, StepSize, RK_stages
import matplotlib.pyplot as plt
from random import random

def Milestone6():
    
    T = 500
    dt = 0.01
    n = int(T/dt) + 1
    time_domain = linspace(0,T,n)
    mu = 3.0039e-7

    U0LP = array([[0.8, 0.6, 0, 0],[0.8, -0.6, 0, 0],[-0.1, 0, 0, 0],[0.1, 0, 0, 0],[1.01, 0, 0, 0]])
    Np = 5      # Lagrange Points
    LPAUX = Lagrange_points(U0LP, Np,mu)
    LP = zeros([5,2])
    LP[0,:] = LPAUX[3,:] #Reordeno
    LP[1,:] = LPAUX[4,:] 
    LP[2,:] = LPAUX[2,:] 
    LP[3,:] = LPAUX[0,:] 
    LP[4,:] = LPAUX[1,:] 

    labelPTot = ['L1','L2','L3','L4','L5'] #ordenados
    ShapeLP = ["<",">","d","^","v"]
    ColorLP = ["yellow","cyan","violet","sienna","lightcoral"]
    print(LP)
    for i in range(5):
        plt.plot(LP[i,0],LP[i,1],ShapeLP[i],color = ColorLP[i],label=labelPTot[i])
    plt.plot(-mu, 0, 'o', color = "g", label = 'Tierra')
    plt.plot(1-mu, 0, 'o', color = "b", label = 'Luna')
    plt.grid()
    plt.title("Puntos de Lagrange del CR3BP Tierra-Luna")
    plt.legend(loc = 'upper left',bbox_to_anchor=(1., 0.95))
    plt.savefig('Milestone 6 media/' + 'G1 '+ str(i) +'.png')
    plt.show()


  
    U0_LP_Sel = zeros(4)
    U0_LP_SelStab = zeros(4)
    eps = 1e-4*random()
    
    for k in range(5):

        sel = k + 1 

        if sel == 1:
            labelP = 'L1'
        elif sel == 2:
            labelP = 'L2'
        elif sel == 3:
            labelP = 'L3'
        elif sel == 4:
            labelP = 'L4'
        elif sel == 5:
            labelP = 'L5'
        
        U0_LP_Sel[0:2] = LP[sel-1,:] + eps
        U0_LP_Sel[2:4] = eps

        U0_LP_SelStab[0:2] = LP[sel-1,:]
        U0_LP_SelStab[2:4] = 0

        Autoval_LP = Stability_LP(U0_LP_SelStab, mu) #estabilidad
        print(around(Autoval_LP.real,8))

        for j in range (size(1)):

            U_LP = Cauchy_Problem(time_domain, ERK, FL, U0_LP_Sel)

            #fig, (ax1, ax2) = plt.subplots(1, 2)
            plt.plot(U_LP[:,0], U_LP[:,1],'-',color = "k", label = 'Orbit')
            plt.plot(-mu, 0, 'o', color = "g", label = 'Tierra')
            plt.plot(1-mu, 0, 'o', color = "b", label = 'Luna')
            for i in range(5):
                plt.plot(LP[i,0],LP[i,1],ShapeLP[i],color = ColorLP[i],label=labelPTot[i])
            plt.xlim(-2,2)
            plt.ylim(-2,2)
            plt.title(f"Simulacion CR3BP Tierra-Luna con esquema {ERK.__name__}. Orbita en {labelP}. t = {T}s, dt = {dt}. Vista completa" )    
            plt.legend(loc = 'upper left',bbox_to_anchor=(1., 0.95))
            plt.xlabel('X')
            plt.ylabel('Y')
            plt.savefig('Milestone 6 media/' + 'G2 '+ str(k) + str(j) +'.png')
            plt.grid()
            plt.show()
                    
            plt.plot(U_LP[:,0], U_LP[:,1],'-',color = "k", label = "Orbit" )
            plt.plot(LP[sel - 1,0], LP[sel - 1,1] , ShapeLP[sel-1],color = ColorLP[sel-1], label = labelPTot[sel-1])
            plt.title(f"Simulacion CR3BP Tierra-Luna con esquema {ERK.__name__}. Detalle de orbita en {labelP}. t = {T}s, dt = {dt}" )
            plt.xlabel('X')
            plt.ylabel('Y')
            plt.legend(loc = 'upper right',bbox_to_anchor=(1, 0.5))
            plt.savefig('Milestone 6 media/' + 'G3 '+ str(k) + str(j) +'.png')
            plt.grid()   
            plt.xlim(LP[sel - 1,0]-0.2,LP[sel - 1,0]+0.2)
            plt.ylim(LP[sel - 1,1]-0.2,LP[sel - 1,1]+0.2)
            plt.legend(loc = 'upper left',bbox_to_anchor=(1., 0.95))
            plt.savefig('Milestone 6 media/' + 'G4 '+ str(k) + str(j) +'.png')
            plt.show()

Milestone6()