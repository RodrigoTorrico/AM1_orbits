from numpy.linalg import eig
from numpy import zeros
from scipy.optimize import fsolve
from Systems_of_Equations.Newton import Jacobian
from Physics.CR3B import CR3B

# Np --> numero de puntos de Lagrange
def Lagrange_points(U_0, Np, mu):
    
    # Matriz para almacenar los puntos de Lagrange
    L_P = zeros([5, 2])

    def F(Y):
        X = zeros(4)
        X[0:2] = Y
        X[2:4] = 0
        FX = CR3B(X, mu)
        return FX[2:4]

    # Calcular puntos de Lagrange usando fsolve
    for i in range(Np):
        L_P[i, :] = fsolve(F, U_0[i, 0:2])

    return L_P


# Verifica la estabilidad de los puntos de Lagrange calculados
def Stability_LP(U_0, mu):

    def F(Y):
        return CR3B(Y, mu)

    # Calcular la matriz jacobiana
    A = Jacobian(F, U_0)
    
    # Calcular valores propios y vectores propios
    valores, vectores = eig(A)

    return valores

def FL(U,t):                 # Wrapped de la función CR3B
    mu = 3.0039e-7
    return CR3B(U,mu)