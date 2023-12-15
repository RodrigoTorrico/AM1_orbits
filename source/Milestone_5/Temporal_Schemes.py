from scipy.optimize import newton

def Euler(U, t1, t2, F):
    
    return U + (t2 - t1) * F(U, t2)

def Crank_Nicolson(U, t1, t2, F):
    
    def ResidualCN(X):
        
        return X - U - (t2 -t1)/2 * (F(U, t1) + F(U, t2)) - (t2 -t1)/2 * (F(X, t1) + F(X, t2 + (t2 - t1)))
    
    return newton(ResidualCN, U)


def Inverse_Euler(U, t1, t2, F):

    def ResidualIE(G):
        
        return G - U - (t2 - t1) *  F(U, t2)

    return newton(func = ResidualIE, x0 = U)

def RK4(U, t1, t2, F):
    
    k1 = F(U,t2)
    k2 = F(U + (t2 - t1) * k1/2, t2 + (t2 -t1)/2)
    k3 = F(U + (t2 - t1) * k2/2, t2 + (t2 -t1)/2)
    k4 = F(U + (t2 - t1) * k3, t2 + (t2 - t1))
    
    return U + (t2 - t1) * (k1 + 2*k2 + 2*k3 + k4)/6
    
