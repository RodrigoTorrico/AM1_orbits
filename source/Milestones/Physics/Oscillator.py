from numpy import zeros, array


def Oscillator(U,t):
    return array([U[1], -U[0]])

def System(U0, F, t):
    eps = 1e-6
    N = len(U0)
    M =  zeros( (N, N))
    delta = zeros(N) 
     
    for j in range(N):  
         
        delta[:] = 0 
        delta[j] = eps 
        M[:, j] = ( F( U0 + delta, t ) - F( U0 - delta, t ) )/(2*eps)
  
    return M

def test_system_matrix(): 

    U0 = array( [ 0, 0 ] ) 
    t = 0
    M = System(U0, Oscillator, t)
    print("M=", M)

test_system_matrix()