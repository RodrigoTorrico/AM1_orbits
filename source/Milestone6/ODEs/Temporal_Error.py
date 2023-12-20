
from ODEs.Cauchy_Problem import Cauchy_Problem
from numpy import array, zeros, log10, vstack, ones
from numpy.linalg import norm, lstsq

def Error_Richardson(time_domain, temporal_scheme, F, U0, q):
    
    T = len(time_domain)
    Nv = len(U0)
    t2 = zeros(2*T)
    
    for i in range(T - 1):
        t2 [2*i] = time_domain [i]
        t2 [2*i + 1] = (time_domain[i] + time_domain[i + 1])/2
    t2 [2*(T - 1)] = time_domain[T - 1] 
     
    U = Cauchy_Problem(time_domain, temporal_scheme, F, U0)
    # Malla fina
    V = Cauchy_Problem(t2, temporal_scheme, F, U0)
    
    # Inicializamos error
    E = zeros((T, Nv))
    for i in range (T):
        E[i, :] = (V[i, :] - U[i, :])/(1 - 1./2**q)
        
    # Inicializamos solucion
    S = zeros((T, Nv))
    S = U + E
   
    return E, S

def convergence_rate(time_domain, temporal_scheme, F, U0,n):
    
    log_E = zeros(n)
    log_n = zeros(n)
    t1 = time_domain
    
    U = Cauchy_Problem(t1, temporal_scheme, F, U0)
    
    for i in range(n):
        T = len(t1)
        t2 = zeros(2*T - 1)

        for j in range(T - 1):
            t2[2*j] = t1[j]
            t2[2*j + 1] = (t1[j] + t1[j + 1]) / 2

        t2[2*(T - 1)] = t1[T - 1]
        
        V = Cauchy_Problem(t2, temporal_scheme, F, U0)
        
        error = norm(V[(T - 1)*2, :] - U[int(((T - 1)*2)/2), :])
        log_E[i] = log10(error)
        log_n[i] = log10(len(t2))
        t1 = t2
        U = V
        
    for j in range(n):
        
        if abs(log_E[j]) > 12:
            break 
        
        j = min(j, n - 1)
        x = log_n[0:j + 1]
        y = log_E[0:j + 1]
        A = vstack( [ x, ones(len(x)) ] ).T
        n_temp, c = lstsq(A, y, rcond=None)[0]
        order = abs(n_temp) 
        log_E2 = log_E - log10( 1 - 1./2**order)
    print(log_E)
    print(log_E2)
    print(log_n)
    return order, log_E, log_E2, log_n