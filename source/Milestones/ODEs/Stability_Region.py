from numpy import linspace, zeros

def Stability_Region(Temporal_Scheme, x0, y0, xf, yf, N):
    x, y = linspace(x0, xf, N), linspace (y0, yf, N)
    rho = zeros((N, N))
    for i in range(N): 
      for j in range(N):

          w = complex(x[i], y[j])
          r = Temporal_Scheme( 1., 1., 0., lambda u, t: w*u )
          rho[i, j] = abs(r) 

    return rho, x, y
