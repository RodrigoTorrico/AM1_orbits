from numpy import array, power

def Kepler(U,t):
    
    x , y, vx, vy = U
    m = (x ** 2 + y ** 2) ** 0.5
    
    if m == 0:
        # Evita la división por cero
        return array([vx, vy, 0, 0])
    
    return array ([vx, vy, -x/m, -y/m])           
