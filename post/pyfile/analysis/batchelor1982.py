def A11(p, lam):
    return 1 - 60*lam**3/(1+lam)**4/p**4 + 32*lam**3*(15-4*lam**2)/(1+lam)**6/p**6 - 192*lam**3*(5-22*lam**2+3*lam**4)/(1+lam)**8/p**8

def A12(p, lam):
    return 1.5/p - 2*(1+lam**2)/((1+lam)**2*p**3) + 1200*lam**3/(1+lam)**6/p**7

def B11(p, lam):
    return 1 - 68*lam**5/(1+lam)**6/p**6 - 32*lam**3*(10-9*lam**2+9*lam**4)/(1+lam)**8/p**8

def B12(p, lam):
    return 0.75/p + (1+lam**2)/(1+lam)**2/p**3