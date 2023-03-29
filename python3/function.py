import numpy as np
def funct(x, A, beta, r):
     return A * np.exp(-(x/r)**beta) 
