import numpy as np 
from scipy.stats import ortho_group

def normal_matrix(n):
    Q = scipy.stats.ortho_group.rvs(dim = n)



    return Q

normal_matrix = normal_matrix(4)
print(normal_matrix)

def SimultaneousDiag(normal_matrix, epsilon):
    A = normal_matrix[0]
    B = normal_matrix[1]
    size = A.size
    print(size)
    Q = np.identity(3)
    
    return None

print(SimultaneousDiag())