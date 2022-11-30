import numpy as np 
from scipy.stats import ortho_group

def normal_matrix(n):
    Q = ortho_group.rvs(dim = n)
    Q_inverse = np.linalg.inv(Q)

    vector_A = np.random.choice([-1,1],n)
    vector_B = np.random.choice([-1,1],n)
   
    diag_A = np.diag(vector_A)
    diag_B = np.diag(vector_B)

    A = Q @ diag_A @ Q_inverse
    B = Q @ diag_B @ Q_inverse

    return A, B

normal_matrix = normal_matrix(4)

def off_2(A, B):
    A_array = np.ravel(A)
    A_diag = np.diag(A)
    
    B_array = np.ravel(B)
    B_diag = np.diag(B)

    square_all = np.dot(A_array, A_array) + np.dot(B_array, B_array)
    square_diag = np.dot(A_diag, A_diag) + np.dot(B_diag, B_diag)

    off = square_all - square_diag

    return off


def SimultaneousDiag(normal_matrix, epsilon):
    A = normal_matrix[0]
    B = normal_matrix[1]

    n = A.shape[0]
    Q = np.identity(n)

    A_norm = np.linalg.norm(A)
    B_norm = np.linalg.norm(B)
    off = off_2(A, B)
    epsilon_norm = epsilon*(A_norm + B_norm)

    while off > epsilon_norm:
        for i in range(0, n):
            for j in range(1, n):
                print(off_2)
                print(epsilon_norm)

    

    return None

print(SimultaneousDiag(normal_matrix, 0.1))