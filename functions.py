import numpy as np 
from scipy.stats import ortho_group
from ortools.linear_solver import pywraplp

def normal_matrix(n):
    Q = ortho_group.rvs(dim = n)
    Q_transpose = np.transpose(Q)

    vector_A = np.random.choice([-1,1],n)
    vector_B = np.random.choice([-1,1],n)
   
    diag_A = np.diag(vector_A)
    diag_B = np.diag(vector_B)

    A = Q @ diag_A @ Q_transpose
    B = Q @ diag_B @ Q_transpose

    return A, B, Q

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

def create_M_ij(A,B,i,j):
    a_ii = A[i,i]
    a_c_ii = np.conjugate(A[i,i])
    a_ij = A[i,j]
    a_c_ij = np.conjugate(A[i,j])
    a_ji = A[j,i]
    a_c_ji = np.conjugate(A[j,i])
    a_jj = A[j,j]
    a_c_jj = np.conjugate(A[j,j])

    b_ii = B[i,i]
    b_c_ii = np.conjugate(B[i,i])
    b_ij = B[i,j]
    b_c_ij = np.conjugate(B[i,j])
    b_ji = B[j,i]
    b_c_ji = np.conjugate(B[j,i])
    b_jj = B[j,j]
    b_c_jj = np.conjugate(B[j,j])

    m_21 = (a_c_ii -a_c_jj)/np.sqrt(2)
    m_22 = (a_ii - a_jj)/np.sqrt(2)
    m_23 = (b_c_ii - b_c_jj)/np.sqrt(2)
    m_24 = (b_ii - b_jj)/np.sqrt(2)

    m_1 = np.array([a_c_ij, a_ji, b_c_ij, b_ji])
    m_2 = np.array([m_21, m_22, m_23, m_24])
    m_3 = np.array([(-1*a_c_ji),(-1*a_ij),(-1*b_c_ji),(-1*b_ij)])

    M = np.matrix([m_1, m_2, m_3])
    M = np.transpose(M)
    
    return M

def minimization(M):
   

    return None
    
def SimultaneousDiag(normal_matrix, epsilon):
    A = normal_matrix[0]
    B = normal_matrix[1]

    n = A.shape[0]
    Q = np.identity(n)

    A_norm = np.linalg.norm(A)
    B_norm = np.linalg.norm(B)

    off_two = off_2(A, B)
    epsilon_norm = epsilon*(A_norm + B_norm)

    while off_two > epsilon_norm:
        for i in range(0, n):
            for j in range(i+1, n):
                M = create_M_ij(A,B,i,j)
                min_problem = minimization(M)

                
                


                
    

    return None

SimultaneousDiag(normal_matrix, 0.01)