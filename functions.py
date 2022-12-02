import numpy as np 
from scipy.stats import ortho_group
from functools import partial
from scipy.optimize import minimize


def normal_matrix(n):
    '''
    By generating a random orthonormal matrix and two random diagonal
    matrix creates two normal matrix that commute
    '''
    Q = ortho_group.rvs(dim = n)
    Q_transpose = np.transpose(Q)

    vector_A = np.random.choice([-1,1],n)
    vector_B = np.random.choice([-1,1],n)
   
    diag_A = np.diag(vector_A)
    diag_B = np.diag(vector_B)

    A = Q @ diag_A @ Q_transpose
    B = Q @ diag_B @ Q_transpose

    return A, B, Q


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

def create_e_vector(n, pos):
    e_vector = []
    
    for k in range(n):
        entry = [0]
        if k == pos:
            entry = [1]
       
        e_vector.append(entry)
    
    e_vector = np.array(e_vector)
    return e_vector

def create_Rij(n, i, j, sols):
    # sols = sols.x
    identity = np.identity(n)

    theta = sols[0]
    phi = sols[1]

    c = np.cos(theta)
    c_c = np.conjugate(c)

    s = np.exp(1j*phi)*np.sin(theta)
    s_c = np.conjugate(s)

    e_i = create_e_vector(n, i)
    e_i_t = np.transpose(e_i) 

    e_j = create_e_vector(n, j)
    e_j_t = np.transpose(e_j) 

    R = identity + (c - 1)*(e_i.dot(e_i_t)) - (s_c)*(e_i.dot(e_j_t)) + (s)*(e_j.dot(e_i_t)) + (c_c - 1)*(e_j.dot(e_j_t))

    return R

def objective_function(M, x):
    '''
    Objective Function
    '''
    theta = x[0]
    phi = x[1]

    c = np.cos(theta)
    s = np.exp(1j*phi)*np.sin(theta)

    z = np.array([[np.square(c)],[np.sqrt(2)*c*s],[np.square(s)]])

    return np.linalg.norm(M.dot(z))

def constraint(x):
    '''
    Constraint
    '''
    theta = x[0]
    phi = x[1]

    c = np.cos(theta)
    s = np.exp(1j*phi)*np.sin(theta)

    return np.square(abs(c)) + np.square(abs(s)) - 1

def bnds():
    '''
    Bounds
    '''
    b_1 = (-1*np.pi/4, np.pi/4)
    b_2 = (-1*np.pi, np.pi)

    bounds = (b_1, b_2)

    return bounds
    return None

np.set_printoptions(suppress=True, precision=4)
    
epsilon = 0.01
matrix = normal_matrix(3)
print(matrix[2])
x_0 = (-1*np.pi/4, -1*np.pi)

A = matrix[0]
B = matrix[1]


n = A.shape[0]
Q = np.identity(n)

A_norm = np.linalg.norm(A)
B_norm = np.linalg.norm(B)

off_two = off_2(A, B)
epsilon_norm = epsilon*(A_norm + B_norm)

bounds = bnds()
constraints = {'type':'eq', 'fun':constraint}
cons = [constraints]

while off_two > epsilon_norm:
    for i in range(0, n):
        for j in range(i+1, n):
            M = create_M_ij(A,B,i,j)
            objective = partial(objective_function, M)
            sols = minimize(objective,x_0, method='SLSQP', bounds=bounds, constraints=cons)
            R = create_Rij(n, i, j, [np.pi/10, np.pi/10])
            R = np.asmatrix(R)
            R_H = R.H

            Q = Q @ R
            Q_transpose = np.transpose(Q)


            A = R_H @A@R
            B = R_H @B@R

            diag_A = Q_transpose @ A @ Q
            print(sols.x)
            print(diag_A)

            A_norm = np.linalg.norm(A)
            B_norm = np.linalg.norm(B)

            off_two = off_2(A, B)
            print(off_two)
            epsilon_norm = epsilon*(A_norm + B_norm)


