import numpy as np
from scipy.sparse import csr_matrix
import scipy
import sys
import matplotlib.pyplot as plt
from numpy.random import rand
import scipy.constants 
import scipy.sparse.linalg
from matplotlib import rc
from numpy.random import default_rng
from numpy import random


'''
This python file executes all of the plots and calculations that were asked for in the quiz on Canvas
'''

'''
Question number 2 asked to plot the memory usage of a matrix as function of the gridsize  for dense and CSR storage format
'''

def plot_memory(N=1000, percent_nonzero = 0.2):
    '''
    this function plots the memory usage in MB for a dense vs a csr matrix, filled with random elements where 
    20% of the matrix elements are non-zero.
    '''
    N_ind = []
    mem_csr = []
    mem_dense = []
    for i in range(N):
        array = np.zeros(((i+1)**2))
        array[np.random.choice(np.arange(array.size), replace=False, size=int(array.size * percent_nonzero))] = random.rand()
        array.shape = ((i+1,i+1))
        N_ind.append(i+1)
        mem_dense.append(array.data.nbytes*10**-6)
        mem_csr.append(csr_matrix(array).data.nbytes*10**-6)
    plt.figure()
    plt.plot(N_ind, mem_csr, color='g', label='CSR')
    plt.title('Memory usage comparison for CSR and Dense NxN Matrix \n ${} \%$ of values are non-zero'.format(percent_nonzero*100))
    plt.plot(N_ind, mem_dense, color = 'r', label='Dense')
    plt.legend()
    plt.ylabel('Memory usage / MB')
    plt.xlabel('N')
    plt.savefig('memory_usage.png')
    # plt.show()


'''
Firstly, the function, that calculates the linear index k for given indices i and j of a N x N matrix 
'''

def calc_k(i,j,N):
    '''
    This function evaluates the linear index k for a NxN matrix with given indices i,j
    '''
    if j>=N or i>=N:
        print('Matrix does not have so many indices')
        return None
    return (i*N)+j


'''
Secondly, we have the function that fills a CSR-Matrix with values
The boundary conditions for both cases; Dirichlet and Neuman are already implemented in this function
'''

def populate_csr(fixed_values_dir, fixed_values_neu, N):
    '''
    This function populates a CSR-Matrix with values corresponding to the five-point-stencil.
    The elements will be set to -4 and 1 corresponding to the five-point-stencil.
    First, a lil matrix is created which will be changed to a csr-format in the end.
    
    -----------------------------------------------------------------------------------------

    The boundary conditions are implemented as follows: 

    the function is given two lists as arguments with the corresponding index for the Voltage
    in the 1-dim array, for which the Voltage itself or its derivative is fixed to a certain 
    value. 

    If a value is found in list_dir (corr. to Dirichlet boundary), the corresponding 
    diagonal element is set to 1.

    If a value is found in list_neu (corr. to Neumann boundary), the corresponding 
    diagonal element is set to 1 and one of the neighbouring elements is set to -1. 
    Which of the neighbounring elemenets it is depends on if the boundary is on the right, left
    side or on the top- or bottom side of the grid. 

    '''
    # the given N corresponds to the dimension of the grid (dim NxN), so the Matrix A in Ax=b
    # will have dimension N**2 x N**2
    A = scipy.sparse.lil_matrix((N**2,N**2))
    A.setdiag(-4)                                   # set diagonal elements to -4
    for i in range(N**2):                     
        for j in range(N**2):
            if i == j:
                if i in fixed_values_dir:           # if the Vij corresponding to the diagonal element is fixed to some value
                    A[i,j]   = 1                    # the diagonal element is set to one. None of the other elements in the row are != 0
                elif i in fixed_values_neu:         # if the derivative of Vij is fixed to a certain value 
                                                    # we have to differentiate between four cases: 
                                                    # the derivative is always fixed with respect to the 
                                                    # normal direction on the corresponding side
                    if i in np.arange(0,N):         # up side
                        A[i,j] = 1
                        A[i,j+N] = -1
                    elif i in np.arange(N**2-N, N**2):             # down side
                        A[i,j] = 1
                        A[i,j-N] = -1
                    elif i % N == 0:                # left side 
                        A[i,j+1]=-1
                        A[i,j] = 1
                    elif i % (N-1) == 0:            # right side
                        A[i,j-1]=-1
                        A[i,j] = 1
                else:                               # we have to pay attention to different cases when filling the matrix 
                    if i==0:                        # here, we don't want the two values on the left of the diagonal                                                   
                        A[i,j+1] = 1                  
                        A[i,j+N] = 1
                    elif i+1-N<=0:                  # here, we don't want the leftmost element in the corresponding row
                        A[i,j-1] = 1
                        A[i,j+1] = 1 
                        A[i,j+N] = 1                
                    elif i+1>=N**2:                 # here, we don't want the two rightmost elements in the corresponding row
                        A[i,j-1] = 1
                        A[i,j-N] = 1
                    elif i+N>=N**2:                 # here, we don't want the rightmost element in the corresponding row
                        A[i,j-1] = 1
                        A[i,j+1] = 1
                        A[i,j-N] = 1
                    else:                           # if none of the cases above apply, we fill every element in the corresponding row
                        A[i,j-1] = 1                # with the value 1
                        A[i,j+1] = 1
                        A[i,j-N] = 1  
                        A[i,j+N] = 1
    return A.tocsr()



'''
Thirdly, we had to plot the matrix for a gridsize of N = 16
'''

def plot_matrix(N=16):
    '''
    this function plots the matrix using pyplot.imshow()
    One can also use pyplot.spy(matrix)
    The advantage is not having to convert the matrix to a dense format, but the pyplot.imshow()
    option creates a more informative plot.
    '''
    result = populate_csr([],[],N)
    plt.figure()
    plt.imshow(result.todense(), cmap = 'binary')   
    plt.title('Matrix for laplacian five-point-stencil, N = 16, dim: (N*N x N*N) \nwhite corresponds to an element with value -4, black to 1' )
    plt.ylabel('$N_{y}$')
    plt.xlabel('$N_{x}$')
    plt.savefig('csr_matrix_N_16.png')
    # plt.show()

'''
The next exercise was to solve the laplace equation for a dipole on a grid with charges +1C and -1C 
Therefore, a vector b is created to make it possible to solve the matrix equation of a form like Ax = b
where x is the vector of the Voltage values at each grid points (the NxN matrix was turned into a vector of size N**2)
'''


def b_dipole(N=21):
    '''
    This function creates the vector b for solving 
    Ax = b -> x = A^-1 b
    for the dipole
    A matrix of shape NxN is created and the elements corresponding to the charge locations (10,5) and (10,6) are set 
    to given values of plus/minus 1C
    The function returns the 1-dim array that is to be multiplied to the inverse of the matrix A. 
    '''
    V_mat = np.zeros(shape=(N,N))
    V_mat[9,14] = -1 * (-1) / scipy.constants.epsilon_0 
    V_mat[9,5]  = 1 * (-1) / scipy.constants.epsilon_0
    V_vec = V_mat.flatten()
    return V_vec





def solve_dipole(N=21):
    '''
    This function solves the Laplace equation for the given dipole 
    The tool scipy.sparse.linalg.cg is used, which returns the vector x in the equation:
    Ax = b
    The potential of the dipol is plotted using pyplot.imshow()
    The electric field of the dipole is evaluated by:
    E = - grad(V)
    using np.gradient. It is plotted on top of the potential using pyplot.streamplot().
    '''
    matrix = populate_csr([],[],N)
    vector = b_dipole()
    result = np.reshape(scipy.sparse.linalg.cg(matrix, vector)[0], (N,N))                   # Here, the equation Ax = b is solved using linalg.cg
    print('Error-code dipole: {}'.format(scipy.sparse.linalg.cg(matrix, vector)[1]))        # if scipy.sparse.linalg.cg()[1] != 0 the algorithm did not converge
    e_field_x = -np.gradient(result)[1]                                                     # computing x- and y- component of electric field
    e_field_y = -np.gradient(result)[0]
    x = np.linspace(0,N-1,N)
    y = np.linspace(0,N-1,N)
    X,Y = np.meshgrid(x,y)
    plt.figure()
    plot = plt.imshow(result, cmap = 'seismic')                                     # plotting the potential of the dipole
    plt.streamplot(X,Y,e_field_x, e_field_y, color = 'black', arrowsize=1, arrowstyle='->', linewidth=1)    # plotting the electric field on top
    plt.scatter(14,9, label='Q = -1', color = 'green', s = 10)                      # plotting charge positions 
    plt.scatter(5,9,label='Q = +1', color = 'black', s = 10)
    plt.ylabel('$N_{y}$')
    plt.xlabel('$N_{x}$')
    plt.colorbar(plot, label = 'Voltage / V')
    plt.title('El. field and potential of a Dipole \n electric field lines correspond to the black lines')
    plt.legend()
    plt.savefig('poisson_dipole.png')
    # plt.show()
    return None


'''
The next question was to plot the potential and the electric field for a plate capacitor. 
'''

def b_capacitor(N):
    '''
    This function creates the vector b as in Ax = b for solving the poisson equation 
    for a plate capacitor. 
    In a similar manner as before, a NxN matrix is created and the points on the grid with fixed voltages (on the plates)
    are fixed. Also a list of indices at points where the voltage is fixed, is created. These indices are indices for a 1-dim
    array. This list will be given the function populate_csr() as an argument so that the corresponding matrix elements can be set to one. 
    '''
    a = np.zeros(shape=(N,N))
    a[round(2/5*(N-1)),round(1/3*(N-1)):round(2/3*(N-1))] = 100             # fixing plates to voltages 100 and -100
    a[round(3/5*(N-1)),round(1/3*(N-1)):round(2/3*(N-1))] = -100

    # dir-boundary conditions in list (stores the indices)
    list = []
    a_linear = a.flatten()
    for i in range(len(a_linear)):
        if a_linear[i] != 0:
            list.append(int(i))

    a = a.flatten()
    return list, a                                                          # returning list of dirichlet indices and the vector itself

def solve_capacitor(N = 64):
    '''
    This function, as the one for the dipole solves the poisson equation for the plate capacitor on a grid of 64x64. 
    scipy.sparse.linalg.gmres is used here as a solver because it turned out to work much better than scipy.sparse.linalg.cg
    '''
    list_dir, vector = b_capacitor(N)
    matrix = populate_csr(list_dir, [], N)
    result = scipy.sparse.linalg.gmres(matrix, vector)              # solve equation
    print('Error-code capacitor: {}'.format(result[1]))                       # if >0 it did not converge
    result = np.reshape(result[0], (N,N))                           # bring in shape (N,N) to apply np.gradient()
    e_field_x = -np.gradient(result)[1]                             # and thus compute the electric field 
    e_field_y = -np.gradient(result)[0]

    x = np.linspace(0,N-1,N)                                        # plotting in a similar manner as above
    y = np.linspace(0,N-1,N)
    X,Y = np.meshgrid(x,y)
    plt.figure()
    plot = plt.imshow(result, cmap='seismic')                       
    stream = plt.streamplot(X,Y,e_field_x, e_field_y, color = 'black', arrowsize=1, arrowstyle='->', linewidth=1)
    plt.legend(['Electric field'])
    plt.title('Voltage distribution and electric field for a plate capacitor \n')
    plt.ylabel('$N_{y}$')
    plt.xlabel('$N_{x}$')
    plt.colorbar(plot, label='Voltage / V')
    plt.savefig('poisson_capacitor.png')
    # plt.show()
    return None


######################################################################################################################################################

######################################################################################################################################################

######################################################################################################################################################


def b_capacitor(N):
    '''
    This function creates the vector b as in Ax = b for solving the poisson equation 
    for a plate capacitor. 
    In a similar manner as before, a NxN matrix is created and the points on the grid with fixed voltages (on the plates)
    are fixed. Also a list of indices at points where the voltage is fixed, is created. These indices are indices for a 1-dim
    array. This list will be given the function populate_csr() as an argument so that the corresponding matrix elements can be set to one. 
    '''
    a = np.zeros(shape=(N,N))
    

    # dir-boundary conditions in list (stores the indices)
    list = []
    a_linear = a.flatten()
    for i in range(len(a_linear)):
        if a_linear[i] != 0:
            list.append(int(i))

    a = a.flatten()
    return list, a                                                          # returning list of dirichlet indices and the vector itself

def solve_capacitor(N = 64):
    '''
    This function, as the one for the dipole solves the poisson equation for the plate capacitor on a grid of 64x64. 
    scipy.sparse.linalg.gmres is used here as a solver because it turned out to work much better than scipy.sparse.linalg.cg
    '''
    list_dir, vector = b_capacitor(N)
    matrix = populate_csr(list_dir, [], N)
    result = scipy.sparse.linalg.gmres(matrix, vector)              # solve equation
    print('Error-code capacitor: {}'.format(result[1]))                       # if >0 it did not converge
    result = np.reshape(result[0], (N,N))                           # bring in shape (N,N) to apply np.gradient()
    e_field_x = -np.gradient(result)[1]                             # and thus compute the electric field 
    e_field_y = -np.gradient(result)[0]

    x = np.linspace(0,N-1,N)                                        # plotting in a similar manner as above
    y = np.linspace(0,N-1,N)
    X,Y = np.meshgrid(x,y)
    plt.figure()
    plot = plt.imshow(result, cmap='seismic')                       
    stream = plt.streamplot(X,Y,e_field_x, e_field_y, color = 'black', arrowsize=1, arrowstyle='->', linewidth=1)
    plt.legend(['Electric field'])
    plt.title('Voltage distribution and electric field for a plate capacitor \n')
    plt.ylabel('$N_{y}$')
    plt.xlabel('$N_{x}$')
    plt.colorbar(plot, label='Voltage / V')
    plt.savefig('poisson_capacitor.png')
    # plt.show()
    return None







'''
Now, a function to execute all these functions.
'''

def main():
    '''
    this function executes all of the functions of interest for the assignment. 
    Once called, all plots and calculations are being carried out. 
    '''
    plot_memory()
    plot_matrix()
    solve_dipole()
    solve_capacitor()
    return None

main()























