'''
This module file will contain the different elements in class formats

Parameters:
    N - shape function matrix
    G - gradient matrix of shape functions
    H - laplacian matrix of shape functions
'''
import numpy as np


##class for element: will have element type. Based off that, we will assign it values for N,B, and H
class Element:
    def __init__(self, type):
        self.type = type

    # done initializing
    # now to set up the N object
    def N(self):
        if self.type == 'Q4':
            #print("N matrix properly loaded")
            N_matrix = lambda x, y: N_Q4(x, y)
            return N_matrix
        else:
            print("Please enter Q4")
            return 0

    def G(self):
        if self.type == 'Q4':
            #print("B matrix properly loaded")
            G_matrix = lambda x, y: G_Q4(x, y)
            return G_matrix
        else:
            print("Please enter Q4")
            return 0

    def G_map(self):
        if self.type == 'Q4':
            # print("B matrix properly loaded")
            G_matrix = lambda x, y: G_Q4_map(x, y)
            return G_matrix
        else:
            print("Please enter Q4")
            return 0
### Shape Functions assuming 1 1 on either side i.e. standard parents elements
#Q4

#B/c of cubit, swap (1 and 3) and (2 and 4)

def N_Q4(x, y):
    N = np.zeros((1, 4))
    N[0, 0] = (1. / 4) * (1 - x) * (1 - y)
    N[0, 1] = (1. / 4) * (1 + x) * (1 - y)
    N[0, 2] = (1. / 4) * (1 + x) * (1 + y)
    N[0, 3] = (1. / 4) * (1 - x) * (1 + y)
    return N
def G_Q4(x, y):
    Q = np.zeros((2, 4))
    Q[0, 0] = -(1. / 4) * (1 - y)
    Q[0, 1] = (1. / 4) * (1 - y)
    Q[0, 2] = (1. / 4) * (1 + y)
    Q[0, 3] = -(1. / 4) * (1 + y)
    Q[1, 0] = -(1. / 4) * (1 - x)
    Q[1, 1] = -(1. / 4) * (1 + x)
    Q[1, 2] = (1. / 4) * (1 + x)
    Q[1, 3] = (1. / 4) * (1 - x)
    return Q


def G_Q4_map(x, y, C):
    Q = np.zeros((2, 4))
    Q[0, 0] = -(1. / 4) * (1 - y)
    Q[0, 1] = (1. / 4) * (1 - y)
    Q[0, 2] = (1. / 4) * (1 + y)
    Q[0, 3] = -(1. / 4) * (1 + y)
    Q[1, 0] = -(1. / 4) * (1 - x)
    Q[1, 1] = -(1. / 4) * (1 + x)
    Q[1, 2] = (1. / 4) * (1 + x)
    Q[1, 3] = (1. / 4) * (1 - x)
    J = Q.dot(C)
    detJ = np.linalg.det(J)
    return Q, detJ
'''
def N_Q4(x, y):
    N = np.zeros((1, 4))
    N[0, 0] = (1./4) * (1 - x) * (1 - y)
    N[0, 1] = (1./4) * (1 + x) * (1 - y)
    N[0, 2] = (1./4) * (1 + x) * (1 + y)
    N[0, 3] = (1./4) * (1 - x) * (1 + y)
    return N

def G_Q4(x, y):
    Q = np.zeros((2, 4))
    Q[0, 0] = -(1. / 4) * (1 - y)
    Q[0, 1] = (1. / 4) * (1 - y)
    Q[0, 2] = (1. / 4) * (1 + y)
    Q[0, 3] = -(1. / 4) * (1 + y)
    Q[1, 0] = -(1. / 4) * (1 - x)
    Q[1, 1] = -(1. / 4) * (1 + x)
    Q[1, 2] = (1. / 4) * (1 + x)
    Q[1, 3] = (1. / 4) * (1 - x)
    return Q
'''

#Q8
def N_Q8(x,y):
    N = np.zeros((1,8))
    N[0, 0] = (1/4.)*(x**2-x)*(y**2-y)
    N[0, 1] = (1/4.)*(x**2+x)*(y**2-y)
    N[0, 2] = (1/4.)*(x**2+x)*(y**2+y)
    N[0, 3] = (1/4.)*(x**2-x)*(y**2+y)
    N[0, 4] = (1/2.)*(1-x**2)*(1+y)
    N[0, 5] = (1/2.)*(1-x)*(1-y**2)
    N[0, 6] = (1/2.)*(1-x**2)*(1-y)
    N[0, 7] = (1/2.)*(1+x)*(1-y**2)


def G_Q8(x,y):
    G = np.zeros((2,8))
    G[0, 0] = (1 / 4.) * (2 * x - 1) * (y** 2 - y)
    G[0, 1] = (1 / 4.) * (2 * x + 1) * (y** 2 - y)
    G[0, 2] = (1 / 4.) * (2 * x + 1) * (y** 2 + y)
    G[0, 3] = (1 / 4.) * (2 * x - 1) * (y** 2 + y)
    G[0, 4] = -x * (1 + y)
    G[0, 5] = -(1 / 2.) * (1 - y ** 2)
    G[0, 6] = -x * (1 - y)
    G[0, 7] = (1 / 2.) * (1 - y ** 2)
    G[1, 0] = (1 / 4.) * (x** 2 - x) * (2 * y - 1)
    G[1, 1] = (1 / 4.) * (x** 2 + x) * (2 * y - 1)
    G[1, 2] = (1 / 4.) * (x** 2 + x) * (2 * y + 1)
    G[1, 3] = (1 / 4.) * (x** 2 - x) * (2 * y + 1)
    G[1, 4] = (1 / 2.) * (1 - x ** 2)
    G[1, 5] = -y * (1 - x)
    G[1, 6] = -(1 / 2.) * (1 - x** 2)
    G[1, 7] = -y * (1 + x)

