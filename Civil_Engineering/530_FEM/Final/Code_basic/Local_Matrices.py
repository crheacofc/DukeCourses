'''
The goal of this module is to contain the building blocks for putting together the elemental matrices
The different integrands for the weak form PDE -- henceforth kernels -- will be individual functions to be called in the main program

Parameters:
    el_type = string for the type of element you want to use i.e. 'Q4'
    k_type = type of WF K-Matrix assuming [K]{u}={f} form of equation i.e.'Basic Diffusion'
    f_type = stype of WF F-matrix assuming [K]{u}={f} form of equation i.e.'Basic Diffusion'


'''
import numpy as np
from numpy.linalg import inv

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
    G = np.zeros((2, 4))
    G[0, 0] = -(1. / 4) * (1 - y)
    G[0, 1] = (1. / 4) * (1 - y)
    G[0, 2] = (1. / 4) * (1 + y)
    G[0, 3] = -(1. / 4) * (1 + y)
    G[1, 0] = -(1. / 4) * (1 - x)
    G[1, 1] = -(1. / 4) * (1 + x)
    G[1, 2] = (1. / 4) * (1 + x)
    G[1, 3] = (1. / 4) * (1 - x)
    J = G.dot(C)
    detJ = np.linalg.det(J)
    B = np.dot(inv(J),G)
    return B, detJ
def Dif_basic_K(Coords): #let k be an optional arguement
    num_nodes=4 # number of nodes
    # Now to step through K_mat and integrate each element
    Ke = np.zeros((num_nodes,num_nodes))
    x_ip = [-np.sqrt(3)/3,np.sqrt(3)/3]

    W_ip = [1,1]
    for i in range(len(x_ip)):
        for j in range(len(x_ip)): # now ready to integrate

            B,Jdet = G_Q4_map(x_ip[i],x_ip[j],Coords)
            B_T = np.transpose(B)
            Ke = Ke + np.matmul(B_T,B)* Jdet * W_ip[i] * W_ip[j]
            #print(Ke)
            #Ke = Ke + ((np.transpose(B)) * k * B * Jdet * W_ip[i] * W_ip[j])
            #Ke = Ke + (np.dot(np.transpose(B_val), k*B_val) * abs(J_d_val)*w_ip[i]*w_ip[j])

    '''
    file = open('../outputs/K_elemental.txt', 'w')
    for i in range(4):
        if i != 0:
            file.write('\n')
        for j in range(4):
            file.write(str(Ke[i, j]))
            file.write(" ")
    file.close()
    '''
    return Ke
# now for the {F} kernel for the basic WF diffusion equation

def Dif_basic_F(S,Coords):

    F = np.zeros((len(Coords),1))
    #print(F)
    x_ip = [-np.sqrt(3) / 3, np.sqrt(3) / 3]
    W_ip = [1, 1]
    for i in range(len(x_ip)):
        N = N_Q4(x_ip[i], x_ip[i])
        N_T = np.transpose(N)
        B, Jdet = G_Q4_map(x_ip[i], x_ip[i], Coords)
        #print("N_T is ",N_T)
        S_val = S(x_ip[i],x_ip[i])
        #print(S_val)
        F = F + N_T * S_val * Jdet * W_ip[i] * W_ip[i]
        #print(F)
    #reference file
    '''
    file = open('../outputs/f_elemental.txt', 'w')
    for i in range(4):
        if i != 0:
            file.write('\n')
        file.write(str(F[i]))
        file.write(" ")
    file.close()
    '''
    return F
