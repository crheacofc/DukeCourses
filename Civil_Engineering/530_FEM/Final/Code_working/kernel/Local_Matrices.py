'''
The goal of this module is to contain the building blocks for putting together the elemental matrices
The different integrands for the weak form PDE -- henceforth kernels -- will be individual functions to be called in the main program

Parameters:
    el_type = string for the type of element you want to use i.e. 'Q4'
    k_type = type of WF K-Matrix assuming [K]{u}={f} form of equation i.e.'Basic Diffusion'
    f_type = stype of WF F-matrix assuming [K]{u}={f} form of equation i.e.'Basic Diffusion'


'''
import numpy as np
from elements import Element



def Dif_basic_K(Coords,k=1): #let k be an optional arguement
    num_nodes=4 # number of nodes
    E = Element("Q4")

   # K  = np.zeros((num_nodes,num_nodes))
    # Now to step through K_mat and integrate each element
    Ke = np.zeros((len(Coords),len(Coords)))

    x_ip = [-np.sqrt(3)/3,np.sqrt(3)/3]

    W_ip = [1,1]
    for i in range(len(x_ip)):
        for j in range(len(x_ip)): # now ready to integrate
            G_m = E.G_map()
            B,Jdet = G_m(x_ip[i],x_ip[j],Coords)
            B_T = np.transpose(B)
            Ke = Ke + B_T.dot(B)*k* Jdet * W_ip[i] * W_ip[j]
            #print(Ke)
            #Ke = Ke + ((np.transpose(B)) * k * B * Jdet * W_ip[i] * W_ip[j])
            #Ke = Ke + (np.dot(np.transpose(B_val), k*B_val) * abs(J_d_val)*w_ip[i]*w_ip[j])

    #file = open('../outputs/K_elemental.txt', 'w')
    #for i in range(4):
    #    if i != 0:
    #        file.write('\n')
    #    for j in range(4):
    #        file.write(str(K[i, j]))
    #        file.write(" ")
    #file.close()
    return Ke

# now for the {F} kernel for the basic WF diffusion equation

def Dif_basic_F(S,Coords):
    E = Element("Q4")
    N = E.N()
    G_map = E.G_map()
    F = np.zeros((len(Coords),1))
    #print(F)
    x_ip = [-np.sqrt(3) / 3, np.sqrt(3) / 3]
    W_ip = [1, 1]
    for i in range(len(x_ip)):

        N_l = N(x_ip[i], x_ip[i]) #N local
        N_T = np.transpose(N_l)
        B, Jdet = G_map(x_ip[i], x_ip[i], Coords)
        #print("N_T is ",N_T)
        S_val = S(x_ip[i],x_ip[i])
        #print(S_val)
        F = F + N_T * S_val * Jdet * W_ip[i] * W_ip[i]
        #print(F)
    #reference file
    file = open('../outputs/f_elemental.txt', 'w')
    for i in range(4):
        if i != 0:
            file.write('\n')
        file.write(str(F[i]))
        file.write(" ")
    file.close()
    return F

def ES_K(Coords):
    num_nodes = 4  # number of nodes
    E = Element("Q4")

    # Now to step through K_mat and integrate each element
    Ke = np.zeros((len(Coords), len(Coords)))
    x_ip = [-np.sqrt(3) / 3, np.sqrt(3) / 3]
    print(len(Coords))
    W_ip = [1, 1]
    for i in range(len(x_ip)):
        for j in range(len(x_ip)):  # now ready to integrate
            G_m = E.G_map()
            B, Jdet = G_m(x_ip[i], x_ip[j], Coords)