'''
The goal of this module is to contain the building blocks for putting together the elemental matrices
The different integrands for the weak form PDE -- henceforth kernels -- will be individual functions to be called in the main program

Parameters:
    el_type = string for the type of element you want to use i.e. 'Q4'
    k_type = type of WF K-Matrix assuming [K]{u}={f} form of equation i.e.'Basic Diffusion'
    f_type = stype of WF F-matrix assuming [K]{u}={f} form of equation i.e.'Basic Diffusion'


'''
import numpy as np
import scipy.integrate as spi
from elements import Element

#lets make classes for K and F
class WF_Matrix:
    def __init__(self,k_type,f_type,el_type):
        self.k_type = k_type
        self.f_type = f_type
        self.el_type = el_type
    def K(self):
        k=1;
        if self.k_type=='Basic Diffusion':
            K = Dif_basic_K(self.el_type, k)
            return K
        else:
            print("Please try 'Basic Diffusion'")
            return 0

    def F(self):
        k = 1;
        if self.k_type == 'Basic Diffusion':
            S = lambda x, y: x+y**2
            F = Dif_basic_F(self.el_type, S) #Grab f vector
            return F
        else:
            print("Please try 'Basic Diffusion'")
            return 0
# First kernel will just be the basic WF matrix for the left hand side of the diffusion equation!

def Dif_basic_K(el_type,Coords,k=1,P_type = 'scalar'): #let k be an optional arguement
    El = Element(el_type)
    if el_type=='Q4':
        num_nodes=4 # number of nodes
    elif el_type =='Q8':
        num_nodes=8
    else:
        print('Please enter a valid element type')
    K  = np.zeros((num_nodes,num_nodes))
    G = El.G() # snag B-matrix for element
    if (P_type == 'scalar'):  # calculates B function based on scalar or vector problem!
        B = G
    elif (P_type == 'vector'):
        B = np.zeros((3,8))
        for i in range(2*num_nodes):
            B[0,2*i] = G[0,i]
            B[1,2*i-1] = G[1,i]
            B[2,2*i] = G[0,i]
            B[2,2*i-1] = G[1,i]
    else:
        print('Please enter valid problem type')
    #calculate jacobian determinant  -- take another look at this!!!!!

    J = lambda x, y: np.dot(G(x, y),Coords)
    J_d = lambda x, y: np.linalg.det(J(x,y))

   # K_mat = lambda x, y: np.dot(np.transpose(B(x,y)), k*B(x,y))*abs(Jd_m(x,y)) # Create K_matrix based on definition in scalar problem - heat
    #K_mat = lambda x, y: np.dot(K_mat1(x,y), J(x,y)) ##if mapped element, compute jacobian and tack on
    #print(K_mat(0,0)[0,0])
    # Now to step through K_mat and integrate each element
    Ke = np.zeros((len(Coords),len(Coords)))
    x_ip = [-np.sqrt(3)/3,np.sqrt(3)/3]
    w_ip = [1,1]
    for i in range(len(x_ip)):
        for j in range(len(x_ip)): # now ready to integrate
            B_val = G(x_ip[i],x_ip[j])
            J_d_val = J_d(x_ip[i],x_ip[j])
            #print(B_val)
            #print(J_d_val)
            Ke = Ke + (np.dot(np.transpose(B_val), k*B_val) * abs(J_d_val)*w_ip[i]*w_ip[j])
    print(Ke)

    file = open('../outputs/K_elemental.txt', 'w')
    for i in range(4):
        if i != 0:
            file.write('\n')
        for j in range(4):
            file.write(str(K[i, j]))
            file.write(" ")
    file.close()
    return Ke
# now for the {F} kernel for the basic WF diffusion equation

def Dif_basic_F(el_type,S,Coords,P_type='scalar'):
    El = Element('Q4')
    if el_type == 'Q4':
        num_nodes = 4
    elif el_type =='Q8':
        num_nodes=8
    else:
        print('Please enter a valid element type')
    G = El.G()  # snag B-matrix for element
    if (P_type == 'scalar'):  # calculates B function based on scalar or vector problem!
        B = G
    elif (P_type == 'vector'):
        B = np.zeros((3, 8))
        for i in range(2 * num_nodes):
            B[0, 2 * i] = G[0, i]
            B[1, 2 * i - 1] = G[1, i]
            B[2, 2 * i] = G[0, i]
            B[2, 2 * i - 1] = G[1, i]
    else:
        print('Please enter valid problem type')
    F = np.zeros(len(Coords))
    N = El.N()  # snag B-matrix for element
    J = lambda x, y: np.dot(G(x, y), Coords)
    J_d = lambda x, y: np.linalg.det(J(x, y))

    x_ip = [-np.sqrt(3) / 3, np.sqrt(3) / 3]
    w_ip = [1, 1]
    for i in range(len(x_ip)):
        N_val = N(x_ip[i], x_ip[i])
        J_d_val = J_d(x_ip[i], x_ip[i])
        S_val = S(x_ip[i],x_ip[i])
        F = F + np.dot(np.transpose(N_val),S_val*J_d_val*w_ip[i])

    #reference file
    file = open('../outputs/f_elemental.txt', 'w')
    for i in range(4):
        if i != 0:
            file.write('\n')
        file.write(str(F[0, i]))
        file.write(" ")
    file.close()
    return F


'''
Kl = Dif_basic_K('Q4',Coords = [[-1, -1], [1, -1], [1, 1], [-1, 1]])
Fl = Dif_basic_F('Q4', lambda x, y: x+y)
file = open('../outputs/K_elemental.txt','w')
for i in range(4):
    if i != 0:
        file.write('\n')
    for j in range(4):
        file.write(str(Kl[i, j]))
        file.write(" ")
file.close()
'''
