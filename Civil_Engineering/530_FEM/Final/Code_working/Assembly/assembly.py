'''
This file will assemble the global matrices in two ways:
Method 1 - Brute force and create a large sparse matrix
Method 2 - Create IJV matrix (sparse)
Parameters:
    Con_Mat - connectivity matrix
    Coord_mat - matrix of coordinates
    K_e = elemental stiffnesss matrices
    f_e = elemental f matrices

'''
from Local_Matrices import Dif_basic_K, Dif_basic_F
from read_mesh import read_in
import numpy as np
# Method 1:

def assemble(con_mat,Coord_mat,el_type='Q4',P_type='scalar'):

    if (el_type=='Q4'):
        num_node = 4 #nodes per element
    elif (el_type=='Q8'):
        num_node = 8
    else:
        print("Please enter a valid element type")


    #now lets get the degrees of freedom (i.e. num_dof)
    if (P_type=='scalar'):
        num_dof = len(Coord_mat)
        num_dof_el = num_node #elemental degrees of freedom
        A = con_mat
    elif (P_type=='vector'):
        num_dof = num_node*2 #elemental dof
        num_dof_total = 2*len(Coord_mat)
        A = np.zeros((len(con_mat),num_dof)) #create assembly matrix which is num_elements x num_dof
        for e in range(len(con_mat)): #for each element
            for i in range(num_node): #all the columns of a row
                for j in range(2): #num of dof per node
                    if j==0:
                        A[e,i*2+j] = 2*con_mat[e,i]-1 #assign x first based of nodal id of certain element hence
                                                      # con_mat[e,i]
                    else:
                        A[e,i*2+j] = 2*con_mat[e,i] #assign y from nodal ID

    else:
        print("Please enter a valid problem type")

    K = np.zeros((num_dof,num_dof))  # K matrix will be num_dof x num_dof
    F = np.zeros(num_dof)
    #Loop over each element
    for e in range(len(con_mat)):
        Coord_mat_el = np.zeros((2,num_node)) #create local element's coordinates of nodes; 2 because x,y
        coord_ids = con_mat[e]
        #print(coord_ids) #check to make sure elements are read in right. It works!
        for i  in range(len(coord_ids)):#now we need to get the local coordinates!
            Coord_mat_el[0, i] = Coord_mat[coord_ids[i]][0] #get x-coordinate
            Coord_mat_el[1, i] = Coord_mat[coord_ids[i]][1] #get y-coordinate
            #now that we have coordinates for each node in the element, lets get the K-mat elemental
        #print(e)
        #print(Coord_mat_el.T)

        K_e = Dif_basic_K( Coords=Coord_mat_el.T)
        #print(K_e)

        F_e = Dif_basic_F( lambda x, y: 0, Coord_mat_el.T)
        #print(Coord_mat_el.T) #check that proper values are being returned (they are)

        for i in range(num_dof_el):
            dof_1 = A[e][i]  # get degrees of freedom
            for j in range(num_dof_el):
                dof_2 = A[e][j]
                K[dof_1, dof_2] += K_e[i,j]
            # end j loop
            #print(F_e[i])
            F[dof_1] += F_e[i] #first element is nothing and second is dof
        #end i loop
    #end elemental loop

    file = open('../outputs/K.txt', 'w')
    for i in range(num_dof):
        if i != 0:
            file.write('\n')
        for j in range(num_dof):
            file.write(str(K[i, j]))
            file.write(" ")
    file.close()
    file = open('../outputs/F.txt', 'w')
    for i in range(num_dof):
        if i != 0:
            file.write('\n')
        file.write(str(F[i]))
        file.write(" ")
    file.close()

    return K,F

def assemble_ijv(con_mat,Coord_mat,el_type='Q4',P_type='scalar'):

    if (el_type=='Q4'):
        num_node = 4 #nodes per element
    elif (el_type=='Q8'):
        num_node = 8
    else:
        print("Please enter a valid element type")


    #now lets get the degrees of freedom (i.e. num_dof)
    if (P_type=='scalar'):
        num_dof = len(Coord_mat)
        num_dof_el = num_node #elemental degrees of freedom
        A = con_mat
    elif (P_type=='vector'):
        num_dof = num_node*2
        A = np.zeros((len(con_mat),num_dof)) #create assembly matrix which is num_elements x num_dof
        for e in range(len(con_mat)): #for each element
            for i in range(num_node): #all the columns of a row
                for j in range(2): #num of dof per node
                    if j==0:
                        A[e,i*2+j] = 2*con_mat[e,i]-1 #assign x first based of nodal id of certain element hence
                                                      # con_mat[e,i]
                    else:
                        A[e,i*2+j] = 2*con_mat[e,i] #assign y from nodal ID

    else:
        print("Please enter a valid problem type")

    #K = np.zeros((num_dof,num_dof))  # K matrix will be num_dof x num_dof
    i_array = [] #row
    j_array = [] #col
    v_array = [] #val
    F = np.zeros(num_dof)
    #Loop over each element
    for e in range(len(con_mat)):
        Coord_mat_el = np.zeros((num_node,2)) #create local element's coordinates of nodes; 2 because x,y
        coord_ids = con_mat[e] #access row
        print(coord_ids)
        for i  in range(len(coord_ids)):#now we need to get the local coordinates!
            Coord_mat_el[i, 0] = Coord_mat[coord_ids[i]][0] #get x-coordinate ;
            Coord_mat_el[i ,1] = Coord_mat[coord_ids[i]][1] #get y-coordinate
        print(Coord_mat_el)
            #now that we have coordinates for each node in the element, lets get the K-mat elemental
        K_e = Dif_basic_K(el_type, Coords=Coord_mat_el)
        F_e = Dif_basic_F(el_type,lambda x,y: 10)
        for i in range(num_dof_el):
            dof_1 = A[e][i]  # get degrees of freedom
            for j in range(num_dof_el):
                dof_2 = A[e][j]
                i_array.append(dof_1)
                j_array.append(dof_2)
                v_array.append(K_e[i,j])
                #K[dof_1,dof_2] += K_e[i,j]

            # end j loop
            F[dof_1] = F_e[0,dof_1]
        #end i loop
    #end elemental loop

    return i_array,j_array,v_array,F


'''
coord_mat = np.array([[0,0],[1,0],[2,0],[1,0],[1,1],[1,2]])
#K = assemble(np.array([[1,2,5,4],[2,3,6,5]]),coord_mat)

i_vals,j_vals,v_vals = assemble_ijv(np.array([[0,1,4,3],[1,2,5,4]]),coord_mat)


file = open('../outputs/ijv.txt','w')
count = 0

while (count < len(i_vals)):
    if count != 0:
        file.write('\n')
    file.write(str(i_vals[count])+" "+str(j_vals[count])+" "+str(v_vals[count]))
    count += 1
file.close()


'''
'''
mesh = '../mesh/mesh_nodes_mini.txt'
connect = '../mesh/mesh_connectivity_mini.txt'
sides = '../mesh/sides_mini.txt'
read_in(mesh,connect,sides)

total_nodes = len(mesh)
K = assemble(connect,mesh)
file = open('../outputs/K.txt','w')
for i in range(total_nodes):
    if i != 0:
        file.write('\n')
    for j in range(total_nodes):
        file.write(str(K[i, j]))
        file.write(" ")
file.close()
'''
