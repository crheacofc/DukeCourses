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

    K = np.zeros((num_dof,num_dof))  # K matrix will be num_dof x num_dof
    F = np.zeros(num_dof)
    #Loop over each element
    for e in range(len(con_mat)):
        Coord_mat_el = np.zeros((2,num_node)) #create local element's coordinates of nodes; 2 because x,y
        coord_ids = con_mat[e]
        print(coord_ids) #check to make sure elements are read in right. It works!
        for i  in range(len(coord_ids)):#now we need to get the local coordinates!
            Coord_mat_el[0, i] = Coord_mat[coord_ids[i]][0] #get x-coordinate
            Coord_mat_el[1, i] = Coord_mat[coord_ids[i]][1] #get y-coordinate
            #now that we have coordinates for each node in the element, lets get the K-mat elemental
        #print(e)
        #print(Coord_mat_el.T)

        K_e = Dif_basic_K( Coords=Coord_mat_el.T)
        #print(K_e)
        def func(x,y):
            if ((x>0.4 and x<0.6)and(y>0.4 and y<0.6)):
                return 10.0
            else:
                return 0.0
        F_e = Dif_basic_F( func, Coord_mat_el.T)
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



    return K,F
