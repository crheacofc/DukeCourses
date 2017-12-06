'''
Main program for running Galekin Least Squares method for the finite element method.
All elemental choices are kept in the "Elements" folder and are to be imported modularly and they are class based
Mesh numbering starts at 0


'''


import matplotlib.pyplot as plt
import numpy.linalg as npl
#imports from other files
from elements import *
from assembly import assemble, assemble_ijv
from read_mesh import read_in
from EBC import Apply_EBC


# Mesh info read in

mesh = '../mesh/mesh_nodes_mini.txt'
connect = '../mesh/mesh_connectivity_mini.txt'
sides = '../mesh/sides_mini.txt'

mesh,con_mat,top,left,bottom,right = read_in(mesh,connect,sides)

#running of code!

dim = 6

K,F = assemble(con_mat,mesh)
#i_vals,j_vals,v_vals,F = assemble_ijv(con_mat,mesh)
BC_vals = np.zeros(dim-1) #set top EBC to zero
K_EBC,F_EBC = Apply_EBC(K,F,mesh,top[:],BC_vals)

#Basic - very basic - solving
Sol = npl.solve(K_EBC,F_EBC)
#Sol = ssl.spsolve(K_EBC,F_EBC)


#print(Sol)
#now to put back IC to solutions
Sol_final = np.array(np.zeros(len(mesh)))
count = 0
for i in range(len(Sol_final)):
    if (i in top):
        pass #leave as zero like boundary...
    else:
        Sol_final[i] = Sol[count] #next number in solution matrix
        count += 1

#print(Sol_final)
#need to turn solution into matrix now

sol_mat = np.zeros((dim,dim))
count  = 0
for i in range(dim):
    for j in range(dim):
        sol_mat[i][j] = Sol_final[count]
        count +=1

fig,ax=plt.subplots()

x = np.linspace(0.5, -0.5, dim)
y = np.linspace(0.5, -0.5, dim)
X,Y = np.meshgrid(x,y, indexing='xy')
z = plt.contourf(X ,Y , sol_mat)
fig.colorbar(z,ax=ax)
plt.show()
