'''
Main program for running Galekin Least Squares method for the finite element method.
All elemental choices are kept in the "Elements" folder and are to be imported modularly and they are class based
Mesh numbering starts at 0


'''


import matplotlib.pyplot as plt
import numpy.linalg as npl
#imports from other files
from elements import *
from assembly import assemble
from read_mesh import read_in
from EBC import Apply_EBC


# Mesh info read in

#mesh = '../mesh/mesh_nodes_mini.txt'
#connect = '../mesh/mesh_connectivity_mini.txt'
#sides = '../mesh/sides_mini.txt'

#mesh = '../mesh/mesh_test.txt'
#connect = '../mesh/connect_test.txt'
#sides = '../mesh/sides_test.txt'

mesh = '../mesh/basic.txt'
connect = '../mesh/basic_con.txt'
sides = '../mesh/basic_sides.txt'


mesh,con_mat,top,left,bottom,right = read_in(mesh,connect,sides)
#print(bottom)
#running of code!

dim = 5

K,F = assemble(con_mat,mesh)
#i_vals,j_vals,v_vals,F = assemble_ijv(con_mat,mesh)
#BC_vals = np.zeros(dim) #set top EBC to zero
#BC_vals = np.zeros(8)
#for i in range(len(BC_vals)):
   # BC_vals[i] +=5
#BC_vals = [1]
#boundary = bottom+top
#BC_vals = [10,10,10,0,0,0]
boundary = top+left+bottom+right
BC_vals = np.zeros(len(boundary))

#print(F)
K_EBC,F_EBC,N_wo = Apply_EBC(K,F,mesh,boundary,BC_vals)
print(F_EBC)
#Basic - very basic - solving
Sol = npl.solve(K_EBC,F_EBC)
#Sol = ssl.spsolve(K_EBC,F_EBC)
print(Sol)

#plt.plot(N_wo,Sol)
#plt.show()

Mat_sol = np.zeros((dim,dim))
count1= 0
count2 = 0
BC = 0
for i in range(dim):
    for j in range(dim):
        if count1 in boundary:
            Mat_sol[i][j] = BC_vals[BC]
            BC += 1
        else:
            Mat_sol[i][j] = Sol[count2]
            count2 += 1
        count1 += 1
fig,ax=plt.subplots()

x = np.linspace(0, 1, dim)
y = np.linspace(0, 1, dim)
X,Y = np.meshgrid(x,y, indexing='xy')
z = plt.contourf(X ,Y , Mat_sol)
fig.colorbar(z,ax=ax)
plt.show()
