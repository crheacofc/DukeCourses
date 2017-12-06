#Calculating the determinant of the jacobain of two elements
from sympy import symbols, diff
from sympy.matrices import Matrix
import numpy as np
x, y = symbols('x y')
a , b = symbols('a b')
#first element
N1_1 = (1/(4*a*b))*(a-x)*(b-y)
N1_2 = (1/(4*a*b))*(a+x)*(b-y)
N1_3 = (1/(4*a*b))*(a+x)*(b+y)
N1_4 = (1/(4*a*b))*(a-x)*(b+y)

N1=[N1_1, N1_2, N1_3, N1_4]
Coord_mat1 = Matrix(([-a,-b],[a,-b],[a,b],[-a,b]))

Coord_mat2 = Matrix(([-a,-b],[-a,b],[a,b],[a,-b]))

def Det_jac(N,Coord_mat):
    G = Matrix(([0,0,0,0],[0,0,0,0]))
    for i in range(4):
        G[0,i] = diff(N[i] ,x)
        #print(G[0,i])
        G[1,i] = diff(N[i] ,y)
        #print(G[1,i])
    ##now jacobian
    J = G*Coord_mat
    return J.det()

print(Det_jac(N1,Coord_mat1))
print(Det_jac(N1,Coord_mat2))

