##Homework 7 for FEM code for problem 1!
import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
def function(x):
    return np.e**(2*x)

def Gauss_quad(func,num_points,points,weights):
    val = 0
    for i in range(num_points):
        val+=weights[i]*func(points[i])
    return val
G_vals=np.zeros(4)
points1 =[0]
weights1 = [2]
G_1 = Gauss_quad(function,1,points1,weights1)

points2 =[-np.sqrt(3)/3,np.sqrt(3)/3]
weights2 = [1,1]
G_2 = Gauss_quad(function,2,points2,weights2)

points3 =[-np.sqrt(3/5),0,np.sqrt(3/5)]
weights3 = [5/9,8/9,5/9]
G_3 = Gauss_quad(function,3,points3,weights3)

points5 =[-(1/3)*np.sqrt(5+2*np.sqrt(10/7)),-(1/3)*np.sqrt(5-2*np.sqrt(10/7)),0,(1/3)*np.sqrt(5-2*np.sqrt(10/7)),(1/3)*np.sqrt(5+2*np.sqrt(10/7))]
weights5 = [(322-13*np.sqrt(70))/900,(322+13*np.sqrt(70))/900,128/225,(322+13*np.sqrt(70))/900,(322-13*np.sqrt(70))/900]
G_5 = Gauss_quad(function,5,points5,weights5)

true_val=spi.quad(function,-1,1)[0]

G_mat = [G_1,G_2,G_3,G_5]
T_mat = np.zeros(4)
for i in range(4):
    T_mat[i] = true_val

abs_err = np.abs(T_mat-G_mat)

plt.plot([1,2,3,5],abs_err)
plt.title('Error Estimate of Gaussian Quadrature')
plt.xlabel('Number of Quad Points')
plt.ylabel('Absolute Error')
#plt.show()
plt.savefig('/Users/crhea/Documents/Grad/Math/530_FEM/HW7/problem1.png')
