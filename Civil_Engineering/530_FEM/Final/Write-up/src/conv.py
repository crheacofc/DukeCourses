##Convergence study for the stabilized adv-diffusion. Values from MOOSE

import matplotlib.pyplot as plt
import numpy as np
# x by x grids
P_1 = 100 # Reference value for L^2 norms
P_2 = 50
P_3 = 10
P_4 = 20
P_5 = 30
P_6 = 40

E_12 = 0.553022
E_13 = 0.3006194
E_14 = 3.995657*10**-1
E_15 = 0.4858466
E_16 = 0.542456

x = [P_2,P_3,P_4,P_5,P_6]
y = [E_12,E_13,E_14,E_15,E_16]

plt.plot(x,y,'ro')
plt.axis([60,0,0,1.0])
plt.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)),'--') ##plot line best fit

plt.title('L^2 error stabilized advection-diffusion equation')
plt.xlabel('Number of Degrees of Freedom (^2)')
plt.ylabel('L^2 error (compared to DOF=100)')
#plt.show()
plt.savefig('/Users/crhea/Documents/Grad/Duke_courses/Math/530_FEM/Final/Write-up/L^2_error_stabilized.pdf')
plt.close()
##now for the non stabilized
P_1 = 100 # Reference value for L^2 norms
P_2 = 50
P_3 = 10
P_4 = 20
P_5 = 30
P_6 = 40

E_12 = 9.690529e-03
E_13 = 0
E_14 = 4.561869e-02
E_15 = 2.358355e-02
E_16 = 1.405533e-02

x = [P_2,P_4,P_5,P_6]
y = [E_12,E_14,E_15,E_16]

plt.plot(x,y,'ro')
plt.axis([60,0,0,1.0])
plt.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)),'--') ##plot line best fit

plt.title('L^2 error non-stabilized advection-diffusion equation')
plt.xlabel('Number of Degrees of Freedom (^2)')
plt.ylabel('L^2 error (compared to DOF=100)')
#plt.show()
plt.savefig('/Users/crhea/Documents/Grad/Duke_courses/Math/530_FEM/Final/Write-up/L^2_error_nonstabilized.pdf')

