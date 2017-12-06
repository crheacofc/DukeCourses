#symbolic integration in several dimensions
from sympy import *
import numpy as np
x=Symbol('x')
y=Symbol('y')
z=Symbol('z')
r=Symbol('r')
t=Symbol('t')
p=Symbol('p')
phi=Symbol('phi')


#double integral!
int1 = (1./7)

#inner integral
i1 = integrate(int1,(y,x**2,x**(1/2)))
print(i1)
#second integral
i2 = integrate(i1,(x,0,1))
print(i2)





## define integrand
int1 = 3

#inner integral
i1 = integrate(int1,(z,0,4))
#print(i1)
#second integral
i2 = integrate(i1,(y,0,1))
#print(i2)
#outermost integral
i3 = integrate(i2,(x,0,2))

#print(i3)

####cylindrical
intp = x*y
intp = intp.subs({x:r*cos(t),y:r*sin(t)})
J = r
intp = intp*J
#inner - z
## No need to do change of coordinates for z -- done automatically
z_in = 0+0*x+0*y
z_out = 2-x**2+y
z_in =z_in.subs({x:r*cos(t),y:r*sin(t)})
z_out =z_out.subs({x:r*cos(t),y:r*sin(t)})
ic1 = integrate(J*intp,(z,z_in,z_out))

#middle - r
r_in = 0
r_out = 2-cos(t)
ic2 = integrate(ic1,(r,r_in,r_out))

#outer - theta
t_in = 0
t_out = 2*np.pi
ic3 = integrate(ic2,(t,t_in,t_out))


####Spherical
intp = 0*x
intp = intp.subs({x:p*sin(phi)*cos(t),y:p*sin(phi)*sin(t),z:p*cos(phi)})
J = p**2*sin(phi)
intp = intp*J
#inner - p
p_in = 0
p_out = 2-cos(t)*sin(phi)
ic1 = integrate(intp,(p,p_in,p_out))

#middle - phi
phi_in = 0
phi_out = 2*np.pi
ic2 = integrate(ic1,(phi,phi_in,phi_out))
#outer - theta
t_in = 0
t_out = 2-cos(phi)
ic3 = integrate(ic2,(t,t_in,t_out))

#print(ic3)
