"""
Routine to use Lane-Emden equations to calculate the structure of a polytrope of index n
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import pdb
import astropy.constants

def le_diffeq(x,y) :
    """
    Calculate derivatives for Lane-Emden equation

    Args:
          x : independent variable xsi
          y : dependent variables: z and y

    Returns:
          derivatives of the variables with respect to x
    """
    z=y[0]
    yy=y[1]
    return np.array([-(yy**n + (2./x)*z), z])
    #return np.array([-(1.0/x**2.0)*(2.0*x*z+(x*2.0)*yy**n), z])

def rk4(x,y,h,func) :
    """
    Runge-Kutte 4th order integrator

    Args:
          x: independent variable
          y: dependent variables (N-dimensional array)
          h: step size
          func: function that evaluats the derivates

    Returns:
          values of dependent variables after one step
    """
    k1 = h * func( x, y)
    k2 = h * func( x+(h/2.), y+(k1/2.))
    k3 = h * func( x+(h/2.), y+(k2/2.))
    k4 = h * func( x+h, y+k3)
    return y + (1./6.)*(k1 + (2.*k2) + (2.*k3) + k4)

ksi = 0.001
n=1.5
h = 0.00001

n=3
ksi = 0.01
h= 0.01

y = np.array([ -1.0/3.0 * ksi + n*4.0/120.0 * ksi**3.0 - 6.0*n*(8.0*n-5.0)/15120.0 * ksi**5,
       1.0 - (ksi**2.0)/6.0 + n/120.0 * ksi**4.0 - n*(8.0*n-5.0)/15120.0 * ksi**6])

ksi_values = ksi
y_values = np.array(y,ndmin=2)

while y[1] >= 0 :
    y = rk4(ksi,y,h,le_diffeq)
    ksi_values = np.append(ksi_values,ksi) 
    y_values = np.append(y_values,[y],axis=0)
    ksi += h      # increment ksi by the step-size value at each iteration
 
theta_n = y_values**n
theta_np1 = y_values**(n+1.0)
q=(y[0]*(ksi_values**2.0))
pdb.set_trace()
q *= (ksi_values[-1]**2 * y_values[-1,0])**(-1)

M = astropy.constants.M_sun.cgs.value
R = astropy.constants.R_sun.cgs.value
G = astropy.constants.G.cgs.value
k = astropy.constants.k_B.cgs.value
rho_c = -(M/((4.0*math.pi/3.0)*(R)**3))*(ksi_values[-1]/3.0)*(1.0/y_values[-1,0])
rho_avg = -3.0/ksi_values[-1]*y_values[-1,0]*rho_c
N_n = ((4*math.pi)**(1.0/n))/(n+1.) * (-ksi_values[-1]**((n+1.)/(n-1.))*y_values[-1,0])**((1.-n)/n)
W_n = 1.0/(4*math.pi*(n+1)*(y_values[-1,0])**2)
P_c = W_n * (G*M**2)/R**4
O_n = 1.0/(-(n+1)*ksi_values[-1]*y_values[-1,0])
T_c = (G*M*0.617*1.67e-24/(k*R))*O_n
print ksi_values[-1],rho_c/rho_avg,W_n,P_c,N_n
print y_values[-1,0]
#print ksi_values[-1],rho_c/rho_avg,N_n,W_n,O_n,round(rho_c,3),round(P_c,3),round(T_c,3)
pdb.set_trace()


fig1 = plt.figure()
ax = plt.subplot(111)
plt.title('Polytrope for n=1.5')
plt.xlabel('ksi')
plt.plot(ksi_values,y_values)
ax.annotate('Theta',xy=(np.median(ksi_values),np.median(y_values)))
plt.plot(ksi_values,theta_n)
ax.annotate('Theta^n',xy=(np.median(ksi_values),np.median(theta_n)))
plt.plot(ksi_values,theta_np1)
ax.annotate('Theta^n+1',xy=(np.median(ksi_values),np.median(theta_np1)))
plt.plot(ksi_values,q)
ax.annotate('q',xy=(np.median(ksi_values),np.median(q)))

plt.show()
