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


def le_solve(xi0=0.001, n=1.5, h=0.001) :
    """
    Routine to solve Lane-Emden equation for given polytropic index

    Args:

    Returns:
         2 arrays with  z=y' and y=....
    """
    xi = xi0

    #set boundary condition
    y = np.array([ -1.0/3.0 * xi + n*4.0/120.0 * xi**3.0 - 6.0*n*(8.0*n-5.0)/15120.0 * xi**5,
       1.0 - (xi**2.0)/6.0 + n/120.0 * xi**4.0 - n*(8.0*n-5.0)/15120.0 * xi**6])

    # initialize arrays to save quantities (really should not be needed! use
    # xi and y themselves?
    xi_values = xi
    y_values = np.array(y,ndmin=2) # 2D array with [Nsteps,2] for two vars

    # integrate until we hit surface
    while y[1] >= 0 :
        y = rk4(xi,y,h,le_diffeq)
        xi_values = np.append(xi_values,xi) 
        y_values = np.append(y_values,[y],axis=0)
        xi += h      # increment xi by the step-size value at each iteration

    return xi_values, y_values

def le_values(xi, y) : 
    """    
    Routine to calculate various variables given polytropic structure

    Args:
          xi, y arrays, with y[0]=y', y[1]=y
    """

    theta_n = y[:,1]**n
    theta_np1 = y[:,1]**(n+1.0)
    q=(y[0]*(xi**2.0))
    pdb.set_trace()
    q *= (xi[-1]**2 * y[-1,0])**(-1)

    M = astropy.constants.M_sun.cgs.value
    R = astropy.constants.R_sun.cgs.value
    G = astropy.constants.G.cgs.value
    k = astropy.constants.k_B.cgs.value
    rho_c = -(M/((4.0*math.pi/3.0)*(R)**3))*(xi[-1]/3.0)*(1.0/y[-1,0])
    rho_avg = -3.0/xi[-1]*y[-1,0]*rho_c
    N_n = ((4*math.pi)**(1.0/n))/(n+1.) * (-xi[-1]**((n+1.)/(n-1.))*y[-1,0])**((1.-n)/n)
    W_n = 1.0/(4*math.pi*(n+1)*(y[-1,0])**2)
    P_c = W_n * (G*M**2)/R**4
    O_n = 1.0/(-(n+1)*xi[-1]*y[-1,0])
    T_c = (G*M*0.617*1.67e-24/(k*R))*O_n
    print xi[-1],rho_c/rho_avg,W_n,P_c,N_n


#def le_plot(xi, y, erase=True,....) :
