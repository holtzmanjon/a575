import matplotlib.pyplot as plt
import numpy as np
import pdb

def norm() :
    """ calculate normalization for IMF """
    return (100**-1.35 - 0.1**-1.35 ) / -1.35

def func(x) :
    """ return IMF for input mass """
    if x < 0.1 or x > 100 :
        return 0
    else :
        return x**-2.35 / norm()

def intfunc(x) :
    """ return cumulative integral of IMF to input mass """
    return (x**-1.35 - 0.1**-1.35) / -1.35 / norm()

def invint(x) :
    """ return inverse of cumulative integral of IMF """
    return (x*norm()*-1.35 + 0.1**-1.35)**(-1/1.35)

def main(n) :
    """ generate n random samples of IMF """

    # get the IMF values
    dndm=invint(np.random.random(n))

    # bin them manually
    bins=np.arange(0.1,100,0.1)
    data=np.histogram(dndm,bins)
    plt.clf()
    plt.subplot(311)

    # log-log plot
    plt.loglog(data[1][0:-1],data[0])
    plt.xlim([0,10])

    # semilog histogram
    plt.subplot(312)
    plt.hist(dndm,log=True,bins=1000)
    plt.xlim([0,10])

    # linear histogram
    plt.subplot(313)
    plt.hist(dndm,log=False,bins=1000)
    plt.xlim([0,10])

