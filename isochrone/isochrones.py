import os
import numpy as np
from astropy.io import ascii
import pdb
import astropy
import math

os.environ['ISOCHRONE_DIR'] = '/Users/holtz/apogee/dist/isochrones/'
def basicread(infile) :
    """
    Routine to read a Padova isochrone file using low-level I/O and return
    a numpy structured array with the contents

    Args: 
        file name : input data file name

    Returns: 
        structured array with isochrone data
    """

    # open the file
    file = open(infile)

    # initialize the lists to hold the file contents
    z = [] ; age = [] ; mini = []; mact = []; logl = []; logte = []; logg = []
    mbol = []; u = []; b = []; v = []; r = []; i = []; j = []; h = []; k = [] 
    intimf = []; stage = []

    # loop through lines in file
    nlines = 0
    for line in file :
       # ignoring comment lines starting with #, add input data into list vars
       if line.startswith('#') is False :
           cols=line.split()
           z.append(cols[0])
           age.append(cols[1])
           mini.append(cols[2])
           mact.append(cols[3])
           logl.append(cols[4])
           logte.append(cols[5])
           logg.append(cols[6])
           mbol.append(cols[7])
           u.append(cols[8])
           b.append(cols[9])
           v.append(cols[10])
           r.append(cols[11])
           i.append(cols[12])
           j.append(cols[13])
           h.append(cols[14])
           k.append(cols[15])
           intimf.append(cols[16])
           stage.append(cols[17])
           nlines+=1

    # declare the output structured array
    data = np.recarray(nlines,dtype=[
                       ('z','f4'),('age','f4'),
                       ('mini','f4'),('mact','f4'),
                       ('logl','f4'),('logte','f4'),('logg','f4'),
                       ('mbol','f4'),('u','f4'),('b','f4'),
                       ('v','f4'),('r','f4'),('i','f4'),
                       ('j','f4'),('h','f4'),('k','f4'),
                       ('intimf','f4'),('stage','i4')
                       ])

    # fill the contents and return
    data['z'] = z
    data['age'] = age
    data['mini'] = mini
    data['mact'] = mact
    data['logl'] = logl
    data['logte'] = logte
    data['logg'] = logg
    data['mbol'] = mbol
    data['u'] = u
    data['b'] = b
    data['v'] = v
    data['r'] = r
    data['i'] = i
    data['j'] = j
    data['h'] = h
    data['k'] = k
    data['intimf'] = intimf
    data['stage'] = stage
    return data

def basicread2(infile) :
    """
    Routine to read a Padova isochrone file using low-level I/O and return
    a numpy structured array with the contents

    Args: 
        file name : input data file name

    Returns: 
        structured array with isochrone data
    """

    # open the file
    file = open(infile)

    # setup up the names and data types of the columns, and initialize 
    #   list of lists
    ncols= 18
    names = [('z','f4'),('age','f4'),
            ('mini','f4'),('mact','f4'),
            ('logl','f4'),('logte','f4'),('logg','f4'),
            ('mbol','f4'),('u','f4'),('b','f4'),
            ('v','f4'),('r','f4'),('i','f4'),
            ('j','f4'),('h','f4'),('k','f4'),
            ('intimf','f4'),('stage','i4')]
    listdata=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

    # loop through lines, ignoring comments, filling listdata
    nlines = 0
    for line in file :
       # ignoring comment lines starting with #, add input data into list vars
       if line.startswith('#') is False :
          cols=line.split()
          for i in range(ncols) :
              listdata[i].append(cols[i])
          nlines+=1

    # define the numpy structrued array and fill the columns
    data = np.recarray(nlines,dtype=names)
    for i in range(ncols) :
       #pdb.set_trace()
       data[names[i][0]] = listdata[i]

    return data


def basicread3(infile) :
    """
    Routine to read a Padova isochrone file using low-level I/O and return
    a numpy structured array with the contents

    Args: 
        file name : input data file name

    Returns: 
        structured array with isochrone data
    """

    #use astropy.io.ascii.read!
    data=ascii.read(infile,names=['z','age','mini','mact','logl','logte','logg','mbol','u','b','v','r','i','j','h','k','intimf','stage'])

    return data


def read(infile,columns=None,age=None) :
    """
    Routine to read a Padova isochrone file using low-level I/O and return
    a numpy structured array with the contents

    Args: 
        file name : input data file name
        columns=[list] : list of columns to extract
        age = age : single age to extract

    Returns: 
        structured array with isochrone data
    """

    if os.getenv('ISOCHRONE_DIR') != "" :
        data=ascii.read(os.getenv('ISOCHRONE_DIR')+'/'+infile,
             names=['z','age','mini','mact','logl','logte','logg',
                    'mbol','u','b','v','r','i','j','h','k','intimf','stage'])
    else :
        data=ascii.read(infile,
             names=['z','age','mini','mact','logl','logte','logg',
                    'mbol','u','b','v','r','i','j','h','k','intimf','stage'])

    # option to select a certain age
    if age is not None:
        gd=np.where(data['age']==age)
        data=data[gd]

    # default columns
    if columns is None:
        # can set default columns here, or keep all quantitites
        #data.keep_columns(['Z','age','logte','logl','intimf','stage'])
        pass
    # option to extract specified columns
    else:
        data.keep_columns(columns)

    return data

def radius(logl,logte) :
    """ 
    Get stellar radius given luminosity, effective temperature
    """
    teff=10.**logte * astropy.units.K
    lum = 10.**logl * astropy.units.Lsun * astropy.constants.L_sun.cgs
    return np.sqrt(lum / (4.*math.pi*astropy.constants.sigma_sb.cgs*teff**4.))

def mkhess(age, file='zp00.dat', xval='logte',yval='logl',xr=[3.0,4.0],yr=[-6.,5],nbins=200) :
    """
    Routine to make a Hess diagram of a particular age from a file
    """
    # read the isochrone
    iso = read(file,age=age)

    # initialize Hess diagram and set bin sizes given limits and nbins
    hess = np.zeros([nbins,nbins])
    dx=(xr[1]-xr[0])/nbins
    dy=(yr[1]-yr[0])/nbins

    # calculate bin locations
    xbin = ((iso[xval]-xr[0])/dx).astype('int')
    ybin = ((iso[yval]-yr[0])/dy).astype('int')

    # loop over each pair of isochrone points
    for i in range(len(xbin)-1) :
        # number of stars in between these points
        nimf = iso['intimf'][i+1]-iso['intimf'][i]

        # get min and max bin numbers
        xmin= np.min(xbin[i:i+2])
        xmax= np.max(xbin[i:i+2])
        ymin= np.min(ybin[i:i+2])
        ymax= np.max(ybin[i:i+2])

        # make sure that bin range is within the output Hess diagram
        if ymax > 0 and xmax > 0 and ymin < nbins-1 and xmin < nbins-1 :
            # number of bins over which stars are spread
            nbin = (xmax-xmin+1)*(ymax-ymin+1)

            # make sure we don't go over the edge of the Hess diagram
            xmin= np.max([xmin,0])
            xmax= np.min([xmax,nbins-1])
            ymin= np.max([ymin,0])
            ymax= np.min([ymax,nbins-1])
 
            # add the stars in!
            hess[ymin:ymax+1,xmin:xmax+1] += nimf / nbin

    return hess
