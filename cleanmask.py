## @module cleanmask.py
# Documentation for cleanmask.py
#
# Tools for excluding regions and cleaning the line mask.

import numpy as np
import astropy.units as u
import astropy.constants as const
import specpolFlow as pol


def default_exclude_regions(velrange):
    '''
    Defines exclusion regions around H-lines and the telleric regions. 
    
    :param velrange: velocity range around the H-line to be excluded.

    Returns a dictionary with the start and end wavelengths of all of 
    the exclusion regions.

    '''

    Halpha = 656.3
    Hbeta = 486.14
    Hgamma = 434.05
    Hdelta = 410.17
    Hepsilon = 397.01
    lines = np.array([Halpha,Hbeta,Hdelta,Hepsilon])*u.nm

    Hjump = np.array([360,392])*u.nm
    tellstart = np.array([587.5,627.5,684.0,717.0,757.0,790.0,809.0])   # *u.nm
    tellend = np.array([592.0,632.5,705.3,735.0,771.0,795.0,990.0])     # *u.nm

    jumpend = (velrange/const.c*Hjump+Hjump).to(u.nm)/u.nm
    linestart = (-velrange/const.c*lines+lines).to(u.nm)/u.nm
    lineend = (velrange/const.c*lines+lines).to(u.nm)/u.nm
    linestart = np.append(linestart,jumpend[0])
    lineend = np.append(lineend,jumpend[1])

    #data = pd.DataFrame(data = {"WLStart": np.concatenate((tellstart,linestart)), 'WLFinish': np.concatenate((tellend,lineend))})
    data = {
        "WLStart": np.concatenate((tellstart,linestart)),
        'WLFinish': np.concatenate((tellend,lineend))
    }

    return data


def clean_model_mask(name_in, name_out, data):
    '''
    Function that cleans a given mask.

    :param name_in: name of input uncleaned mask
    :param name_out: name of the output cleaned mask
    :param data: python dictionary containing the start and stop wavelength of regions to be excluded
    
    '''
    # Place the regions into a matrix
    L = data['WLStart'].size

    regions = np.zeros((L,2))

    regions[:,0] = data['WLStart'][:]
    regions[:,1] = data['WLFinish'][:]


    name = name_in
    mask_i = pol.iolsd.mask(fname = name)

    #Identifying if a line is inside or outside the regions
    dcard = []
    for lines in mask_i.wl:
      for j in range(0,np.size(regions[:,0])):
        if regions[j,0] < lines < regions[j,1]:
          dcard.append(lines)

    #Creating a list with the mask.wl elements
    mask_list = mask_i.wl.tolist()

    for line in dcard:
      index = (mask_list).index(line)
      mask_i.iuse[index] = 0
    
    l = np.size(mask_i.wl)
    spec_lines = []
    for j in range(0,l):
      if mask_i.iuse[j] == 1:
        spec_lines.append(mask_i.wl[j])
    mask_i.save(name_out)
    return()
