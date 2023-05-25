## @module cleanmask.py
# Documentation for cleanmask.py
#
# Tools for excluding regions and cleaning the line mask.

import numpy as np
import astropy.units as u
import astropy.constants as const
import specpolFlow.iolsd
import matplotlib.pyplot as plt

def default_exclude_regions(velrange):
    '''
    Defines exclusion regions around H-lines and the telleric regions. 
    
    :param velrange: velocity range around the H-line to be excluded.

    Returns a dictionary with the start and end wavelengths of all of 
    the exclusion regions.

    '''
    velrange=velrange*u.km/u.s
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
    mask_i = specpolFlow.iolsd.mask(fname = name)

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
    print('Masks cleaned!')
    return

def split_order(data):
  '''
  Split an observation object into a list of observation objects with one order per item
  '''
  # one order is where the wavelength backtracks. 
  ind = np.where((data.wl[1:]-data.wl[0:-1]) < 0)[0]
  norder = ind.size+1
  ind = np.append(-1,ind)
  ind = np.append(ind,data.wl.size)
  print('{} orders'.format(norder))

  list_order=[]
  for i in range(0,norder):
    list_order.append(data[ind[i]+1:ind[i+1]])
 
  return(list_order)

def Excluded_Regions_Visual(input_file, input_mask, WLRegions):
    # reading the observed spectrum
    file_obs = input_file
    data_obs = specpolFlow.iolsd.read_spectrum(file_obs)
    # splitting the observed spectrum by order
    list_order = split_order(data_obs)

    # # shifting the model spectrum for its radial velocity
    # # (note the rshift function asks for numpy unit quantities)
    # mod_wave_shift = rshift(mod_wave*u.nm, float(obs["vradCorrected"])*u.km/u.s)

    # read in the mask 
    file_mask=input_mask
    mask = specpolFlow.iolsd.mask(fname=file_mask)

    mask_used = mask.wl[mask.iuse==1]
    mask_not_used = mask.wl[mask.iuse==0]


    for i, order in enumerate(list_order):
        fig, ax = plt.subplots(1,1, figsize=(10,2))
        #Setting limits to axes
        ax.set_xlim([order.wl[0],order.wl[-1]])
        ax.set_ylim(0.8,1.2)
        ax.set_title('Order: {}'.format(i))

        ## Plotting the spectrum
        p = ax.plot(order.wl, order.specI, label='Observation')

        # plotting the lines from the masks
        p = ax.vlines(x=mask_used, ymin=0, ymax=1.2,color='green',linestyle='--',lw=1)
        p = ax.vlines(x=mask_not_used,ymin=0, ymax=np.max(order.specI),color='r',linestyle='--',lw=1)


        ## Plotting the excluded regions in shaded grey
        for i in range(0,np.shape(WLRegions)[0]):
            #if WLRegions.loc[i,'Type']=='Telluric':
                #color = '0.75'
            color = '0.5'
            Xr = np.arange(WLRegions['WLStart'][i],WLRegions['WLFinish'][i],0.01)
            p = ax.fill_between(Xr,y1=0,y2=1.2,facecolor =color)
        return


