## @module bz.py
# Documentation for bz.py
#
# More information can come here.

#import specpolFlow as pol  #There is an implicit dependency on specpolFlow.iolsd
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import astropy.constants as const
import pandas as pd
import copy
import scipy.special as specialf


def cog_I(lsd, Ic):
    
    nominator = np.trapz(lsd.vel * (Ic-lsd.specI), x=lsd.vel)
    denominator = np.trapz( Ic-lsd.specI, x=lsd.vel )
    return(nominator/denominator)
           
def cog_IV(lsd, Ic):
    nominator = np.trapz(lsd.vel * np.abs( lsd.specV * (Ic-lsd.specI) ), x=lsd.vel )
    denominator = np.trapz( np.abs( lsd.specV * (Ic-lsd.specI) ), x=lsd.vel )
    return(nominator/denominator)

def cog_V(lsd):
    nominator = np.trapz(lsd.vel * np.abs(lsd.specV), x=lsd.vel )
    denominator = np.trapz( np.abs(lsd.specV), x=lsd.vel )
    return(nominator/denominator)

def cog_min(lsd):
    cog_min = lsd.vel[lsd.specI.argmin()]
    if cog_min.size > 1:
        cog_min = cog_min[0]
    return(cog_min)

def FAP(lsd):
    '''Returns the V, null1, and null2 FAP for a given LSD object
    
    If you would like a specific range in velocity, simply slice the LSD object beforehand. 
    Note that the calcBz function also returns the FAP inside the spectral line, 
    over the same velocity range as used in the Bz calculation. 

    :param lsd: lsd object (input)

    :return: FAP V, FAP N1, FAP N2. 
    '''

    #'fitting' the flat line (essentially an average weighted by 1/sigma^2)
    contV = np.sum(lsd.specV/lsd.specSigV**2) / np.sum(1./lsd.specSigV**2)
    contN1 = np.sum(lsd.specN1/lsd.specSigN1**2) / np.sum(1./lsd.specSigN1**2)
    contN2 = np.sum(lsd.specN2/lsd.specSigN2**2) / np.sum(1./lsd.specSigN2**2)
    
    chi2V = np.sum(((lsd.specV - contV)/lsd.specSigV)**2)
    chi2N1 = np.sum(((lsd.specN1 - contN1)/lsd.specSigN1)**2)
    chi2N2 = np.sum(((lsd.specN2 - contN2)/lsd.specSigN2)**2)
   
    approxDOF = (lsd.npix-1.)
    
    probV = 1.0-specialf.gammainc(approxDOF/2., chi2V/2.)
    probN1 = 1.0-specialf.gammainc(approxDOF/2., chi2N1/2.)
    probN2 = 1.0-specialf.gammainc(approxDOF/2., chi2N2/2.)

    return(probV, probN1, probN2)




def calcBz(lsd, cog='I', norm='auto', lambda0=500*u.nm, geff=1.2, velrange=None, bzwidth=None, plot=True):
    '''Calculate the Bz of an LSD profile
    
    :param lsd: lsd object (input). It is assumed that the lsd.vel is in km/s.
    :param cog: choice of value, or calculation method for the center of gravity. The choices are:
                'I': center of gravity of I, 'V': center of gravity of V, 'IV': center of gravity of I*V,
                a float: a user defined value in km/s.
    :param norm: calculation method for the continuum. The choices are:
                    'auto': the median of I outside of velrange (if defined) or the full range (if velrange not defined)
                    float: a used defined value to use for Ic.
    :param lambda0: wavelength of the transition (default=500 nm or 5000 AA).
                    For a LSD profile, this is the lambda value the LSD profile shape was scaled with.
                    Need to be a astropy unit object, with length units.
                    The astropy unit package will take care of the units conversion to give the Bz in Gauss.
    :param geff: effective lande factor of the transition.
                 For an LSD profile, this is the geff value the LSD profile shape was scaled with.
    :param velrange: range of velocity to use for the determination of the
                     line center and the continnum. If not defined, the whole range
                     will be used. If bzwidth is not defined, this range will also be
                     used for the line Bz calculation.
    :param bzwidth: distance from the line center to use in the Bz calculation.
                    One element = same on each side of line center.
                    Two elements, left and right of line center.
                    Not defined: using velrange.
    :param plot: whether or not a graph is displayed.
    :return: a dictionary with Bz and FAP calculations
    '''

    # Velrange is used to identify the position of the line,
    # for calulating the cog, and for calculating the position
    # of the continnum.
    # If Velrange is not defined, it will use the whole range.
    # The range for calculating Bz itself is controlled by bzwidth below
    if velrange != None:
        inside = np.logical_and(lsd.vel>=velrange[0], lsd.vel<=velrange[1])
        lsd_in = lsd[inside]
        lsd_out = lsd[np.logical_not(inside)]
    else:
        lsd_in=copy.copy(lsd)
        
    
    # Check if norm is a string.
    if isinstance(norm, str):
        print('using AUTO method for the normalization')
        if velrange != None:
            print('  using the median of the continuum outside of the line')
            norm_val = np.median(lsd_out.specI)
        else:
            print('  no range in velocity given, using the median of the whole specI to determine continnum')
            norm_val = np.median(lsd_in.specI)
    else:
        norm_val = copy.copy(norm)
        print('using given norm value')
          

    # Check ig cog is a string.
    if isinstance(cog, str):
        # calculate the cog for the chosen method
        if cog == 'I':
            cog_val = cog_I(lsd_in, norm_val)
        if cog == 'min':
            cog_val = cog_min(lsd_in)
        if cog == 'IV':
            cog_val = cog_IV(lsd_in, norm_val)
        if cog == 'V':
            cog_val = cog_V(lsd_in)
    else:
        cog_val=copy.copy(cog)
        

    # Now we define the position of the line for the Bz calculation itself.
    # If the keyword bzwidth is defined, we use that range from the chosen cog.
    if bzwidth == None:
        # No bzwidth. using vrange if defined
        if velrange != None:
            print('no bzwidth defined, using velrange to calculate Bz')
            lsd_bz = copy.copy(lsd_in)
            # saving the range for plotting later.
            p_bzwidth = np.copy(velrange)
        else:
            print('no bzwidth nor velrange defined, using full velocity range to calculate Bz')
            p_bzrange = [lsd.vel.min(), lsd.vel.max()]
            lsd_bz = copy.copy(lsd)
    else:
        # Check whether it is a numpy array
        if isinstance(bzwidth, list) or isinstance(bzwidth, tuple):
            if len(bzwidth) == 1:
                #print('list with one element')
                # keeping the actual bz calculation range for plotting later.
                p_bzwidth = [cog_val-bzwidth, cog_val+bzwidth]
                lsd_bz = lsd[ np.logical_and(lsd.vel >= p_bzwidth[0], lsd.vel <= p_bzwidth[1]) ]
            elif len(bzwidth) == 2:
                #print('list with two elements')
                p_bzwidth = [cog_val-bzwidth[0], cog_val+bzwidth[1]]
                lsd_bz = lsd[ np.logical_and(lsd.vel >= p_bzwidth[0], lsd.vel <= p_bzwidth[1]) ]
            else:
                print('bzwidth has too many elements (need one or two)')
                ## Add an error handling here?
        else:
            p_bzwidth = [cog_val-bzwidth, cog_val+bzwidth]
            lsd_bz = lsd[ np.logical_and(lsd.vel >= p_bzwidth[0], lsd.vel <= p_bzwidth[1]) ]

    # Actual calculation of the Bz:
    
    # u.G is not properly defined in astropy.
    # It is defined in Tesla base units.
    # So here we are creating the right cgs base units for a Gauss.
    G_cgs = u.def_unit('G_cgs', 1 * u.g**0.5/u.cm**0.5/u.s)
    # This is the constant for the Zeeman splitting
    # Lambda_B = constant * lambda0**2 B
    lambda_B_constant = const.e.esu / (4 * np.pi * const.m_e.cgs * const.c.cgs**2)
    #print('Lambda_B constant: ', lambda_B_constant.to(u.AA / (u.AA**2 * G_cgs)))
    
    # set the velocity step for error propagation
    # and set to km/s (because the lsd profile objects don't have units associated with them)
    deltav = (lsd_bz[1].vel - lsd_bz[0].vel)*u.km/u.s # This is in km/s
    
    # Calculation of the integral in the denominator of the Bz function with a trapeze numerical integral
    # For the square error caculation, we propagate like we would for sommation numerical integral.
    ri0v = np.trapz(norm_val-lsd_bz.specI, x=lsd_bz.vel )*u.km/u.s # This is in km/s
    si0v = np.sqrt(np.sum(lsd_bz.specSigI**2 )* deltav**2) # This will naturaly be in km/s

    # Calculation of the integral in the numerator of the Bz function with a trapeze numerical integral
    # For the error caculation, we propagate like we would for sommation numerical integral.
    vf = np.trapz( (lsd_bz.vel - cog_val) * lsd_bz.specV, x=lsd_bz.vel-cog_val )*(u.km/u.s)**2 # This is in (km/s)^2
    svf = np.sqrt(np.sum( (lsd_bz.vel - cog_val )**2 * lsd_bz.specSigV**2 )*(u.km/u.s)**2 * deltav**2)
    
    # Same as above, but for the null profiles.
    nf1 = np.trapz( (lsd_bz.vel - cog_val) * lsd_bz.specN1, x=lsd_bz.vel-cog_val )*(u.km/u.s)**2
    snf1 = np.sqrt(np.sum( (lsd_bz.vel - cog_val )**2 * lsd_bz.specSigN1**2 )*(u.km/u.s)**2* deltav**2)

    nf2 = np.trapz( (lsd_bz.vel - cog_val) * lsd_bz.specN2, x=lsd_bz.vel-cog_val )*(u.km/u.s)**2
    snf2 = np.sqrt(np.sum( (lsd_bz.vel - cog_val )**2 * lsd_bz.specSigN2**2 )*(u.km/u.s)**2* deltav**2)

    # Making the actual Bz calculation.
    # for the units to work out, lambda0 needs to be passed
    # as a unit quantity (e.g. u.nm or u.AA)
    blv = (-1*vf / ( ri0v*geff*lambda0*const.c*lambda_B_constant)).to(G_cgs)
    bln1 = (-1*nf1 / ( ri0v*geff*lambda0*const.c*lambda_B_constant)).to(G_cgs)
    bln2 = (-1*nf2 / ( ri0v*geff*lambda0*const.c*lambda_B_constant)).to(G_cgs)

    # Get the FAP in the same range as the one used for Bz
    FAP_V, FAP_N1, FAP_N2 = FAP(lsd_bz)
        
    result = {
            'Ic': norm_val,
            'cog': cog_val,
            'Bzwidth min': p_bzwidth[0],
            'Bzwidth max': p_bzwidth[1],
            'V bz (G)': blv.value,
            'V bz sig (G)': ( np.abs(blv * np.sqrt( (svf/vf)**2 + (si0v/ri0v)**2 ))).to(G_cgs).value,
            'V FAP': FAP_V,
            'N1 bz (G)': bln1.value,
            'N1 bz sig (G)': ( np.abs(bln1 * np.sqrt( (snf1/nf1)**2 + (si0v/ri0v)**2 ))).to(G_cgs).value,
            'N1 FAP': FAP_N1,
            'N2 bz (G)': bln2.value,
            'N2 bz sig (G)': ( np.abs(bln2 * np.sqrt( (snf2/nf2)**2 + (si0v/ri0v)**2 ))).to(G_cgs).value,
            'FAP_N2': FAP_N2
            }

    df = pd.DataFrame(data=[result])

    if plot:
        fig, ax = lsd.plot(sameYRange=False)
 
        for item in ax:
            if velrange != None:
                item.axvline(x=velrange[0], ls='--', label='velrange')
                item.axvline(x=velrange[1], ls='--')
            item.axvline(x=p_bzwidth[0], ls='dotted', label='bzwidth')
            item.axvline(x=p_bzwidth[1], ls='dotted')
            
        ax[-1].axhline(y=norm_val, ls='--', c='pink', label='Ic')
            
        # for the plot, calculate and display all of the possible methods
        # for calculating the cog.
        for item in ax:
            item.axvline(x=cog_min(lsd_in), label='cog min I', lw=3, alpha=0.5, c='blue')
            item.axvline(x=cog_I(lsd_in, norm_val), label='cog I',lw=3, alpha=0.5, c='red')
            item.axvline(x=cog_IV(lsd_in, norm_val), label='cog I*V',lw=3, alpha=0.5, c='orange')
            item.axvline(x=cog_V(lsd_in), label='cog V',lw=3, alpha=0.5, c='green')
            item.axvline(x=cog_val, label='chosen cog: {}'.format(cog), ls='--', c='k')
           
        ax[-1].legend(loc=0)
        
        red = lsd_bz.vel > cog_val
        blue = lsd_bz.vel < cog_val
        if(len(ax) > 2):
            ax[0].fill_between(lsd_bz.vel[red], lsd_bz.specV[red], step='mid', color='red')
            ax[0].fill_between(lsd_bz.vel[blue], lsd_bz.specV[blue], step='mid', color='blue')
            ax[1].fill_between(lsd_bz.vel[red], lsd_bz.specN1[red], step='mid', color='red')
            ax[1].fill_between(lsd_bz.vel[blue], lsd_bz.specN1[blue], step='mid', color='blue')
        if(len(ax) > 3):
            ax[2].fill_between(lsd_bz.vel[red], lsd_bz.specN2[red], step='mid', color='red')
            ax[2].fill_between(lsd_bz.vel[blue], lsd_bz.specN2[blue], step='mid', color='blue')
        
        return(result,fig)
    else:
        return(result)
