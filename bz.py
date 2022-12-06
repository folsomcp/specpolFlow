#!/usr/bin/python3
## @module bz.py
# Documentation for bz.py
#
# Calculate a longitudinal magnetic field (Bz) from an LSD profile.

#import specpolFlow as pol  #implicitly depends on specpolFlow.iolsd
#import matplotlib.pyplot as plt #Implicitly depends on matplotlib.pyplot
import numpy as np
import astropy.units as u
import astropy.constants as const
import copy
import scipy.special as specialf

#Routines for calculating the center of gravity of a line/LSD profile
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

    The False Alarm Probability (FAP) is the probability that the observed
    data are consistent with the null hypothesis of no magnetic field.
    In this case, the null hypothesis is a flat line in Stokes V
    (an offset from 0 is allowed to account for continuum polarization).
    The probability is evaluated from chi^2.
    
    If you would like a specific range in velocity, simply slice the LSD object beforehand. 
    Note that the calcBz function also returns the FAP inside the spectral line, 
    over the same velocity range as used in the Bz calculation. 

    :param lsd: lsd object (input)

    :return: FAP V, FAP N1, FAP N2. 
    '''

    approxDOF = (lsd.npix-1.)

    #'fitting' the flat line (essentially an average weighted by 1/sigma^2)
    contV = np.sum(lsd.specV/lsd.specSigV**2) / np.sum(1./lsd.specSigV**2)
    chi2V = np.sum(((lsd.specV - contV)/lsd.specSigV)**2)
    probV = 1.0-specialf.gammainc(approxDOF/2., chi2V/2.)
    #repeat for the Null1 and Null2 profiles (if they exist)
    probN1 = 0.
    probN2 = 0.
    
    if lsd.numParam > 2:
        contN1 = np.sum(lsd.specN1/lsd.specSigN1**2) / np.sum(1./lsd.specSigN1**2)
        chi2N1 = np.sum(((lsd.specN1 - contN1)/lsd.specSigN1)**2)
        probN1 = 1.0-specialf.gammainc(approxDOF/2., chi2N1/2.)
    if lsd.numParam > 3:
        contN2 = np.sum(lsd.specN2/lsd.specSigN2**2) / np.sum(1./lsd.specSigN2**2)
        chi2N2 = np.sum(((lsd.specN2 - contN2)/lsd.specSigN2)**2)
        probN2 = 1.0-specialf.gammainc(approxDOF/2., chi2N2/2.)
    
    return(probV, probN1, probN2)


def calcBz(lsd, cog='I', norm='auto', lambda0=500*u.nm, geff=1.2, velrange=None, bzwidth=None, plot=True):
    '''Calculate the Bz of an LSD profile
    
    :param lsd: lsd object (input). It is assumed that the lsd.vel is in km/s.
    :param cog: The value, or calculation method for the center of gravity. The choices are:
                'I': center of gravity of I,
                'V': center of gravity of V,
                'IV': center of gravity of I*V,
                a float: a user defined value in km/s.
    :param norm: calculation method for the continuum. The choices are:
                    'auto': the median of I outside of velrange (if defined) or the full range (if velrange is not defined)
                    float: a user defined value to use for Ic.
    :param lambda0: wavelength of the transition (default=500 nm or 5000 AA).
                    For an LSD profile, this is the lambda value the LSD profile shape was scaled with.
                    Needs to be a astropy unit object, with length units.
                    The astropy unit package will take care of the units conversion to give the Bz in Gauss.
    :param geff: effective Lande factor of the transition.
                 For an LSD profile, this is the geff value the LSD profile shape was scaled with.
    :param velrange: range of velocity to use for the determination of the
                     line center and the continuum. If not defined, the whole range
                     will be used. If bzwidth is not defined, this range will also be
                     used for the line Bz calculation.
    :param bzwidth: distance from the line center to use in the Bz calculation.
                    One element = same on each side of line center.
                    Two elements, left and right of line center.
                    Not defined: using velrange.
    :param plot: whether or not a graph is generated, and returned.
    :return: a dictionary with Bz and FAP calculations,
             optionally a matplotlib figure.
    '''

    # Velrange is used to identify the position of the line,
    # for calculating the cog, and for calculating the position
    # of the continuum.
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
            print('  no range in velocity given, using the median of the whole specI to determine continuum')
            norm_val = np.median(lsd_in.specI)
    else:
        norm_val = copy.copy(norm)
        print('using given norm value')

    # Check if cog is a string.
    if isinstance(cog, str):
        # calculate the cog for the chosen method
        if cog == 'I':
            cog_val = cog_I(lsd_in, norm_val)
        elif cog == 'min':
            cog_val = cog_min(lsd_in)
        elif cog == 'IV':
            cog_val = cog_IV(lsd_in, norm_val)
        elif cog == 'V':
            cog_val = cog_V(lsd_in)
        else:
            raise ValueError('calcBz got unrecognized value for cog: {:}'.format(cog))
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
                raise ValueError('bzwidth has too many elements {:} (need 1 or 2)'.format(len(bzwidth)))
        else:
            p_bzwidth = [cog_val-bzwidth, cog_val+bzwidth]
            lsd_bz = lsd[ np.logical_and(lsd.vel >= p_bzwidth[0], lsd.vel <= p_bzwidth[1]) ]

    # u.G is not properly defined in astropy.
    # It is defined in Tesla base units.
    # So here we are creating the right cgs base units for a Gauss.
    G_cgs = u.def_unit('G_cgs', 1 * u.g**0.5/u.cm**0.5/u.s)
    
    #Evaluate the equation for Bz; use V, Null1, and Null2 if they exist
    blv, blvSig = integrateBz(lsd_bz.vel, lsd_bz.specV, lsd_bz.specSigV, 
                geff, lambda0, cog_val, lsd_bz.specI, lsd_bz.specSigI, norm_val)
    if lsd.numParam > 2:
        bln1, bln1Sig = integrateBz(lsd_bz.vel, lsd_bz.specN1, lsd_bz.specSigN1,
                geff, lambda0, cog_val, lsd_bz.specI, lsd_bz.specSigI, norm_val)
    else: bln1, bln1Sig = (0*G_cgs, 0*G_cgs)
    if lsd.numParam > 3:
        bln2, bln2Sig = integrateBz(lsd_bz.vel, lsd_bz.specN2, lsd_bz.specSigN2,
                geff, lambda0, cog_val, lsd_bz.specI, lsd_bz.specSigI, norm_val)
    else: bln2, bln2Sig = (0*G_cgs, 0*G_cgs)
    
    # Get the FAP in the same range as the one used for Bz
    FAP_V, FAP_N1, FAP_N2 = FAP(lsd_bz)
    
    result = {
            'Ic': norm_val,
            'cog': cog_val,
            'Bzwidth min': p_bzwidth[0],
            'Bzwidth max': p_bzwidth[1],
            'V bz (G)': blv.value,
            'V bz sig (G)': blvSig.value,
            'V FAP': FAP_V,
            'N1 bz (G)': bln1.value,
            'N1 bz sig (G)': bln1Sig.value,
            'N1 FAP': FAP_N1,
            'N2 bz (G)': bln2.value,
            'N2 bz sig (G)': bln2Sig.value,
            'N2 FAP': FAP_N2
            }

    if plot:
        fig  = plotBzCalc(lsd, lsd_in, lsd_bz, velrange,
                          p_bzwidth, norm_val, cog_val, cog)
        return(result,fig)
    else:
        return(result)


def integrateBz(vel, spec, specSig, geff, lambda0, cog_val, specI, specSigI, norm_val):
    #Evaluate the integral equation for Bz.  Only intended for use in calcBz.
    
    # u.G is not properly defined in astropy.
    # It is defined in Tesla base units.
    # So here we are creating the right cgs base units for a Gauss.
    G_cgs = u.def_unit('G_cgs', 1 * u.g**0.5/u.cm**0.5/u.s)
    # This is the constant for the Zeeman splitting
    # Lambda_B = constant * lambda0**2 B
    lambda_B_constant = const.e.esu / (4 * np.pi * const.m_e.cgs * const.c.cgs**2)
    
    # set the velocity step for error propagation
    # and set to km/s (because the lsd profile objects don't have units associated with them)
    deltav = (vel[1] - vel[0])*u.km/u.s # This is in km/s
    
    # Calculation of the integral in the numerator of the Bz function with a trapezoidal numerical integral
    # For the error calculation, we propagate like we would for summation numerical integral.
    fnum = np.trapz( (vel - cog_val) * spec, x=vel-cog_val )*(u.km/u.s)**2 # This is in (km/s)^2
    sfnum = np.sqrt(np.sum( (vel - cog_val )**2 * specSig**2 )*(u.km/u.s)**2 * deltav**2)

    # Calculation of the integral in the denominator of the Bz function with a trapezoidal numerical integral
    # For the square error calculation, we propagate like we would for summation numerical integral.
    ri0v = np.trapz(norm_val-specI, x=vel )*u.km/u.s # This is in km/s
    si0v = np.sqrt(np.sum(specSigI**2 )* deltav**2) # This will naturally be in km/s

    # Make the actual Bz calculation.
    # for the units to work out, lambda0 needs to be passed
    # as a unit quantity (e.g. u.nm or u.AA)
    bl = (-1*fnum / ( ri0v*geff*lambda0*const.c*lambda_B_constant)).to(G_cgs)
    blSig = ( np.abs(bl * np.sqrt( (sfnum/fnum)**2 + (si0v/ri0v)**2 ))).to(G_cgs)
    return bl, blSig


def plotBzCalc(lsd, lsd_in, lsd_bz, velrange, p_bzwidth, norm_val, cog_val, cog):
    """
    Generate a plot showing the center of gravity and integration ranges used in the calculation of Bz from and LSD profile.
    Called by the main calcBz function.

    :param lsd: the full LSD profile to plot
    :param lsd_in: slice of an LSD profile with the range used for COG calculation
    :param lsd_bz: slice of an LSD profile with the range used for integration of Bz
    :param velrange: the velocity range used for COG calculation
    :param p_bzwidth: the velocity range used for integration of Bz
    :param norm_val: the continuum level used for normalization
    :param cog_val: the final COG used for Bz calculation
    :param cog: the input COG flag/value given by the user
    :return: a matplotlib figure object
    """
    #This function relies on the plot method of the LSD profile class
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
    ax[0].fill_between(lsd_bz.vel[red], lsd_bz.specV[red], step='mid', color='red')
    ax[0].fill_between(lsd_bz.vel[blue], lsd_bz.specV[blue], step='mid', color='blue')
    if(len(ax) > 2):
        ax[1].fill_between(lsd_bz.vel[red], lsd_bz.specN1[red], step='mid', color='red')
        ax[1].fill_between(lsd_bz.vel[blue], lsd_bz.specN1[blue], step='mid', color='blue')
    if(len(ax) > 3):
        ax[2].fill_between(lsd_bz.vel[red], lsd_bz.specN2[red], step='mid', color='red')
        ax[2].fill_between(lsd_bz.vel[blue], lsd_bz.specN2[blue], step='mid', color='blue')
    
    return fig


###############################################################
#For running calcBz as a terminal program
if __name__ == "__main__":
    
    #Take input file names and velocity ranges as command line arguments,
    #with some additional optional control parameters.
    import argparse
    import iolsd
    parser = argparse.ArgumentParser(description="""
    Calculate the longitudinal magnetic field Bz from LSD profiles.
    Prints Bz from Stokes V and null profiles.
    Also prints the detection 'false alarm probability' (FAP)
    (1 - detection probability) for the profile.""")
    parser.add_argument("fileList", nargs='*',
                        help='LSD profile file(s), to calculate Bz for.  '
                        +'Can be more than one file.')
    parser.add_argument("-v", "--velRange", nargs=2, type=float,
                        metavar=('VEL1', 'VEL2'), required=True, 
                        help='Starting and ending velocity for the range used '
                        +'to calculate line centre and integrate.')
    parser.add_argument("-g", "--Lande", type=float,
                        help='The effective Lande factor used to normalize '
                        +'(scale) the LSD profile.')
    parser.add_argument("-l", "--wavelength", type=float,
                        help='The wavelength, in nm, used to normalize '
                        +'(scale) the LSD profile.')
    parser.add_argument("-p", "--plotFit", action='store_true',
                        help='Optional, plot information about the line range '
                        +'and centre used.')
    args = parser.parse_args()
    #Process the command line parameters
    fileList = []
    for fileName in args.fileList:
        fileList += [fileName]
    velRange = args.velRange
    lande = args.Lande
    wl0 = args.wavelength
    plotFit = args.plotFit
    if fileList == []:
        print('Provide an LSD profile file!\n(try -h for more info)\n')

    results=[]
    #Run the Bz calculation on any files provided
    for fileName in fileList:
        lsd = iolsd.read_lsd(fileName)
        #If there isn't a Lande factor and wavelength provided,
        #see if there is a value in the LSD profile header
        if lande == None:
            indLande = lsd.header.find('lande=')
            if indLande >= 0:
                lande = float(lsd.header[indLande+6:].split()[0])
            else:
                print('A reference effective Lande factor is needed')
        if wl0 == None:
            indWl0 = lsd.header.find('wl=')
            if indWl0 >= 0:
                wl0 = float(lsd.header[indWl0+3:].split()[0])
            else:
                print('A reference effective Wavelength is needed\n')
        
        #The main Bz calculation
        res = calcBz(lsd, cog='I', lambda0=wl0*u.nm, geff=lande,
                     velrange=velRange, plot=plotFit)
        
        #If we want to plot this figure, first unpack the two values returned,
        #then we need to run the matplotlib show function.
        if plotFit:
            res, fig  = res
            import matplotlib.pyplot as plt
            plt.show()
        results += [res]

    #Print the results out
    nPar = lsd.numParam
    txtline = ('file                             Bz_V(G)  sigma_V  ratio FAP_V')
    if nPar >= 3: txtline += '      Bz_N1(G) sigma_N1  ratio FAP_N1'
    if nPar >= 4: txtline += '     Bz_N2(G) sigma_N2  ratio FAP_N2'
    print(txtline)
    for i, res in enumerate(results):
        txtline = '{:30s} {:+9.2f} {:8.2f} {:6.2f} {:9.3e}'.format(
                       fileList[i], res['V bz (G)'], res['V bz sig (G)'],
                       res['V bz (G)']/res['V bz sig (G)'], res['V FAP'])
        if nPar >= 3: txtline += ' {:+9.2f} {:8.2f} {:6.2f} {:9.3e}'.format(
                res['N1 bz (G)'], res['N1 bz sig (G)'],
                res['N1 bz (G)']/res['N1 bz sig (G)'], res['N1 FAP'])
        if nPar >= 4: txtline += ' {:+9.2f} {:8.2f} {:6.2f} {:9.3e}'.format(
                res['N2 bz (G)'], res['N2 bz sig (G)'],
                res['N2 bz (G)']/res['N2 bz sig (G)'], res['N2 FAP'])
        print(txtline)
