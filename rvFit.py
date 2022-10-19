#!/usr/bin/python3
#
# Tools for measuring radial velocities from LSD profiles.

import numpy as np
from scipy.optimize import curve_fit
try:
    import specpolFlow.iolsd as iolsd
except ModuleNotFoundError:
    import iolsd

def fitGaussianRV(lsd, velrange=None, plotFit=False, fullOutput=False):
    """
    Fit a Gaussian function to an LSD profile, to determine a radial velocity.

    :param lsd:      The LSD profile for fitting to.
    :param velrange: Range of velocity to fit include in the fit
                     (a tuple or list with 2 elements).
                     If not given, the whole range will be used.
    :param plotFit: If True, Plot the fit to the LSD profile with matplotlib.
    :param fullOutput: If True, Return all the Gaussian best fit parameters,
                       rather than just the RV and RV error
    :rtype: Tuple of best fit RV and RV error. If fullOutput is True then also
            includes continuum & error, amplitude & error, and width & error.
    """

    if velrange == None:
        velrange = (lsd.vel[1], lsd.vel[-2])
    elif not(isinstance(velrange, list) or isinstance(velrange, tuple)):
        raise TypeError('velrange in fitGaussianRV should be a list or tuple'
                        +', with two elements')
    indVelUse = (lsd.vel >= velrange[0]) & (lsd.vel <= velrange[1])

    #Set the initial estimates for the fitting parameters
    indMinI = np.argmin(lsd.specI[indVelUse])
    initVel = lsd.vel[indVelUse][indMinI]
    initCont = 1.0
    initAmpl = 0.1
    initWidth = 10.
    param0 = [initCont, initAmpl, initVel, initWidth]

    #run the fitting routine
    fitVal, covarVal = curve_fit(gaussProf, lsd.vel[indVelUse],
                                 lsd.specI[indVelUse], p0=param0,
                                 sigma=lsd.specSigI[indVelUse],
                                 absolute_sigma=True)
    fitSigma = np.sqrt(np.diag(covarVal))
    
    fitCont = fitVal[0]
    fitAmpl = fitVal[1]
    fitVel = fitVal[2]
    fitWidth = fitVal[3]
    fitContErr = fitSigma[0]
    fitAmplErr = fitSigma[1]
    fitVelErr = fitSigma[2]
    fitWidthErr = fitSigma[3]
    #print('cont {:.6f} {:.6f} ampl {:.6f} {:.6f} vel {:.4f} {:.4f} width {:.6f} {:.6f}'.format(
    #    fitCont, fitContErr, fitAmpl, fitAmplErr, fitVel, fitVelErr,
    #    fitWidth, fitWidthErr))

    #Optionally plot the fit to the observed LSD profile
    if plotFit == True:
        import matplotlib.pyplot as plt
        plt.errorbar(lsd.vel, lsd.specI,  yerr=lsd.specSigI, c='grey')
        plt.errorbar(lsd.vel[indVelUse], lsd.specI[indVelUse],
                      yerr=lsd.specSigI[indVelUse], c='k')
        plt.plot(lsd.vel[indVelUse], gaussProf(
            lsd.vel[indVelUse], fitCont, fitAmpl, fitVel, fitWidth), 'b.')
        plotSynVels = np.linspace(velrange[0], velrange[1], 1000)
        plt.plot(plotSynVels, gaussProf(
            plotSynVels, fitCont, fitAmpl, fitVel, fitWidth), 'b')
        plt.xlabel('Velocity (km/s)')
        plt.ylabel('I')
        plt.show()

    #Optionally return all fit parameters
    if fullOutput == True:
        return fitVel, fitVelErr, fitCont, fitContErr, fitAmpl, fitAmplErr, \
            fitWidth, fitWidthErr
    return fitVel, fitVelErr


## Define line profile (simple Gaussian)
def gaussProf(x, cont, ampl, center, width):
    y = cont - ampl*np.exp(-(x-center)**2/(2.*width**2))
    return y

def chi2_for_gaussProf(x, yobs, ysig, cont, ampl, center, width):
    y = gaussProf(x, cont, ampl, center, width)
    chi2 = np.sum( ((y - yobs)/ysig)**2 )
    return chi2


###############################################################
#For running fitGaussianRV as a terminal program
if __name__ == "__main__":
    #Take input file names as command line arguments,
    #with some additional optional control parameters.
    import argparse
    parser = argparse.ArgumentParser(description='Fit a Gaussian to an LSD line profile, to determine radial velocity.  Prints the radial velocity in the same velocity units as the LSD profile.')
    parser.add_argument("fileList", nargs='*',
                        help='LSD profile files, to determine RVs for.  Can be more than one file.')
    parser.add_argument("-v", "--velRange", nargs=2, type=float, metavar=('VEL1', 'VEL2'),
                        help='Optional starting and ending velocity for the range used in the fit')
    parser.add_argument("-p", "--plotFit", action='store_true',
                        help='Optional plot the fit to the LSD profile')
    args = parser.parse_args()
    #Process the command line parameters
    fileList = []
    for fileName in args.fileList:
        fileList += [fileName]
    velRange = args.velRange
    plotFit = args.plotFit

    print('fileName                               RV    error')
    #Loop over the files provided, and fit a radial velocity for each
    for fileName in fileList:
        lsd = iolsd.read_lsd(fileName)
        
        #Run the radial velocity fitting routine
        fitVel, fitVelErr = fitGaussianRV(lsd, velRange, plotFit=plotFit)

        print('{:30} {:10.4f} {:8.4f}'.format(fileName, fitVel, fitVelErr))
