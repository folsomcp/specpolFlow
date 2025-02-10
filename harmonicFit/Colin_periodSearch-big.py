#!/usr/bin/python3
#
# Fit sinusoids to data, to generate a periodogram.

import numpy as np
import matplotlib.pyplot as plt

def periodogram(times, vals, errs, order=1,
                pStart=0.1, pEnd=100., pNum=10000, jd0=None,
                save=False, plot=False):
    """
    Calculate a periodogram by fitting sinusoids to data.

    This generates a periodogram by fitting sinusoids and their harmonics
    (from calculated from a Fourier series) to observed data, for a
    range of periods.  The result a periodogram in reduced chi^2 and period.

    :param times: an array of times of the observation, typically
                  barycentric Julian Dates
    :param vals: an array of the observed values (e.g. Bz measurements)
    :param errs: an array of the error bars for the observed values
    :param order: the order of the Fourier series used to fit the data
                  (1 = fundamental only, 2 = fundamental + first harmonic)
    :param pStart: the starting period of the grid searched
    :param pEnd: the ending period of the grid searched
    :param pNum: the number of periods in the grid (logarithmically spaced)
    :param jd0: the reference Julian Date (time) for rotation cycle 0.0.
                If this is None the first observation is used.
    :param save: a Boolean flag to save the periodogram to a file (periods.dat)
    :param plot: a Boolean flag to display a plot of the periodogram

    :return: an array containing the periods used, and an array containing
             the reduced chi^2 for those periods.
    """
    #Check data for consistency
    if times.size != vals.size or times.size != errs.size:
        print('ERROR: difference in sizes if input data arrays: '
              +'times {:} vals {:} errs {:}'.format(
                  times.size, vals.size, errs.size))
        return np.zeros(0), np.zeros(0)
    if times.size <= 2*order+1:
        print('ERROR: too few observations for the number of free parameters'
              +' {:}, {:}\n'.format(times.size, 2*order+1))
    elif times.size <= 2*order+3:
        print('Warning: few observations for to the number of free parameters'
              +' {:}, {:}, this may be unstable'.format(times.size, 2*order+1))
    #Adopt a default reference Julian Date, if one wasn't explicitly given
    if jd0 == None: jd0 = np.min(times)
    time0 = times-jd0

    #Setup a grid of periods to search
    periods = np.logspace(np.log10(pStart), np.log10(pEnd), pNum)

    #We want to fit the data using sinusoid/Fourier series with fixed period,
    #and run through a grid of periods.  The function is then
    # y = a0 + sum_n( a_n sin(n*2*pi*T/P) + b_n cos(n*2*pi*T/P) )
    #This can be turned into a linear problem like M.x = y (. == dot product)
    #by using x = [a0, a1, b1, a2, b2, ... a_npar, b_npar]
    #Then M has dimension (num observed times, free Fourier parameters) and
    # M[T,:] has rows [1, sin(1*2*pi*T/P), cos(1*2*pi*T/P), 
    #                     sin(2*2*pi*T/P), cos(2*2*pi*T/P),
    #                     sin(3*2*pi*T/P), cos(3*2*pi*T/P) ...,
    #                     sin(npar*2*pi*T/P), cos(npar*2*pi*T/P]

    ##This is the relatively clear but slow version,
    ##which loops over periods to fill in the chi^2 array for the periodogram
    #chi2 = np.zeros(pNum)
    #for n, P in enumerate(periods):
    #    #build the design matrix M
    #    omega = 2*np.pi*time0/P
    #    M1 = np.zeros((times.size, 2*order+1))
    #    M1[:,0] = 1.0
    #    for i in range(1, order+1):
    #        M1[:, 2*i-1] = np.sin(omega*i)
    #        M1[:, 2*i] = np.cos(omega*i)
    #    ##Using scipy.optimize.lsq_linear
    #    ##This find the x that minimizes: 0.5 * ||A x - b||**2
    #    ##We want to minimize chi2 = ((obs-y)/sigma)**2)
    #    ##                         = ||x*M/sigma - obs/sigma||**2
    #    ##for lsq_linear so use A = M/sigma, b = obs/sigma
    #    #fitRes = lsq_linear(M1/errs, vals/errs)
    #    #chi2[n] = np.sum((np.dot(M1,fitRes.x) - vals/errs)**2)
    #    #
    #    #Calculating the linear least squares solution in a few lines here 
    #    #seems to be 2x faster than scipy.optimize.lsq_linear
    #    #For linear function y we have a matrix M such that y = M.x
    #    #In a linear chi^2 fit, we want the x that minimizes chi^2:
    #    # chi^2 = (O - M.x)^T.S^2.(O - M.x)
    #    #for an observation O and diagonal matrix of observed errors S,
    #    #here ^T is a matrix transpose. 
    #    #Taking the derivative of chi^2 wrt x, and setting it to zero,
    #    #the minimum of chi^2 should be at:
    #    # x = (M^T.S^2.M)^(-1) . M^T.S^2.O   here ^(-1) is matrix inversion.
    #    #One can also distribute S^2 over M and O, essentially multiplying
    #    #M and O by S to get Ms and Os
    #    # a = (Ms^T.Ms)^(-1) . Ms^T.Os    which may be a bit faster.
    #    S2 = np.diag(errs**2)
    #    fit2 = np.linalg.inv(M1.T @ S2 @ M1) @ (M1.T @ S2 @ vals)
    #    chi2[n] =  (vals - M1@fit2).T @ S2 @ (vals - M1@fit2)

    #This is the relatively slow version, using function calls in a loop
    chi2 = np.zeros(pNum)
    for n, P in enumerate(periods):
        #build the design matrix M
        cycles = time0/P
        fitCoeff, chi2_P, chi2nu_P = fourier_fit_data(times, vals, errs, P, jd0, order)
        chi2[n] = chi2nu_P
    
    ##This version is faster (~10x) but uses more memory and is less clear.
    ##Here we add another dimension to all the arrays, for the grid of 
    ##periods. This allows us to avoid a time consuming loop, but makes the
    ##array operations more confusing since we need to preserve the periods
    ##dimension.
    ##
    ##Make an array of rotation cycles for each observation at each test period
    #cycles = np.outer(1./periods, time0)
    #omegas = 2*np.pi*cycles
    ##an array of integers i for the i-th value in the Fourier series
    #iorders = np.arange(1, order+1)
    #
    ##build the design matrix
    #MM = np.zeros((pNum, times.size, 2*order+1))
    #MM[:,:,0] = 1.0
    #MM[:,:,1::2] = np.sin(np.einsum('ij,k->ijk', omegas, iorders))
    #MM[:,:,2::2] = np.cos(np.einsum('ij,k->ijk', omegas, iorders))
    #MM /= errs[np.newaxis, :, np.newaxis]
    #
    ##Evaluate the solution to the linear least squares problem
    #valsOverSig = vals/errs
    #auto = (MM.transpose((0,2,1)) @ valsOverSig)
    #covar = MM.transpose((0,2,1)) @ MM
    #covar = np.linalg.inv(covar)
    #fit3 = np.einsum('ijk,ik->ij', covar, auto)
    #
    ##Evaluate the reduced chi^2 for the optimal solution at each period.
    #chi = np.einsum('ijk,ik->ij', MM, fit3) - valsOverSig
    #chi2 = np.einsum('ij,ij->i', chi, chi)
    dof = times.size-(2*order+1)
    #chi2 /= dof
    print('num obs points {:}, num fit params {:}, degrees of freedom {:}'.format(
        times.size, (2*order+1), dof))

    #Optionally save or plot the periodogram
    if save == True:
        fout = open('periods.dat', 'w')
        for i in range(pNum):
            fout.write('{:14.7e} {:13.6e}\n'.format(periods[i], chi2[i]))
        fout.close()

    if plot == True:
        plt.plot(periods, chi2, c='k', linewidth=1.0)
        plt.xlabel('period')
        plt.ylabel(r'reduced $\chi^2$')
        plt.show()

    indMin = np.argmin(chi2)
    print('best P', periods[indMin], 'at reduced chi^2', chi2[indMin])
    return periods, chi2


def plot_phased(times, vals, errs, period, jd0=None, order=1):
    """
    Generate a plot of phased observations with a best fit sinusoid

    This phases the observed data with a period and reference Julian Date,
    then fits a sinusoid (Fourier series) to the data, and generates
    a plot of the phased observations and best fit curve.

    :param times: an array of times of the observation, typically
                  barycentric Julian Dates
    :param vals: an array of the observed values (e.g. Bz measurements)
    :param errs: an array of the error bars for the observed values
    :param period: the period to phase the data with
    :param jd0: the reference Julian Date (time) for rotation cycle 0.0.
                If this is None the first observation is used.
    :param order: the order of the Fourier series used to fit the data
                  (1 = fundamental only, 2 = fundamental + first harmonic)
    :return: A 2D array containing phases of observations, the observed 
             values, and their errors.  A 2D array containing phases and 
             the best fit sinusoid at those phases.  A matplotlib figure
             object. The associated axes object, with a plot of the phased 
             data and best fit curve.
    """
    #Check data for consistency
    if times.size != vals.size or times.size != errs.size:
        print('ERROR: difference in sizes if input data arrays: '
              +'times {:} vals {:} errs {:}'.format(
                  times.size, vals.size, errs.size))
        return 0, 0, np.zeros((3,0)), np.zeros((2,0))
    
    ##Adopt a default reference Julian Date, if one wasn't explicitly given
    if jd0 == None: jd0 = np.min(times)

    #Get the best fit sinusoid/Fourier series for this period
    fitCoeff, chi2, chi2nu = fourier_fit_data(times, vals, errs, period, jd0, order)
    print('best fit reduced chi^2 {:.6f} for P={:.7} JD0={:}'.format(chi2nu,
                                                               period, jd0))
    
    #Generate the best fit curve for the plot
    phases = np.linspace(-0.5, 1.5, 200)
    yFourrier = phase_curve_from_coeff(fitCoeff, phases)

    #It's nice to plot from phase -0.5 to 1.5 to illustrate periodicity
    #So generate a copy of the data for plotting
    phaseO = np.mod((times-jd0)/period, 1)
    #(note: mod(-1.1, 1) = 0.9, instead of -0.1)
    phaseO2 = np.zeros_like(phaseO)
    phaseO2[phaseO <= 0.5] = phaseO[phaseO <= 0.5] + 1.0
    phaseO2[phaseO > 0.5] = phaseO[phaseO > 0.5] - 1.0
    phaseO = np.concatenate((phaseO, phaseO2))
    vals2 = np.tile(vals,2)
    errs2 = np.tile(errs,2)
                         
    #Plot the phased data and best fit curve
    fig, ax = plt.subplots(1,1)
    #ax.errorbar(phaseO, vals2, yerr=errs2, fmt='o', c='k', zorder=1.1)
    ax.plot(phaseO, vals2, 'k.', alpha=0.1, zorder=1.1)
    ax.plot(phases, yFourrier, zorder=2)
    ax.set_xlabel('Phase')
    ax.set_xlim(-0.5,1.5)
    ax.set_title('P = {:.5f}   JD0 = {:.4f}'.format(period, jd0))

    #More elaborate ticks (include minor ticks, ticks on inside of frame)
    from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
    #ax.xaxis.set_major_locator(MultipleLocator(0.5)) #ticks spaced 0.5 in phase
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(axis='both', which='both', direction='in',
                   top=True, right=True)
    
    return np.stack((phaseO, vals2, errs2)), np.stack((phases, yFourrier)), fig, ax

def fourier_fit_data_phased(cycles, vals, errs, order):
    """
    Fit a Fourier series (sinusoidal curve) to a set of observations
    with known rotation phases.
    
    :param cycles: an array of rotation cycles or phases of the observations
    :param vals: an array of the observed values (e.g. Bz measurements)
    :param errs: an array of the error bars for the observed values
    :param order: the order of the Fourier series used to fit the data
    :return: the best fit Fourier coefficients, the chi^2, and
             the reduced chi^2 of the fit
    """
    omegas = 2*np.pi*cycles
    #an array of integers i for the i-th value in the Fourier series
    iorders = np.arange(1, order+1)
    #Build the design matrix
    MM = np.zeros((cycles.size, 2*order+1))
    MM[:,0] = 1.0
    MM[:,1::2] = np.sin(np.outer(omegas, iorders))
    MM[:,2::2] = np.cos(np.outer(omegas, iorders))
    MM /= errs[:, np.newaxis]
    #Evaluate the solution to the linear least squares problem
    valsOverSig = vals/errs
    fitCoeff = np.linalg.inv(MM.T @ MM) @ (MM.T @ valsOverSig)

    #Get the reduced chi^2 for this optimal fit
    chi2 = (valsOverSig - MM@fitCoeff).T @ (valsOverSig - MM@fitCoeff)
    chi2nu = chi2/(cycles.size - fitCoeff.size)
    return fitCoeff, chi2, chi2nu

def fourier_fit_data(times, vals, errs, period, jd0=None, order=1):
    """
    Fit a Fourier series (sinusoidal curve) with a fixed period
    to a set of observations
    
    :param times: an array of times of the observation, typically
                  barycentric Julian Dates
    :param vals: an array of the observed values (e.g. Bz measurements)
    :param errs: an array of the error bars for the observed values
    :param period: the period to phase the data with
    :param jd0: the reference Julian Date (time) for rotation cycle 0.0.
                If this is None the first observation is used.
    :param order: the order of the Fourier series used to fit the data
    :return: the best fit Fourier coefficients, the chi^2, and
             the reduced chi^2 of the fit
    """
    #Get the best fit sinusoid/Fourier series for this period
    #Adopt a default reference Julian Date, if one wasn't explicitly given
    if jd0 == None: jd0 = np.min(times)
    cycles = (times-jd0)/period
    fitCoeff, chi2, chi2nu = fourier_fit_data_phased(cycles, vals, errs, order)
    return fitCoeff, chi2, chi2nu


def phase_curve_from_coeff(fitCoeff, phases):
    """
    Generate a phased sinusoidal curve from a set of set of Fourier coefficients

    :param fitCoeff: an array of Fourier coefficients
                     (n=0, sin_n=1, cos_n=1, sin_n=2, cos_n=2, ...)
    :param phases: an array of rotation phases to evaluate the curve at
    :return: an array of the Fourier series evaluated at the input phases
    """
    order = (fitCoeff.size - 1)//2
    #an array of integers i for the i-th value in the Fourier series
    iorders = np.arange(1, order+1)
    #Build the design matrix
    M2 = np.zeros((phases.size, 2*order+1))
    M2[:,0] = 1.0
    M2[:,1::2] = np.sin(np.outer(2*np.pi*phases, iorders))
    M2[:,2::2] = np.cos(np.outer(2*np.pi*phases, iorders))
    yFourrier = M2 @ fitCoeff
    return yFourrier



###############################################################
#For running periodSearch as a terminal program
if __name__ == "__main__":
    #Take input file names as command line arguments,
    #with some additional optional control parameters.
    import argparse
    pStart = 0.5
    pEnd = 50.
    pNum = 20000
    
    parser = argparse.ArgumentParser(description='Generate a periodogram by fitting a set of sinusoids to observed data. ')
    parser.add_argument("dataFile",
                        help='File of observed data. This text file should have columns of Julian Date, observed value, error on that value.')
    parser.add_argument("-o", "--order", type=int, default=1,
                        help="Order of the sinusoid (Fourier series) used for "
                        +"the fit, default: {:}".format(1))
    parser.add_argument("-s", "--periodStart", type=float, default=pStart,
                        help="Starting period for the search range, default: "
                        +"{:}".format(pStart))
    parser.add_argument("-e", "--periodEnd", type=float, default=pEnd,
                        help="Ending period for the search range, default: "
                        +"{:}".format(pEnd))
    parser.add_argument("-n", "--numPeriods", type=int, default=pNum,
                        help="Number of pixels in the periodogram, default: "
                        +"{:}".format(pNum))
    parser.add_argument("-j", "--JD0", type=float,
                        help='Reference Julian Date for rotation cycle = 0.0, '
                        +'defaults to the earliest observation')
    parser.add_argument("-p", "--plotPeriodogram", action='store_true',
                        help="Plot the periodogram of reduced chi^2")
    parser.add_argument("-f", "--fixedPeriod", type=float,
                        help='Optional, instead of searching a grid of periods,'
                        +' just calculate the phase curve for this period')
    args = parser.parse_args()

    inFile = args.dataFile
    order = args.order
    pStart = args.periodStart
    pEnd = args.periodEnd
    pNum = args.numPeriods
    jd0 = args.JD0
    bPlot = args.plotPeriodogram
    pFixed = args.fixedPeriod

    #Read the observed data points
    jds, vals, errs = np.loadtxt(inFile, unpack=True, usecols=(0,1,2))

    #Print the read in data values as a confirmation
    print("read:                    dates"
          +"               values"
          +"               errors")
    for i in range(jds.size):
        print("{:30} {:20} {:20}".format(jds[i], vals[i], errs[i]))

    #Generate the periodogram
    periods, chi2 = periodogram(jds, vals, errs, order,
                                pStart, pEnd, pNum, jd0=jd0,
                                save=True, plot=bPlot)

    #Get a best period if we need one
    if pFixed == None:
        pFixed = periods[np.argmin(chi2)]
    
    #If we want make a plot of the data phased with a period & best fit curve
    phasedVals, phasedCurve, fig, ax = plot_phased(jds, vals, errs, pFixed, 
                                                   jd0=jd0, order=order)
    if bPlot: plt.show()

    #Save the phased data and best fit curve, for further plotting later
    foutPD = open('phasedData.dat', 'w')
    for i in range(phasedVals.shape[1]):
        foutPD.write('{:9.6f} {:10} {:10}\n'.format(phasedVals[0,i],
                                    phasedVals[1,i], phasedVals[2,i]))
    foutPD.close()
    foutPC = open('phasedFitCurve.dat', 'w')
    for i in range(phasedCurve.shape[1]):
        foutPC.write('{:9.6f} {:13.6e}\n'.format(phasedCurve[0,i],
                                              phasedCurve[1,i]))
    foutPC.close()
