#!/usr/bin/python3
#
# Fit sinusoids to data, to generate a periodogram.

import numpy as np
import matplotlib.pyplot as plt
#for confidence level calculations also include:
import scipy.special as specialf
import scipy.optimize as optimize

def periodogram_fourier(times, vals, errs, order=1,
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

    #If this will use more than ~2 GB of ram, be memory efficent but slower
    if pNum*times.size*(2*order+1) > 100000000:
        #This is the relatively slow version, using function calls in a loop
        print('using slower memory efficient mode')
        chi2 = np.zeros(pNum)
        for n, P in enumerate(periods):
            #build the design matrix M
            cycles = time0/P
            fitCoeff, fitErr, chi2_P, chi2nu_P = fourier_fit_data(
                times, vals, errs, P, jd0, order)
            chi2[n] = chi2nu_P
        dof = times.size-(2*order+1)

    else:
        #This version is faster (~10x) but uses more memory and is less clear.
        #Here we add another dimension to all the arrays, for the grid of 
        #periods. This allows us to avoid a time consuming loop, but makes the
        #array operations more confusing since we need to preserve the periods
        #dimension.
        #
        #array of rotation cycles for each observation at each test period
        cycles = np.outer(1./periods, time0)
        omegas = 2*np.pi*cycles
        #array of integers i for the i-th value in the Fourier series
        iorders = np.arange(1, order+1)
        
        #build the design matrix
        MM = np.zeros((pNum, times.size, 2*order+1))
        MM[:,:,0] = 1.0
        MM[:,:,1::2] = np.sin(np.einsum('ij,k->ijk', omegas, iorders))
        MM[:,:,2::2] = np.cos(np.einsum('ij,k->ijk', omegas, iorders))
        MM /= errs[np.newaxis, :, np.newaxis]
        
        #Evaluate the solution to the linear least squares problem
        valsOverSig = vals/errs
        auto = (MM.transpose((0,2,1)) @ valsOverSig)
        covar = MM.transpose((0,2,1)) @ MM
        covar = np.linalg.inv(covar)
        fit3 = np.einsum('ijk,ik->ij', covar, auto)
        
        #Evaluate the reduced chi^2 for the optimal solution at each period.
        chi = np.einsum('ijk,ik->ij', MM, fit3) - valsOverSig
        chi2 = np.einsum('ij,ij->i',chi, chi)
        dof = times.size-(2*order+1)
        chi2 /= dof
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
    fitCoeff, fitErr, chi2, chi2nu = fourier_fit_data(
        times, vals, errs, period, jd0, order)
    print('best fit reduced chi^2 {:.6f} for P={:.7} JD0={:}'.format(chi2nu,
                                                               period, jd0))
    
    print('  Fourier coefficents {:.6e} +/- {:.6e}'.format(
                                          fitCoeff[0], fitErr[0]))
    for i in range(1,len(fitCoeff)):
        print('                      {:.6e} +/- {:.6e}'.format(
                                          fitCoeff[i], fitErr[i]))
    
    #Generate the best fit curve for the plot
    phases = np.linspace(-0.5, 1.5, 200)
    yFourier = phase_curve_from_coeff(fitCoeff, phases)

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
    ax.errorbar(phaseO, vals2, yerr=errs2, fmt='o', c='k', zorder=1.1)
    ax.plot(phases, yFourier, zorder=2)
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
    
    return np.stack((phaseO, vals2, errs2)), np.stack((phases, yFourier)), fig, ax

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
    #fitCoeff = np.linalg.inv(MM.T @ MM) @ (MM.T @ valsOverSig)
    covar = np.linalg.inv(MM.T @ MM)
    auto = (MM.T @ valsOverSig)
    fitCoeff = covar @ auto
    fitErr = np.sqrt(np.diag(covar))

    #Get the reduced chi^2 for this optimal fit
    chi2 = (valsOverSig - MM@fitCoeff).T @ (valsOverSig - MM@fitCoeff)
    chi2nu = chi2/(cycles.size - fitCoeff.size)
    return fitCoeff, fitErr, chi2, chi2nu

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
    fitCoeff, fitErr, chi2, chi2nu = fourier_fit_data_phased(
        cycles, vals, errs, order)
    return fitCoeff, fitErr, chi2, chi2nu

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
    yFourier = M2 @ fitCoeff
    return yFourier

def time_curve_from_coeff(times, period, time0, fitCoeff):
    """
    Generate a sinusoidal curve, as a function of time, from a set of set of Fourier coefficients

    :param times: an array of times to evaluate the curve at, typically Julian Dates are used.
    :param period: the rotation period
    :param time0: the reference time (Julian Date) for rotation cycle 0.0.
    :param fitCoeff: an array of Fourier coefficients
                     (n=0, sin_n=1, cos_n=1, sin_n=2, cos_n=2, ...)
    :return: an array of the Fourier series evaluated at the input times
    """
    cycles = (times-time0)/period
    curve = phase_curve_from_coeff(fitCoeff, cycles)
    return curve

def confidence_from_periodogram(periods, chi2, numObs, order,
                                probLevels=[0.6827, 0.9545, 0.9973],
                                printDisagreement=False):
    """
    Estimate a range of periods for a given confidence level.

    This uses the variation in chi^2 about the minimum in chi^2 to find
    a range of periods that fall within a confidence level corresponding to
    that change in chi2.  This uses the order of the Fourier polynomial
    and the number of observed points to define the degrees of freedom for
    the problem.  When multiple local chi^2 minima fall within the chi^2
    contour, the period range spanning those minima will be reported.

    :param periods: an array of periods from the periodogram
    :param chi2: an array of reduced chi^2 values corresponding to the periods
    :param numObs: the number of observed datapoints
    :param order: the order of the Fourier polynomial used for calculating the periodogram
    :param probLevels: a list of probability/confidence levels to report
                       period ranges for (defaults to 1, 2, and 3 sigma)
    :param printDisagreement: flag to print (and return) the probability that
                              the model is inconsistent with the observations.
    :return: A dictionary containing the probability level as a key and the
             corresponding period range (in a list).  If printDisagreement
             is True there is also a 'P of  disagreement' entry in
             the dictionary.
    """
    dofChi2 = numObs - (2*order + 1)
    #for joint confidence on all fitting parameters including period
    dofProb = 2*order + 1 + 1 
    chi2Full = chi2*dofChi2
    indChi2Min = np.argmin(chi2)
    chi2Min = chi2Full[indChi2Min]
    periodBest = periods[indChi2Min]
    periodRanges = {}

    if printDisagreement:
        # Probability that the model is inconsistent with the data
        probOfMin = specialf.gammainc(dofChi2/2., chi2Min/2.)
        print('prob. of dissagreement for chi^2 min: '
              +'{:} (at chi2 {:} reduced chi2 {:}'.format(
                  probOfMin, chi2Min, chi2Min/dofChi2))
        periodRanges['P of  disagreement'] = probOfMin
    
    # Probability of a change in chi^2 (deltaChi2) from chi^2 minimum:
    #     P = gammainc(nParam/2, deltaChi2/2)
    # where nParam is the number of parameter's whose joint confidence
    # region you want, and deltaChi2 is the change in chi^2 from the minimum.
    # gammainc is the regularized lower incomplete Gamma function: 
    #     gammainc = 1 / gamma(a) * integral(exp(-t) * t**(a-1), t=0..x)
    # To find the change in chi^2 corresponding to confidence level P,
    # we need the root of:
    #     gammainc(nParam/2, deltaChi2/2) - P
    def _rootForDeltaChi2(deltaChi2, nParam, targetP):
        return specialf.gammainc(nParam/2., deltaChi2/2.) - targetP

    #range limits for root finding of chi^2 probability function
    rangeStart = 1e-3
    rangeEnd = 1e9

    for targetP in probLevels:
        #Find the root of the function (brentq seems to be relatively efficent)
        #to get the target change in chi^2
        rootDeltaChi = optimize.brentq(_rootForDeltaChi2, rangeStart, rangeEnd,
                                       args=(dofProb, targetP) )
        print('P', targetP, 'deltaChi2', rootDeltaChi, 'dof', dofProb)
        targetChi2 = chi2Min + rootDeltaChi
        inRange = chi2Full <= targetChi2
        periodsInRange = periods[inRange]
        pmin = np.min(periodsInRange)
        pmax = np.max(periodsInRange)
        periodRanges['{:6.4f}'.format(targetP)] = [pmin, pmax]

        print('prob {:6.4f} : period {:.6f} +{:.6f} / -{:.6f} '.format(
              targetP, periodBest, pmax-periodBest, periodBest-pmin)
              +'(chi2 level {:.5f}, reduced chi2 {:.6f})'.format(
              targetChi2, targetChi2/dofChi2))

    return periodRanges

def estimate_dipole_from_Bz_curve(inclination, Bmin, Bmax, limbDark):
    """
    Estimate a dipole magnetic field strength and obliquity

    This uses a simplistic oblique rotator model for the magnetic field
    (Stibs 1950; Preston 1967) and estimates the strength at the magnetic pole
    and the obliquity angle (beta, the angle between the rotation axis and
    dipole axis). This assumes a purely dipolar field. In cases where
    the real field departs from a pure dipole this will likely overestimate
    the dipole field strength.  Also, if either the inclination or obliquiy
    approach 0 or 90, the equations may behave poorly.

    :param inclination: the inclination of the stellar rotation axis in degrees
    :param Bmin: the minimum value of the longitudinal magnetic field curve
    :param Bmax: the maximum value of the longitudinal magnetic field curve
    :param limbDark: a limb darkening coefficient
    :return: The dipole field strength, and the dipole obliquity angle
    """
    # In theory for a dipole:
    # Bz(p) = Bp*(15 + u)/(20*(3 - u))*(cos(beta)*cos(i)
    #                                   + sin(beta)*sin(i)*cos(2pi*p))
    # for rotation phase p, limb darkening coefficient u,
    # dipole field strength Bp and obliquity beta.  
    # With the ratio of Bz_min / Bz_max = r:
    # r = (cos(beta)*cos(i) - sin(beta)*sin(i))
    #    /(cos(beta)*cos(i) + sin(beta)*sin(i))
    #   = cos(i + beta) / cos(i - beta)
    # tan(beta) = (1 - r)/(1 + r)*cot(i)
    incRad = inclination*np.pi/180.
    rB = Bmin/Bmax
    betaRad = np.arctan((1. - rB)/(1. + rB)/np.tan(incRad))
    beta = betaRad*180./np.pi

    limbCoeff = (15. + limbDark)/(20.*(3. - limbDark))
    #using the above Bz equation (for the phase of max Bz) and a trig identity
    Bp = Bmax/(limbCoeff*np.cos(incRad - betaRad))

    #Keep results in the desired domain
    #(here beta points to the +ve pole: Bp >= 0; 0 <= beta <= 180)
    Bp = abs(Bp)
    if beta < 0.0:
        beta += 180
    return Bp, beta


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
    periods, chi2 = periodogram_fourier(jds, vals, errs, order,
                                        pStart, pEnd, pNum, jd0=jd0,
                                        save=True, plot=bPlot)

    #Get a best period if we need one
    if pFixed == None:
        pFixed = periods[np.argmin(chi2)]

    confidence_from_periodogram(periods, chi2, jds.size, order,
                                printDisagreement=True)
    
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

    ### extras ##
    #
    #if jd0 == None: jd0 = np.min(jds)
    #fitCoeff, fitErr, tmpchi2, tmpchi2nu = fourier_fit_data(jds, vals, errs,
    #                                                pFixed, jd0, order)
    #times = np.linspace(2456000., 2461000.0, 5000)
    #timeCurve = time_curve_from_coeff(times, pFixed, jd0, fitCoeff)
    ##save the time curve
    #with open('pTimeFitCurve.dat', 'w') as foutTC:
    #    for i in range(timeCurve.shape[0]):
    #        foutTC.write('{:9.6f} {:13.6e}\n'.format(times[i],
    #                                            timeCurve[i]))
    #print('time for phase 0:', times[np.argmax(np.abs(timeCurve))],
    #      'at Bl', timeCurve[np.argmax(np.abs(timeCurve))])
    #
    #Bp, beta = estimate_dipole_from_Bz_curve(60., np.min(phasedCurve),
    #                                         np.max(phasedCurve), 0.5)
    #print(Bp, beta)
    #
    #incs = np.linspace(1.0, 89., 88)
    #Bps = np.zeros(88)
    #betas = np.zeros(88)
    #Bps2 = np.zeros(88)
    #betas2 = np.zeros(88)
    #Bps3 = np.zeros(88)
    #betas3 = np.zeros(88)
    #for i in range(incs.size):
    #    #Bps[i], betas[i] = estimate_dipole_from_Bz_curve(incs[i], -39.0, 48.0, 
    #    #                                                 0.5)
    #    #Bps2[i], betas2[i] = estimate_dipole_from_Bz_curve(incs[i], -29.0, 48.0, 
    #    #                                                   0.5)
    #    #Bps3[i], betas3[i] = estimate_dipole_from_Bz_curve(incs[i], -49.0, 48.0, 
    #    #                                                   0.5)
    #    ##Bps2[i], betas2[i] = estimate_dipole_from_Bz_curve(incs[i], -39.0, 43.0, 
    #    ##                          1.0)
    #    ##Bps3[i], betas3[i] = estimate_dipole_from_Bz_curve(incs[i], -39.0, 53.0, 
    #    ##                          1.0)
    #    Bps[i], betas[i] = estimate_dipole_from_Bz_curve(incs[i],
    #                                np.min(phasedCurve), np.max(phasedCurve), 0.5)
    #
    #plt.show()
    #plt.plot(incs, Bps)
    #plt.plot(incs, Bps2)
    #plt.plot(incs, Bps3)
    #plt.show()
    #plt.plot(incs, betas)
    #plt.plot(incs, betas2)
    #plt.plot(incs, betas3)
    #plt.show()
