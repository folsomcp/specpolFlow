#!/usr/bin/python3
#
# Calculate a longitudinal magnetic field (Bz) from an LSD profile.

import argparse
import matplotlib.pyplot as plt
try:
    import specpolFlow as pol
except ModuleNotFoundError:
    #If specpolFlow is not installed or in the python path
    #try guessing the location and adding it temporarily
    import sys, os
    loc0 = os.path.dirname(__file__) 
    sys.path.insert(0, os.path.join(loc0, '..', '..'))
    import specpolFlow as pol

def bz_cli():
    """Main function for running LSD.calc_bz() as a terminal program

    Takes no arguments, but uses command line parameters instead.
    """
    #Take input file names and velocity ranges as command line arguments,
    #with some additional optional control parameters.
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
        lsd = pol.read_lsd(fileName)
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
        res = lsd.calc_bz(cog='I', lambda0=wl0, geff=lande,
                          velrange=velRange, plot=plotFit)
        
        #If we want to plot this figure, first unpack the two values returned,
        #then we need to run the matplotlib show function.
        if plotFit:
            res, fig  = res
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
        
    return

#For running this Python script as a terminal program
if __name__ == "__main__":
    bz_cli()
