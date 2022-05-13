#!/usr/bin/python3
#
# Generate an LSD line mask from a VALD line list.

import numpy as np
import iolsd


def make_mask(lineListFile, maskFile, depthCutoff = 0.0):
    """
    Generate an LSD line mask from a VALD line list and save it.
    The main interface for generating new line masks in Python scrips.

    This generates an LSD line mask from valid atomic lines in the list,
    and estimates Lande factors for lines where possible.  Lines where
    no Lande factor can be estimated are omitted.
    
    :param lineListFile: The name of the file containing the line list
    :param maskFile: The name of the to write the mask to
    :param depthCutoff: Only include lines in the mask deeper than this value (defaults to 0, all lines included)
    :rtype: The mask object, in case further processing is desired
    """
    lineList = iolsd.read_VALD(lineListFile)
    mask = convert_list_to_mask(lineList, depthCutoff = depthCutoff)
    mask.save(maskFile)
    
    return mask

    
def convert_list_to_mask(lineList, depthCutoff = 0.0, atomsOnly = True,
                         includeNoLande = False, defaultLande = 1.0):
    """
    Convert a VALD line list to an LSD line mask.

    This estimates Lande factors when VALD doesn't have a known value
    This also can automatically reject any H lines and any molecular lines.
    
    :param lineList: The line list object to be converted to a mask
    :param depthCutoff: Only include lines deeper than this depth cutoff (defaults to 0, all lines included)
    :param atomsOnly: Flag to use only atomic lines (no molecular lines) and exclude H lines (defaults to True)
    :param includeNoLande: Flag to use include lines in the mask even if no Lande factor can be estimated for them (defaults to False)
    :param defaultLande: The value of the effective Lande factor to use for lines where no value can be estimated (if includeNoLande = True)
    :rtype: The line mask object created from this data
    """
    
    elements=('H' ,'He','Li','Be','B' ,'C' ,'N' ,'O' ,'F' ,'Ne','Na','Mg',
              'Al','Si','P' ,'S' ,'Cl','Ar','K' ,'Ca','Sc','Ti','V' ,'Cr',
              'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
              'Rb','Sr','Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
              'In','Sn','Sb','Te','I' ,'Xe','Cs','Ba','La','Ce','Pr','Nd',
              'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf',
              'Ta','W' ,'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po',
              'At','Rn','Fr','Ra','Ac','Th','Pa','U' )

    mask = iolsd.mask(nLines = lineList.nLines)

    skippedIon = []
    j = 0
    #Loop through the line list making sure to use good data
    for i in range(lineList.nLines):
        element = lineList.ion[i].split()[0]
        ion = lineList.ion[i].split()[1]

        if (lineList.depth[i] >= depthCutoff) and ((atomsOnly == False)
            or (element in elements and element != 'H')):
            #Estimate an effective Lande factor if necessary, VALD flags 
            #levels without a known Lande factor using the value 99.0
            landeEff = lineList.landeEff[i]
            if lineList.landeEff[i] >= 90.:
                landeLo = lineList.landeLo[i]
                if landeLo >= 90.:
                    landeLo = estimateLande(lineList.Jlo[i], lineList.configLo[i])
                landeUp = lineList.landeUp[i]
                if landeUp >= 90.:
                    landeUp = estimateLande(lineList.Jup[i], lineList.configUp[i])
                #If we have managed to estimate Lande factors
                if landeLo < 90. and landeUp < 90.:
                    landeEff = getEffectiveLande(landeLo, landeUp, 
                                             lineList.Jlo[i], lineList.Jup[i])

            #If we can't calculate a Lande factor but still want to use the line
            if landeEff >= 90. and includeNoLande == True:
                landeEff = defaultLande

            #For lines with a Lande factor
            if landeEff < 90.:
                #Save the data for good lines
                mask.wl[j]     = lineList.wl[i]
                mask.depth[j]  = lineList.depth[i]
                mask.excite[j] = lineList.Elo[i]
                mask.lande[j]  = landeEff
                mask.iuse[j]   = 1
                #The element code here is the atomic number plus 
                #the ionization state/100, so Fe II becomes 26.01
                try:
                    mask.element[j] = float(elements.index(element)+1) \
                        + (float(ion)-1)*0.01
                except ValueError:
                    mask.element[j] = 0.00
                
                j += 1
        else:
            if not element in skippedIon:
                skippedIon += [element]

    if len(skippedIon) > 0:
        print('skipped lines for species:')
        print(skippedIon)
            
    #Convert wavelengths from A to values in nm
    mask.wl = mask.wl*0.1
    #Remove blank entries in the line mask
    mask.prune()
    
    return mask


def estimateLande(J, config):
    """
    Estimate a Lande factor for a level (if VALD doesn't have one)
    
    :param J: the J quantum number for this level
    :param config: the text string from VALD with the electron configuration
     and term symbol for the level
    :rtype: the Lande factor estimated for this level
    """
    
    #Letters for L quantum numbers (note J is omitted)
    Lval = ('S','P','D','F','G','H','I','K','L','M','N','O','Q','R','T','U','V')
    lval = ('s','p','d','f','g','h','i','k','l','m','n','o','q','r','t','u','v')
    coupling = config.split()[0]

    #For levels in LS coupling
    if coupling == 'LS':
        #The last set of characters should be the term symbol in LS coupling.
        #We want the total orbital angular momentum quantum number L
        #and the total spin quantum number S.
        #This should be in the form (2S+1)L, where L is a capital letter
        #and (2S+1) is a number, possibly followed by * to indicate odd parity.
        #L value letters go like S=0, P=1, D=2, F=3, G=4, and so on
        L = -1
        S = -1.
        termSym = config.split()[-1]
        for i in range(len(termSym)):
            if termSym[i] in Lval:
                L = int(Lval.index(termSym[i]))
                try:
                    S = (float(termSym[i-1]) - 1.)/2.
                except ValueError:
                    print('Failed to parse S in LS coupling, from:\n', config)
                    return 99.0
        if L < 0 or S < 0:
            print('Failed to parse L and S in LS coupling, from:\n', config)
            return 99.0
        
        #The Lande factor in L-S coupling
        #(e.g. Degl'Innocenti & Landolfi 2004 eq. 3.8)
        #note, levels with J = 0 do not split at all (but break this equation)
        if J == 0.:
            lande = 0.0
        else:
            lande = 1.0 + (J*(J + 1.) - L*(L + 1.) + S*(S + 1.)) \
                / (2.*J*(J + 1.))
        
    #For levels in JJ coupling
    elif coupling == 'JJ':
        lande = 99.0
    
    #For levels in JK coupling
    elif coupling == 'JK':
        lande = 99.0
    
    else:
        #LK coupling is not supported, since I am not aware of a good, simple
        #approximation for calculating Lande factors in this coupling scheme.
        #Molecular lines are not supported here either.
        lande = 99.0
        print('unsupported coupling scheme: ', coupling)
    return lande


def getEffectiveLande(landeLo, landeUp, Jlo, Jup):
    """
    Calculate an effective Lande factor for a transition.
    
    This uses the upper and lower level Lande factors and J quantum numbers.
    The effective Lande factor is similar to an average weighted by the
    splitting pattern of the line.  It parameterizes the displacement of
    the centre of gravity of the sigma components of a Zeeman split line.
    
    :param landeLo: Lower level Lande factor
    :param landeUp: Upper level Lande factor
    :param Jlo: Lower level J quantum number
    :param Jup: Upper level J quantum number
    :rtype: Effective Lande factor for the transition
    """

    #From, e.g., Degl'Innocenti & Landolfi 2004, eq. 3.44 
    #(originally from Shenstone and Blair (1929) I think)
    landeEff = 0.5*(landeUp + landeLo) + 0.25*(landeUp - landeLo)*(Jup*(Jup + 1) - Jlo*(Jlo + 1))
    
    return landeEff


###############################################################
#For running this as a terminal program
if __name__ == "__main__":
    #Take input and output file names as command line arguments
    import argparse
    parser = argparse.ArgumentParser(description='Generate an LSD line mask from a VALD line list. When Lande factors are not in the line list they will be estimated if possible.')
    parser.add_argument("lineList", nargs='?', default='lines.dat', help="Line list from VALD. This should be an 'extract stellar' line list in the 'long' format.")
    parser.add_argument("mask", nargs='?', default='mask.dat', help='Save the line mask to this file.')
    parser.add_argument("-d", "--depth", default='0', help='Optional number, only include lines deeper than this depth cutoff in the mask (defaults to 0, all lines included)')
    args = parser.parse_args()

    llistName = args.lineList
    maskName = args.mask
    depthCutoff = float(args.depth)
    
    #Run the line mask generating code
    print('reading line list from', llistName)
    mask = make_mask(llistName, maskName, depthCutoff = depthCutoff)
    print('mask saved to', maskName)
