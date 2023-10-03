#!/usr/bin/python3
#
# Generate an LSD line mask from a VALD line list.

import numpy as np
from . import iolsd

def make_mask(lineListFile, maskFile, depthCutoff=0.0, 
              wlStart=0.0, wlEnd=1e10, landeStart=-1e10, landeEnd=1e10,
              elementsUsed=[], elementsExclude=[],
              atomsOnly=True, includeNoLande=False, defaultLande=1.0):
    """
    Generate an LSD line mask from a VALD line list and save it.
    The main interface for generating new line masks with Python scrips.
    
    This generates an LSD line mask from atomic lines in the list,
    and estimates Lande factors for lines where possible.  Lines that 
    no Lande factor can be estimated for are omitted by default.
    This assumes input line list wavelengths in A, and returns a mask in nm.
    
    :param lineListFile: The name of the file containing the line list
    :param maskFile: The name of the to write the mask to
    :param depthCutoff: Only include lines in the mask deeper than this value (defaults to 0, all lines included)
    :param wlStart: Optionally, only use lines with wavelengths above this (note: value in nm!)
    :param wlEnd: Optionally, only use lines with wavelengths below this (note: value in nm!)
    :param landeStart:  Optionally, only use lines with effective Lande factors above this
    :param landeEnd: Optionally, only use lines with effective Lande factors below this
    :param elementsUsed: Optionally, provide a list of elements to include in the line mask. This should be a list of strings with the element symbols, e.g. ['C', 'O', 'Si', 'Fe'] (an empty list or a list starting with 'all' will include all elements)
    :param elementsExclude: Optionally, provide a list of elements to exclude from the line mask. Also a list of strings with the element symbols.
    :param atomsOnly: Only include atomic lines (no molecular lines) and exclude H lines (defaults to True)
    :param includeNoLande: Include lines in the mask even if no Lande factor can be estimated for them (defaults to False)
    :param defaultLande: The value of the effective Lande factor to use if no value can be estimated (only used if includeNoLande = True)
    :rtype: The mask object, in case further processing is desired
    """
    lineList = iolsd.read_VALD(lineListFile)
    
    mask = convert_list_to_mask(lineList, depthCutoff=depthCutoff,
                                atomsOnly=atomsOnly,
                                includeNoLande=includeNoLande,
                                defaultLande=defaultLande)
    
    mask = filter_mask(mask, depthCutoff=depthCutoff,
                       wlStart=wlStart, wlEnd=wlEnd,
                       landeStart=landeStart, landeEnd=landeEnd,
                       elementsUsed=elementsUsed,
                       elementsExclude=elementsExclude)
    
    mask.save(maskFile)
    return mask

    
def convert_list_to_mask(lineList, depthCutoff=0.0, atomsOnly = True,
                         includeNoLande = False, defaultLande=1.0):
    """
    Convert a VALD line list to an LSD line mask.
    
    This estimates Lande factors when VALD doesn't have a known value
    This also can automatically reject any H lines and any molecular lines.
    Note, this assumes an input line list in A, and returns a mask in nm.
    
    :param lineList: The line list object to be converted to a mask
    :param depthCutoff: Only include lines deeper than this depth cutoff (defaults to 0, all lines included)
    :param atomsOnly: Flag to use only atomic lines (no molecular lines) and exclude H lines (defaults to True)
    :param includeNoLande: Flag to include lines in the mask even if no Lande factor can be estimated for them (defaults to False)
    :param defaultLande: The value of the effective Lande factor to use for lines where no value can be estimated (only used if includeNoLande = True)
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
    nLines = lineList.nLines
    mask = iolsd.Mask(np.zeros(nLines), np.zeros(nLines), np.zeros(nLines),
                      np.zeros(nLines), np.zeros(nLines),
                      np.zeros(nLines,dtype=int))

    skippedIon = []
    missingLandeIon = []
    missingLandeNum = 0
    j = 0
    #Loop through the line list making sure to use good data
    for i in range(lineList.nLines):
        element = lineList.ion[i].split()[0]
        ion = lineList.ion[i].split()[1]

        if (lineList.depth[i] >= depthCutoff) and (
            (atomsOnly == False) or (element in elements and element != 'H')):
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
                #Track lines not used due to missing Lande factors
                missingLandeNum += 1
                if not lineList.ion[i] in missingLandeIon:
                    missingLandeIon += [lineList.ion[i]]
        elif (lineList.depth[i] >= depthCutoff):
            #Track lines not used due unsupported species
            if not lineList.ion[i] in skippedIon:
                skippedIon += [lineList.ion[i]]

    if missingLandeNum > 0:
        print('missing Lande factors for {:} lines (skipped) from:'.format(
            missingLandeNum))
        print(missingLandeIon)
    if len(skippedIon) > 0:
        print('skipped all lines for species:')
        print(skippedIon)
    
    #Convert wavelengths from A to values in nm
    mask.wl = mask.wl*0.1
    #Remove blank entries in the line mask
    mask = mask.prune()
    
    return mask


def estimateLande(J, config, verbose=False):
    """
    Estimate a Lande factor for a level (if VALD doesn't have one)

    This supports levels in LS coupling, JJ coupling, and JK coupling.
    Levels in LK coupling are not supported.  This relies on VALD3's text 
    format for electron term symbols and configurations.
    
    :param J: The J quantum number for this level
    :param config: The text string from VALD with the electron configuration
     and term symbol for the level
    :param verbose: Report warnings for levels in unsupported coupling schemes if True
    :rtype: The Lande factor estimated for this level.  If the Lande factor can't be estimated this will be 99.0.
    """
    
    #Letters for L quantum numbers (note J is omitted in this notation)
    Lval = ('S','P','D','F','G','H','I','K','L','M','N','O','Q','R','T','U','V')
    lval = ('s','p','d','f','g','h','i','k','l','m','n','o','q','r','t','u','v')

    #In the rare case of a blank string for the electron configuration, abort
    if(len(config.split()) < 2): return 99.00
    
    coupling = config.split()[0]

    #For levels in LS coupling
    if coupling == 'LS':
        #Parse the text string
        L, S = get_LS_numbers(config)
        #and return 99 if that failed.
        if L == -1 or S == -1.0: return 99.00
        
        #The Lande factor in L-S coupling
        #(e.g. Degl'Innocenti & Landolfi 2004 eq. 3.8)
        #note, levels with J = 0 do not split at all (but break this equation)
        if J == 0.0:
            lande = 0.0
        else:
            lande = 1.0 + (J*(J + 1.) - L*(L + 1.) + S*(S + 1.)) \
                / (2.*J*(J + 1.))
        
    #For levels in JJ coupling
    #This covers 'J1-J2' 'J1-j' and (I think) 'j-j' cases
    #In this coupling scheme, individual electron l and s values couple
    #to make j  then the individual j couple to the total angular momentum J.
    #Also, subgroups of electrons can couple together to produce a combined
    #angular momentum J1 or J2, and those then combined to produce he total J.
    #
    #For this we want the J1, J2, L1 and L2 values, and possibly 
    #S1, S2 if they are given (if there are are subgroups of electrons).
    #This requires parsing an electron configuration as well as a term symbol.
    elif coupling == 'JJ':
        #Parse the text string
        L1, S1, J1, L2, S2, J2 = get_JJ_numbers(config)
        #and return 99 if that failed.
        if L1 < 0 or L2 < 0 or S1 < 0.0 or S2 < 0.0 or J1 < 0.0 or J2 < 0.0:
            return 99.0
        
        if J == 0.0:
            lande = 0.0
        else:
            #Landi Degl'Innocenti & Landolfi 2004 Polarization in Spectral Lines
            #Sect 3.1 (page 77, second equation, unnumbered after eq 3.8).
            #Also: W.C. Martin, R. Zalubas, and L. Hagan,
            #Atomic Energy Levels - The Rare-Earth Elements, 
            #Natl. Stand. Ref. Data Ser., Natl. Bur. Stand.
            #(U.S. Government Printing Office, Washington, 1978), No. 60.,
            #Sect 2.4
            gam1 = (J*(J+1.) + J1*(J1+1.) - J2*(J2+1.))/(2.*J*(J+1.))
            gam3 = (J*(J+1.) + J2*(J2+1.) - J1*(J1+1.))/(2.*J*(J+1.))
            if J1 == 0:
                gam2 = 0.
            else:
                gam2 = (J1*(J1+1.) + S1*(S1+1.) - L1*(L1+1.))/(2.*J1*(J1+1.))
            if J2 == 0: gam4 = 0.
            else:
                gam4 = (J2*(J2+1.) + S2*(S2+1.) - L2*(L2+1.))/(2.*J2*(J2+1.))
            lande = 1.0 + gam1*gam2 + gam3*gam4
        ##print('JJ lande {:6.3f}'.format(lande))
    
    #For levels in JK coupling
    elif coupling == 'JK':
        #Parse the text string
        L1, S1, J1, L2, S2, K = get_JK_numbers(config)
        #and return 99 if that failed.
        if L1 < 0 or L2 < 0 or S1 < 0.0 or S2 < 0.0 or J1 < 0.0 or K < 0.0:
            return 99.0
        
        #From W.C. Martin, R. Zalubas, and L. Hagan, 1978,
        #Atomic Energy Levels - The Rare-Earth Elements, 
        #Natl. Stand. Ref. Data Ser., Natl. Bur. Stand., Sect 2.4
        #(They cite: Wybourne, B. G., 1965, Spectroscopic properties of Rare
        #Earths, John Wiley & Sons, New York, N.Y., pg. 100),
        #which is for the more general J1-L2 case of JK coupling.
        #Landi Degl'Innocenti & Landolfi 2004 also provide a formula,
        #but only for the J1-l case of a subgroup coupling with one electron,
        #i.e. case when S2 = 1/2. We use the more general form here.
        if J == 0.0:
            lande = 0.0
        else:
            if J1 == 0.0:
                gam1 = 0.0
            else:
                gam1 = (J1*(J1+1.) - L1*(L1+1.) + S1*(S1+1.))/(2.*J1*(J1+1.))
            gam2 = (K*(K+1.) + J1*(J1+1.) - L2*(L2+1.))/((2.*J+1.)*(2.*K+1.))
            gam3 = (3.*J*(J+1.) - K*(K+1.) + S2*(S2+1.))/(2.*J*(J+1.))
            lande = 2.*gam1*gam2 + gam3
        ##print('JK lande {:6.3f}'.format(lande))
    
    else:
        #LK coupling is not supported, since I am not aware of a good, simple
        #approximation for calculating Lande factors in this coupling scheme.
        #Molecular lines are not supported here either.
        lande = 99.0
        if verbose:
            print('unsupported coupling scheme: ', config)
    return lande


def get_LS_numbers(config):
    """
    Extract the L and S quantum numbers from a VALD term symbol text string,
    for a level in LS coupling.
    
    :param config: the electron configuration and term symbol text string
    :rtype: The L and S quantum numbers (if a parsing error occurs, negative values are returned)
    """
    #Letters for L quantum numbers (note J is omitted in this notation)
    Lval = ('S','P','D','F','G','H','I','K','L','M','N','O','Q','R','T','U','V')
    
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
            #in case the multiplicity for S is 2 digits
            ii = i - 1
            if termSym[0:i].isnumeric(): ii = 0
            try:
                S = (float(termSym[ii:i]) - 1.)/2.
            except ValueError:
                print('Failed to parse S in LS coupling, from:\n', config)
                return -1, -1.
    if L < 0 or S < 0:
        print('Failed to parse L and S in LS coupling, from:\n', config)
        return -1, -1.
    return L, S


def get_JJ_numbers(config):
    """
    Extract the L1, S1, J1, L2, S2, and J2 quantum numbers for a level in JJ 
    coupling from VALD's electron configuration and term symbol text string.

    This should cover the J1-J2 case (where two subgroups combine to produce
    J1 and J2, which combine to produce the total J), the J1-j case (where
    one subgroup makes J1, which combines with an additional electron with j
    combines to produce J), and it should work for the j-j case (for two
    equivalent electrons with their own j, which combine to make J).
    
    :param config: the electron configuration and term symbol text string
    :rtype: The L1, S1, J1, L2, S2, and J2 quantum numbers (if a parsing error occurs, negative values are returned)
    """
    #Letters for L quantum numbers (note J is omitted in this notation)
    Lval = ('S','P','D','F','G','H','I','K','L','M','N','O','Q','R','T','U','V')
    lval = ('s','p','d','f','g','h','i','k','l','m','n','o','q','r','t','u','v')
    #The J1 and J2 values are in brackets () at the end as the term symbol
    termSym = config.split()[-1]
    econfig = config.split()[-2]
    if not ',' in termSym:
        print('Failed to parse J1, J2 in JJ coupling, from:\n', termSym)
        return -1, -1., -1., -1, -1., -1.
    sJ1, sJ2 = termSym.strip('*').strip('()').split(',')
    #Often J1 and J2 are fractions and need to be converted
    J1 = _parseFraction(sJ1, errMsg='Failed to parse J1 in JJ coupling:')
    J2 = _parseFraction(sJ2, errMsg='Failed to parse J2 in JJ coupling:')
    if np.isnan(J1) or np.isnan(J2): return -1, -1., -1., -1, -1., -1.
    #The L and possibly S values are given in the configuration.
    #If there is a subgroup of electrons there should be a capital letter
    #in (), giving the L value, and before that the number 2S+1.
    #If there is a single electron there should be a lower case letter
    #giving l (the last lower case letter in the string for l2) and s = 1/2.
    #The J_i values for a subgroup of electrons is given here in <>,
    #but we already have that from the term symbol so don't need to
    #extract it here.
    L1 = -1
    L2 = -1
    S1 = -1.
    S2 = -1.
    #The L2 (or l2) number should be the last letter in the string (upper or
    # lower case).  If the letter is a capital, try the proceeding number as
    #a multiplicity to get an S2 value, otherwise for a single electron
    #just use s2 = 1/2.
    for i2 in range(len(econfig)-1, 0, -1):
        #for a subgroup with L2
        if econfig[i2] in Lval:
            L2 = Lval.index(econfig[i2])
            ii = i2-1
            if econfig[i2-2:i2].isnumeric(): ii = i2-2
            try:
                S2 = (float(econfig[ii:i2])-1.)/2.
            except ValueError:
                print('Failed to parse S2 in JJ coupling, from:\n', config)
                return -1, -1., -1., -1, -1., -1.
            break
        #for a single e- with l2
        elif econfig[i2] in lval:
            L2 = lval.index(econfig[i2])
            S2 = 0.5
            break
    #Try to find L1, a capital (in brackets) before the L2 (or l2)
    for i1 in range(i2-1, 0, -1):
        if econfig[i1] in Lval:
            L1 = Lval.index(econfig[i1])
            ii = i1-1
            if econfig[i1-2:i1].isnumeric(): ii = i1-2
            try:
                S1 = (float(econfig[ii:i1])-1.)/2.
            except ValueError:
                print('Failed to parse S1 in JJ coupling, from:\n', config)
                return -1, -1., -1., -1, -1., -1.
            break
    #If there is no capital L1 for a subgroup, try to find l1 for a single e-
    #(the j-j case).
    if L1 == -1:
        for i1 in range(i2-1, 0, -1):
            if econfig[i1] in lval:
                L1 = lval.index(econfig[i1])
                S1 = 0.5
                break

    if L1 < 0 or L2 < 0 or S1 < 0.0 or S2 < 0.0:
        print('Failed to parse JJ coupling, from:\n', config)

    return L1, S1, J1, L2, S2, J2


def get_JK_numbers(config):
    """
    Extract the L1, S1, J1, L2, S2, and K quantum numbers for a level in JK
    coupling from VALD's electron configuration and term symbol text string.

    This is sometimes called J1-l or J1-L2 coupling.

    :param config: the electron configuration and term symbol text string
    :rtype: The L1, S1, J1, L2, S2, and K quantum numbers (if a parsing error occurs, negative values are returned)
    """
    #Letters for L quantum numbers (note J is omitted in this notation)
    Lval = ('S','P','D','F','G','H','I','K','L','M','N','O','Q','R','T','U','V')
    lval = ('s','p','d','f','g','h','i','k','l','m','n','o','q','r','t','u','v')
    
    termSym = config.split()[-1]
    econfig = config.split()[-2]

    K = -1.
    #The K value is in square brackets [] at the end as the term symbol
    if (not '[' in termSym) or (not ']' in termSym):
        print('Failed to parse K in JK coupling, from:\n', termSym)
        return -1, -1., -1., -1, -1., -1.
    sK = termSym.split('[')[1].split(']')[0]
    #K may be a fraction and need to be converted
    K = _parseFraction(sK, errMsg = 'Failed to parse K in JK coupling, from:')
    if np.isnan(K): return -1, -1., -1., -1, -1., -1.
    #The multiplicity (2*s+1) of the second electron or subgroup is usually
    #given just before the [], but since it could theoretically be omitted
    #we will infer it from the configuration below.
    
    L1 = -1
    L2 = -1
    S1 = -1.
    S2 = -1.
    J1 = -1.
    #The L2 (or l2) number should be the last letter in the string (upper or
    # lower case).  If the letter is a capital, try the proceeding number as
    #a multiplicity to get an S2 value, otherwise for a single electron
    #just use s2 = 1/2.
    for i2 in range(len(econfig)-1, 0, -1):
        #for a subgroup with L2
        if econfig[i2] in Lval:
            L2 = Lval.index(econfig[i2])
            ii = i2-1
            if econfig[i2-2:i2].isnumeric(): ii = i2-2
            try:
                S2 = (float(econfig[ii:i2])-1.)/2.
            except ValueError:
                print('Failed to parse S2 in JK coupling, from:\n', config)
                return -1, -1., -1., -1, -1., -1.
            break
        #for a single e- with l2
        elif econfig[i2] in lval:
            L2 = lval.index(econfig[i2])
            S2 = 0.5
            break
    #Try to find L1, a capital (in brackets) before the L2 (or l2)
    for i1 in range(i2-1, 0, -1):
        if econfig[i1] in Lval:
            L1 = Lval.index(econfig[i1])
            ii = i1-1
            if econfig[i1-2:i1].isnumeric(): ii = i1-2
            try:
                S1 = (float(econfig[ii:i1])-1.)/2.
            except ValueError:
                print('Failed to parse S1 in JK coupling, from:\n', config)
                return -1, -1., -1., -1, -1., -1.
            #There should be a value after L1 in <> giving the J1 value
            if ('<' in econfig[i1:i2]) and ('>' in econfig[i1:i2]):
                sJ1 = econfig[econfig.find('<')+1 : econfig.find('>')]
                #J1 may be a fraction and need to be converted
                J1 = _parseFraction(sJ1, errMsg =
                                    'Failed to parse J1 in JK coupling:')
                if np.isnan(J1): return -1, -1., -1., -1, -1., -1.
            break
    
    #If there is no capital L1 for a subgroup, try to find l1 for a single e-
    #(the case with two valence electrons, essentially a j-l case).
    if L1 == -1:
        for i1 in range(i2-1, 0, -1):
            if econfig[i1] in lval:
                L1 = lval.index(econfig[i1])
                S1 = 0.5
                #In this case the J1 value is in the term symbol, before the [K]
                #this should be a fraction since L1 is integer, S1 is 1/2,
                #and J1 = L1 +/- S1.
                sJ1 = termSym.split('[')[0]
                #J1 likely is a fraction
                J1 = _parseFraction(sJ1, errMsg =
                                    'Failed to parse J1 in JK coupling:')
                if np.isnan(J1): return -1, -1., -1., -1, -1., -1.
                break

    if L1 < 0 or L2 < 0 or S1 < 0.0 or S2 < 0.0 or J1 < 0.0 or K < 0.0:
        print('Failed to parse JK coupling, from:\n', config)
    
    return L1, S1, J1, L2, S2, K


def _parseFraction(string, errMsg='failed to parse as a fraction:'):
    """
    Parse a text string fraction and convert it to a floating point number.
    
    This only works for fractions with no whole part in front
    (i.e. a numerator/denominator, or only an integer).
    
    :param string: Text string to convert
    :param errMsg: Error message if this operation fails
    :rtype: The converted value, will be NaN if the conversion fails.
    """
    val = float('NaN')
    if '/' in string:
        try:
            val = float(string.split('/')[0])/float(string.split('/')[1])
        except ValueError:
            print(errMsg, string)
            return float('NaN')
    else:
        try:
            val = float(string)
        except ValueError:
            print(errMsg, string)
            return float('NaN')
    return val


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


def filter_mask(mask, depthCutoff = 0.0, wlStart = 0.0, wlEnd = 1e10,
               landeStart = -1e10, landeEnd = 1e10,
               elementsUsed = [], elementsExclude = []):
    """
    Remove lines from an LSD line mask, based on their depth, wavelength,
    effective Lande factor, or element.

    Usually, only one of elementsUsed and elementsExclude should be used.
    Both can be used, in that case only elements in elementsUsed are included,
    and then any elements in elementsExclude are removed from that.
    (Filtering by element currently only works for atomic lines,
    not molecular lines)

    :param mask: The line mask object to filter lines out of.
    :param depthCutoff: Only include lines in the mask deeper than this value (defaults to 0, all lines included)
    :param wlStart: Optionally, only use lines with wavelengths above this (note: value in nm!)
    :param wlEnd: Optionally, only use lines with wavelengths below this (note: value in nm!)
    :param landeStart:  Optionally, only use lines with effective Lande factors above this
    :param landeEnd: Optionally, only use lines with effective Lande factors below this
    :param elementsUsed: Optionally, provide a list of elements to include in the line mask. This should be a list of strings with the element symbols, e.g. ['C', 'O', 'Si', 'Fe'] (an empty list or a list starting with 'all' will include all elements)
    :param elementsExclude: Optionally, provide a list of elements to exclude from the line mask. Also a list of strings with the element symbols.
    :rtype: The line mask object, with the requested lines removed.
    """
    
    elements=('H' ,'He','Li','Be','B' ,'C' ,'N' ,'O' ,'F' ,'Ne','Na','Mg',
              'Al','Si','P' ,'S' ,'Cl','Ar','K' ,'Ca','Sc','Ti','V' ,'Cr',
              'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
              'Rb','Sr','Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
              'In','Sn','Sb','Te','I' ,'Xe','Cs','Ba','La','Ce','Pr','Nd',
              'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf',
              'Ta','W' ,'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po',
              'At','Rn','Fr','Ra','Ac','Th','Pa','U' )
    
    #Remove lines that are outside the requested range in wavelength, Lande factor, or depth.
    indRemove = (mask.depth < depthCutoff) | (mask.wl < wlStart) | (mask.wl > wlEnd) \
        | (mask.lande < landeStart) | (mask.lande > landeEnd)
    
    #If there is also a list of specific elements to use
    if type(elementsUsed) != type([]): print('ERROR: list of elements to include not of type list')
    if elementsUsed != [] and elementsUsed[0].lower() != 'all':
        indKeep = np.zeros_like(mask.iuse, dtype=bool)
        #Set a flag to only keep elements in the requested list
        for elUse in elementsUsed:
            try:
                elNumber = elements.index(elUse.strip())+1
            except ValueError:
                print('ERROR: Unrecognized element ', elUse, ' aborting!')
                import sys
                sys.exit()
            indKeep = indKeep | (mask.element.astype(int) == elNumber)
        indRemove = indRemove | np.logical_not(indKeep)
    
    #Or if there is also a list of specific elements to exclude from the mask
    if type(elementsExclude) != type([]): print('ERROR: list of elements to exclude not of type list')
    if elementsExclude != [] :
        #if elementsUsed != [] and elementsUsed[0].lower() != 'all':
        #    print('Warning: got both a list of elements to include and to exclude. ')
        indExclude = np.zeros_like(mask.iuse, dtype=bool)
        #Set a flag to remove elements in the exclude list
        for elCut in elementsExclude:
            try:
                elNumber = elements.index(elCut.strip())+1
            except ValueError:
                print('ERROR: Unrecognized element ', elCut, ' aborting!')
                import sys
                sys.exit()
            indExclude = indExclude | (mask.element.astype(int) == elNumber)
        indRemove = indRemove | indExclude
    
    #Apply the filtering
    mask.iuse[indRemove] = 0
    
    #Remove blank entries in the line mask
    mask = mask.prune()
    
    return mask



###############################################################
#For running makeMask as a terminal program
if __name__ == "__main__":
    #Take input and output file names as command line arguments,
    #with some additional optional control parameters.
    import argparse
    parser = argparse.ArgumentParser(description='Generate an LSD line mask from a VALD line list. When Lande factors are not in the line list they will be estimated if possible.')
    parser.add_argument("lineList", help="Line list from VALD. This file should be an 'extract stellar' line list in the 'long' format.")
    parser.add_argument("mask", nargs='?', default='mask.dat', help='Save the line mask to this file.')
    parser.add_argument("-d", "--depth", default='0', help='Optional number, only include lines deeper than this depth cutoff in the mask (defaults to 0, all lines included)')
    parser.add_argument("-w1", default='0', help='Optional number, starting wavelength for lines included in the mask (in nm)')
    parser.add_argument("-w2", default='1e10', help='Optional number, ending wavelength for lines included in the mask (in nm)')
    parser.add_argument("-g1", default='-1e10', help='Optional number, lower cutoff in effective Lande factor for lines included in the mask')
    parser.add_argument("-g2", default='1e10', help='Optional number, upper cutoff in effective Lande factor')
    parser.add_argument("-e", "--elements", default='', help="Optional, elements to include in the mask, as list symbols. For multiple elements enclose in quotes, e.g. 'C O Fe Ni'. By default all elements are used.")
    parser.add_argument("-ex", "--excludeEl", default='', help="Optional, elements to exclude from the mask, as a list symbols. For multiple elements enclose in quotes, e.g. 'C O Fe Ni'. By default all elements are used.")
    args = parser.parse_args()
    
    llistName = args.lineList
    maskName = args.mask
    depthCutoff = float(args.depth)
    wl1 = float(args.w1)
    wl2 = float(args.w2)
    lande1 = float(args.g1)
    lande2 = float(args.g2)
    elementsUsed = args.elements.split()
    elementsExclude = args.excludeEl.split()
    
    #Run the line mask generating code
    print('reading line list from', llistName)
    mask = make_mask(llistName, maskName, depthCutoff = depthCutoff,
                     wlStart = wl1, wlEnd = wl2,
                     landeStart = lande1, landeEnd = lande2,
                     elementsUsed = elementsUsed, elementsExclude = elementsExclude,
                     includeNoLande = False)
    print('mask saved to', maskName)
