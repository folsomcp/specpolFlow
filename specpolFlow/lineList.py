"""
Tools for manipulating lists of spectral line data.
Includes tools for reading and parsing VALD version 3 line lists 
(typically 'long format', 'extract stellar' requests).
"""

import numpy as np

###################################

class LineList:
    """
    Container for a set of spectral line data, usually from VALD.

    This usually contains: 
     
    * nLines - number of lines in the line list
    * ion - list of species identifiers (element or molecule and ionization)
    * wl - array of wavelengths
    * loggf - array of oscillator strengths (log gf)
    * Elo - array of excitation potentials for the lower level in
      the transition (in eV)
    * Jlo - array of J quantum numbers for the lower level
    * Eu - array of excitation potentials for the upper level in
      the transition (in eV)
    * Jup - array of J quantum numbers for the upper level
    * landeLo - array of Lande factors for the lower level
    * landeUp - array of Lande factors for the upper level
    * landeEff - array of effective Lande factors for the transition
    * rad - array of radiative damping coefficients
    * stark - array of quadratic Stark damping coefficients
    * waals - array of van der Waals damping coefficients
    * depth - depth at the centre of the spectral line, as estimated by VALD
    * configLo - list of strings with the electron configuration and
      term symbols for the lower level
    * configUp - list of strings with the electron configuration and
      term symbols for the upper level
    * refs - list of references for the sources of the line data (optional)
    """
    def __init__(self, ion, wl, loggf, Elo, Jlo, Eup, Jup, landeLo, landeUp,
                 landeEff, rad, stark, waals, depth, configLo, configUp, refs):
        self.ion      = ion
        self.wl       = wl
        self.loggf    = loggf
        self.Elo      = Elo
        self.Jlo      = Jlo
        self.Eup      = Eup
        self.Jup      = Jup
        self.landeLo  = landeLo
        self.landeUp  = landeUp
        self.landeEff = landeEff
        self.rad      = rad
        self.stark    = stark
        self.waals    = waals
        self.depth    = depth
        self.configLo = configLo
        self.configUp = configUp
        self.refs     = refs
        self.nLines   = self.wl.size

    def __getitem__(self, key):
        """
        Returns a LineList with only the values at the specified index(s).

        :param key: the index or slice being checked
        :rtype: LineList
        """
        ion      = self.ion[key]
        wl       = self.wl[key]
        loggf    = self.loggf[key]
        Elo      = self.Elo[key]
        Jlo      = self.Jlo[key]
        Eup      = self.Eup[key]
        Jup      = self.Jup[key]
        landeLo  = self.landeLo[key]
        landeUp  = self.landeUp[key]
        landeEff = self.landeEff[key]
        rad      = self.rad[key]
        stark    = self.stark[key]
        waals    = self.waals[key]
        depth    = self.depth[key]
        configLo = self.configLo[key]
        configUp = self.configUp[key]
        refs     = self.refs[key]
        lList =  LineList(ion, wl, loggf, Elo, Jlo, Eup, Jup, landeLo,
                          landeUp, landeEff, rad, stark, waals, depth,
                          configLo, configUp, refs)
        return lList

    def __setitem__(self, key, newval):
        """
        Sets all values of the LineList at the specified location
        equal to the input LineList values.

        :param key: the index or slice being overwritten
        :param newval: LineList used to replace the values given by key
        """
        if not(isinstance(newval, LineList)):
            raise TypeError()
        else:
            self.ion[key]      = newval.ion
            self.wl[key]       = newval.wl
            self.loggf[key]    = newval.loggf
            self.Elo[key]      = newval.Elo
            self.Jlo[key]      = newval.Jlo
            self.Eup[key]      = newval.Eup
            self.Jup[key]      = newval.Jup
            self.landeLo[key]  = newval.landeLo
            self.landeUp[key]  = newval.landeUp
            self.landeEff[key] = newval.landeEff
            self.rad[key]      = newval.rad
            self.stark[key]    = newval.stark
            self.waals[key]    = newval.waals
            self.depth[key]    = newval.depth
            self.configLo[key] = newval.configLo
            self.configUp[key] = newval.configUp
            self.refs[key]     = newval.refs
            self.nLines   = self.wl.size

    def __len__(self):
        return self.nLines

    def __str__(self):
        """
        Generate a nicely formatted string of line data for printing
        """
        strList = []
        outStr = '\n'
        lines = self
        #check if this LineList contains arrays or just single values
        if isinstance(self.ion, str) and isinstance(self.wl, float):
            lines = [self]
        for line in lines:
            fmt = ("'{:s}',{:17.4f},{:7.3f},{:8.4f},{:5.1f},{:8.4f},{:5.1f},"
                   "{:7.3f},{:7.3f},{:7.3f},{:6.3f},{:6.3f},{:8.3f},{:6.3f},")
            if len(line.ion) == 3: fmt = fmt[:7] + ' ' + fmt[7:]
            strList.append(fmt.format(line.ion, line.wl, line.loggf,
                                line.Elo, line.Jlo, line.Eup, line.Jup,
                                line.landeLo, line.landeUp, line.landeEff,
                                line.rad, line.stark, line.waals, line.depth))
            strList.append("'  {:>88}'".format(line.configLo))
            strList.append("'  {:>88}'".format(line.configUp))
            strList.append("'{:}'".format(line.refs))
        return outStr.join(strList)

    def insert(self, i, newval):
        """
        Insert lines from a LineList at a specific index.  Returns a new 
        LineList combining the two line lists (does not operate in place)

        :param i: index to insert at
        :param newval: LineList of new values to insert
        :rtype: LineList
        """
        if not(isinstance(newval, LineList)):
            raise TypeError('LineList.insert: can only insert a LineList')

        #Re-allocating these arrays, particularly arrays of strings, is slow!
        ion      = np.insert(self.ion     , i, newval.ion     )
        wl       = np.insert(self.wl      , i, newval.wl      )
        loggf    = np.insert(self.loggf   , i, newval.loggf   )
        Elo      = np.insert(self.Elo     , i, newval.Elo     )
        Jlo      = np.insert(self.Jlo     , i, newval.Jlo     )
        Eup      = np.insert(self.Eup     , i, newval.Eup     )
        Jup      = np.insert(self.Jup     , i, newval.Jup     )
        landeLo  = np.insert(self.landeLo , i, newval.landeLo )
        landeUp  = np.insert(self.landeUp , i, newval.landeUp )
        landeEff = np.insert(self.landeEff, i, newval.landeEff)
        rad      = np.insert(self.rad     , i, newval.rad     )
        stark    = np.insert(self.stark   , i, newval.stark   )
        waals    = np.insert(self.waals   , i, newval.waals   )
        depth    = np.insert(self.depth   , i, newval.depth   )
        configLo = np.insert(self.configLo, i, newval.configLo)
        configUp = np.insert(self.configUp, i, newval.configUp)
        refs     = np.insert(self.refs    , i, newval.refs    )
        lList =  LineList(ion, wl, loggf, Elo, Jlo, Eup, Jup, landeLo,
                          landeUp, landeEff, rad, stark, waals, depth,
                          configLo, configUp, refs)
        return lList

    def save(self, fname):
        """
        Write a line list to a text file.
        This outputs using the VALD version 3 'extract stellar' 'long' format.
        
        A few details (e.g. references at the end of the file)
        are omitted since they are not saved in the LineList class.
        
        :param fname: the file name to save the output to
        """
        
        fOut = open(fname, 'w')
        fOut.write(("{:11.5f},{:12.5f},{:6d},{:7d},{:4.1f}, Wavelength "
                    "region, lines selected, lines processed, Vmicro\n"
                    ).format(self.wl[0], self.wl[-1], self.nLines, 999999, 0.))
        fOut.write("                                                 "
                   "                    Lande factors       "
                   "Damping parameters   Central\n")
        fOut.write("Spec Ion       WL_air(A)  log gf* E_low(eV) J lo "
                   "E_up(eV)  J up  lower   upper    mean   Rad.  "
                   "Stark   Waals   depth\n")
        fOut.write(str(self)+'\n')
        fOut.close()
        return
    
def line_list_zeros(nLines):
    """
    Generate a line list of zeros and blank text.
    
    Used by read_VALD
    (It can be a bit faster to allocate all the array space at once.)
    
    :param nLines: the number of lines in the LineList of zeros
    :rtype: LineList
    """
    ion      = np.tile(np.array([''], dtype='U6'), nLines)
    wl       = np.zeros(nLines)
    loggf    = np.zeros(nLines)
    Elo      = np.zeros(nLines)
    Jlo      = np.zeros(nLines)
    Eup      = np.zeros(nLines)
    Jup      = np.zeros(nLines)
    landeLo  = np.zeros(nLines)
    landeUp  = np.zeros(nLines)
    landeEff = np.zeros(nLines)
    rad      = np.zeros(nLines)
    stark    = np.zeros(nLines)
    waals    = np.zeros(nLines)
    depth    = np.zeros(nLines)
    configLo = np.tile(np.array([''], dtype='U128'), nLines)
    configUp = np.tile(np.array([''], dtype='U128'), nLines)
    refs = np.tile(np.array(['_          unknown source'],dtype='U180'), nLines)
    lList = LineList(ion, wl, loggf, Elo, Jlo, Eup, Jup, landeLo,
                     landeUp, landeEff, rad, stark, waals, depth,
                     configLo, configUp, refs)
    return lList

def read_VALD(fname):
    """
    Read a list of spectral line data from VALD and return a LineList.

    This expects VALD version 3 line list, in an 'extract stellar'
    'long' format.

    :param fname: the file name for the VALD line list.
    :rtype: LineList
    """
    fVald = open(fname, 'r')
    i = 0
    j = 0
    nLines = 0
    for txtLine in fVald:
        if i == 0:
            nLines = int(txtLine.split(',')[2])
            llist = line_list_zeros(nLines)
        #There should be 3 lines of header,
        #then 4 file lines for each set of spectra line data.
        if i > 2 and i < nLines*4+3:
            ii = (i-3)%4
            if ii == 0:
                #The line data are: Spec Ion, Wl(A), log gf, E_lower(eV), 
                #J_lower, E_upper(eV), J_upper, Lande factor lower, 
                #Lande factor upper, effective Lande factor, Rad. damping, 
                #Stark damping, van der Waals damping, central depth
                vals = txtLine.split(',')
                llist.ion[j]      = vals[0].strip('\'')
                llist.wl[j]       = float(vals[1])
                llist.loggf[j]    = float(vals[2])
                llist.Elo[j]      = float(vals[3])
                llist.Jlo[j]      = float(vals[4])
                llist.Eup[j]      = float(vals[5])
                llist.Jup[j]      = float(vals[6])
                llist.landeLo[j]  = float(vals[7])
                llist.landeUp[j]  = float(vals[8])
                llist.landeEff[j] = float(vals[9])
                llist.rad[j]      = float(vals[10])
                llist.stark[j]    = float(vals[11])
                llist.waals[j]    = float(vals[12])
                llist.depth[j]    = float(vals[13])
            if ii == 1:
                llist.configLo[j] = txtLine.strip(' \n\'')
            if ii == 2:
                llist.configUp[j] = txtLine.strip(' \n\'')
            if ii == 3:
                llist.refs[j] = txtLine.strip(' \n\'')
                j += 1
        i += 1
    
    fVald.close()
    return llist

###################################

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
    :return: Effective Lande factor for the transition
    """

    #From, e.g., Degl'Innocenti & Landolfi 2004, eq. 3.44 
    #(originally from Shenstone and Blair (1929) I think)
    landeEff = 0.5*(landeUp + landeLo) + 0.25*(landeUp - landeLo)*(Jup*(Jup + 1) - Jlo*(Jlo + 1))
    
    return landeEff


def estimateLande(J, config, verbose=False):
    """
    Estimate a Lande factor for a level (if VALD doesn't have one).

    This supports levels in LS coupling, JJ coupling, and JK coupling.
    Levels in LK coupling are not supported.  This relies on VALD3's text 
    format for electron term symbols and configurations.
    
    :param J: The J quantum number for this level
    :param config: The text string from VALD with the electron configuration
                   and term symbol for the level
    :param verbose: Report warnings for levels in unsupported coupling schemes
                    if True
    :return: The Lande factor estimated for this level.  If the Lande factor
             can't be estimated this will be 99.0.
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
    :return: The L and S quantum numbers (if a parsing error occurs,
             negative values are returned)
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
    to produce J), and it should work for the j-j case (for two
    equivalent electrons with their own j, which combine to make J).
    
    :param config: the electron configuration and term symbol text string
    :return: The L1, S1, J1, L2, S2, and J2 quantum numbers
             (if a parsing error occurs, negative values are returned)
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
    :return: The L1, S1, J1, L2, S2, and K quantum numbers
             (if a parsing error occurs, negative values are returned)
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
    :return: The converted value, will be NaN if the conversion fails.
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
