"""
General purpose support tools for other modules.
These are mostly intended to be used inside other functions,
but they could occasionaly still be useful on their own.
"""
import copy
import warnings

# constants
# speed of light in m/s
c_ms = 299792458.  
# speed of light in km/s
c_kms = c_ms/1000.

# atomic symbols of elements up to atomic number 99
# (elements with atomic number above that are unsupported,
#  but also highly unstable, and unlikely to be observed)
_atomicSym=('H' ,'He','Li','Be','B' ,'C' ,'N' ,'O' ,'F' ,'Ne','Na','Mg',
            'Al','Si','P' ,'S' ,'Cl','Ar','K' ,'Ca','Sc','Ti','V' ,'Cr',
            'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
            'Rb','Sr','Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
            'In','Sn','Sb','Te','I' ,'Xe','Cs','Ba','La','Ce','Pr','Nd',
            'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf',
            'Ta','W' ,'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po',
            'At','Rn','Fr','Ra','Ac','Th','Pa','U' ,'Np','Pu','Am','Cm',
            'Bk','Cf','Es' )

def doppler_shift_kms(wl, velocity):
    '''
    Doppler shift a value, or array of values, according to an input
    radial velocity in km/s.
    
    :param wl: the wavelengths to doppler shift
    :param velocity: the radial velocity in km/s
    :return: the Doppler shifted wavelengths
    '''
    _wl = copy.deepcopy(wl) #work on a copy
    _wl = _wl + _wl*velocity/c_kms
    return _wl

def vacuum_to_air(wl):
    '''
    Convert wavelengths in vacuum into wavelengths in air.

    The wavelenths must be in units of Angstroms.

    This uses the formula from Donald Morton (2000, ApJ. Suppl., 130, 403)
    for the refraction index, which is also the IAU standard:
    n = 1 + 0.0000834254 + 0.02406147 / (130 - s^2) + 0.00015998 / (38.9 - s^2)
    where s = 10^4 / lambda_vac, and lambda_vac is in Angstroms.
    The conversion is then: lambda_air = lambda_vac / n. 
    This formula comes from Birch and Downs (1994, Metrologia, 31, 315)
    and applies to dry air at 1 atm pressure and 15 C with 0.045% CO2
    by volume. The corrections to Edlen (1953, J. Opt. Soc. Am., 43, 339)
    are less than 0.0001 A at 2000 A and less than 0.001 A at 30000 A.

    This is the vacuum to air conversion that VALD3 uses, see their website:
    http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion

    :param wl: an array of wavelengths to convert (in A)
    :return: the converted wavelengths (in A)
    '''
    swl = 1e4/wl
    refraction = (1.0 + 0.0000834254 + 0.02406147/(130.0 - swl**2)
                  + 0.00015998/(38.9 - swl**2))
    wlConv = wl/refraction
    return wlConv

def air_to_vacuum(wl):
    '''
    Convert wavelengths in air into wavelengths in vacuum

    The wavelenths must be in units of Angstroms.

    This is based formula from Donald Morton (2000, ApJ. Suppl., 130, 403)
    for the refraction index of air, and applies to dry air at 1 atm pressure
    and 15 C with 0.045% CO2 by volume. This uses the inversion of that formula
    derived by N. Piskunov and documented on the VALD3 website
    http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
    specifically:
    # n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s^2)
          + 0.0001599740894897 / (38.92568793293 - s^2), 
    # where s = 10^4 / lambda_air and the conversion is:
    lambda_vac = lambda_air * n.  
    
    :param wl: an array of wavelengths to convert (in A)
    :return: the converted wavelengths (in A)
    '''
    swl = 1e4/wl
    refraction = (1.0 + 0.00008336624212083
                  + 0.02408926869968/(130.1065924522 - swl**2)
                  + 0.0001599740894897/(38.92568793293 - swl**2))
    wlConv = wl*refraction
    return wlConv

def atomic_symbol_to_number(symbol):
    '''
    Convert an atomic symbol (letters) into an atomic number

    (e.g. convert: "Fe" into 26)
    
    :param symbol: a text string for the atomic symbol (e.g. "Fe")
    :return: an integer with the atomic number
    '''
    if not isinstance(symbol, str):
        raise ValueError('In atomic_symbol_to_number: input value must be a '
                         'string')
    _symbol = symbol.strip()
    num = _atomicSym.index(_symbol) + 1
    return

def atomic_number_to_symbol(number):
    '''
    Convert an atomic number into a (text) symbol

    (e.g. convert: 26 into "Fe")
    
    :param number: an integer of the atomic number
    :return: an text string with the atomic symbol (e.g. "Fe")
    '''
    _num = number
    if isinstance(number, float):
        _num = int(number)
    elif not isinstance(number, int):
        raise ValueError('In atomic_number_to_symbol: input atomic number '
                         'must be an int')
    if _num < 1 or _num > 99:
        raise ValueError('In atomic_number_to_symbol: input atomic number '
                         'outside supported range')
    symbol = _atomicSym[_num - 1]
    return copy.deepcopy(symbol)

def ion_symbol_to_number(symbol):
    '''
    Convert an ion symbol (atom and ionization state) into a number

    Uses text strings in the form "Fe 1" for neutral iron, and
    "Fe 2" for singly ionized iron. Generates ion numbers in the form
    atomic number + ionization state/100, e.g. "Fe 1" is 26.00,
    "Fe 2" is 26.01, "C 3" is 6.02.
    
    :param symbol: a text string for the ion symbol
    :return: a float with the ion number
    '''
    if not isinstance(symbol, str):
        raise ValueError('In ion_symbol_to_number: input value must be a '
                         'string')
    _symbol = symbol.strip()
    element = _symbol.split()[0]
    ion = _symbol.split()[1]
    try:
        num = float(_atomicSym.index(element) + 1) + (float(ion) - 1)*0.01
    except ValueError:
        num = 0.00
    return num

def ion_number_to_symbol(number):
    '''
    Convert an ion symbol (atom and ionization state) into a number

    Uses ion numbers in the form atomic number + ionization state/100
    Generates text strings in the form "Fe 1" for neutral iron, and
    "Fe 2" for singly ionized iron e.g. 26.00 is "Fe 1",
    26.01 is "Fe 2", 6.02 is "C 3".
    
    :param number: a float with the ion number
    :return: a string with the ion symbol
    '''
    _num = number
    if isinstance(number, int):
        _num = float(int)
    elif not isinstance(number, float):
        raise ValueError('In ion_number_to_symbol: input ion number '
                         'must be a float')
    if _num < 1.0 or _num >= 100.0:
        raise ValueError('In ion_number_to_symbol: input atomic number '
                         'outside supported range')
    
    try:
        strIon = _atomicSym[int(_num) - 1]
        strIon = strIon + ' {:d}'.format(round(100.*(_num % 1.0) + 1.0))
    except:
        warnings.warn('\nIn ion_number_to_symbol: Could not '
                      'identify the element from the ion string,\n'
                      'trying to just use the number '
                      '{:.2f}'.format(number),
                      stacklevel=2)
        strIon = '{:.2f}'.format(_num)
    return strIon
