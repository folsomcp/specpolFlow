Some convention:

filename : fname
reading function read_objectName()
use save in classes for the functions that save


[ ] Rename ioLSD to something better, and arrange the init for direct access to these functions?


Class: lsd_obj
    * _init_
    * def save(self, fname)
    * def __len__(self)
    * def __getitem__(self, key):
    * def __setitem__(self, key, newval):
    
    * def norm(self, normValue):
    * def vshift(self, velShift):  ??? name consistency in parameters?

    * def set_weights(self, wint_old, wpol_old, wint_new, wpol_new):
    * def scale(self, scale_int, scale_pol):

    * def plot(self, figsize=(10,10), sameYRange=True, plotZeroLevel=True, **kwargs):
        [ ] add fig=None, ax=None so that another call can overplot another LSD profile. 
        [ ] maybe add an option for scatter with errors and just lines. 

-> read_lsd()

Class: mask_obj

-> read_mask()

[ ] Class: region_obj()

-> read_region()
-> default_region()


Class: spectrum_obj

-> read_spectrum()

Class: line_list_obj

-> read_VALD()


[ ] Put the rvfit into LSD?
[ ] Add a rvfit to spectrum (putting the line into a lsd object and calling rvfit?)


[ ] Clean Mask needs some more structure


# Command line tools:


## GUI only tools

- cleanMaskUI
[ ] Can we run this from an open python shell
[ ] Make the dymanic use of LSDpy better ? (it wants to write file...)


## CommandLineScripts/ 
- Bz
- rvfit
- makeMask


## Structure of the repository at the end
[ ] Check first how to make pip work for installation. 

repo / specpolFlow / all of the source files
     / CommandLineScripts
     / Tutorials

[ ] arrange the __init__ file the way we want it. 
