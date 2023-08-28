[ ] Rename ioLSD to something better, and arrange the init for direct access to these functions?



Class: lsd_obj

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
