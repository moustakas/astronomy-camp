#!/usr/bin/env python

from pyraf import iraf

# Load the packages we might need.

iraf.noao(_doprint=0)
iraf.onedspec(_doprint=0)
iraf.twodspec(_doprint=0)
iraf.twodspec.apextract(_doprint=0)

# Make sure that the dispersion axis is in the header.

iraf.hedit("science_exposure_filename", DISPAXIS=1, add+)

# TODO: Bias subtraction? Flatfielding?

# Run the spectral extraction program.

iraf.twodspec.apextract.apall("science_exposure_filename", find=False, \
        recenter=False, resize=False)

# TODO: Make an interactive image of the 2d spectrum to use to identify the lines with the arc lamps?

# Now identify the spectral lines in the arc lamps. Need to replace l1 l2 with the range of rows that have the spectral lines.

iraf.identify("science_exposure_filename", section="line l1 l2")

# Tell the extracted spectrum what the wavelength solutions are.

iraf.hedit("science_exposure_filename.ms", "REFSPEC1", \
        "science_exposure_filename", add+)

iraf.dispcor("science-exposure-filename", \
        "calibrated-science-exposure-filename")

# TODO: Plot the extracted spectrum?
