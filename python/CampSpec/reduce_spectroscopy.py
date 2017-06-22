#!/usr/bin/env python

import os

# Delete some directories/files from previous runs.

os.system("rm -rf login.cl pyraf database uparm")
os.system("mkiraf")

# Now load IRAF

import pyraf.iraf as iraf

# Load the packages we might need.

iraf.noao(_doprint=0)
iraf.onedspec(_doprint=0)
iraf.twodspec(_doprint=0)
iraf.apextract(_doprint=0)
iraf.unlearn(iraf.apall)

# The name of the science file.

filename = 'vega_9.3narrow.fit'
extracted_filename = 'vega_9.3narrow.ms.fits'
calibrated_filename = 'vega_9.3narrow.calib.fits'

# Delete previous results.

os.system("rm "+extracted_filename+" "+calibrated_filename)

# Make sure that the dispersion axis is in the header.

iraf.hedit(images=[filename], fields=["DISPAXIS"], value=["1"], add="Yes")

# Run the spectral extraction program.

iraf.apextract.setParam("dispaxis", "1")

iraf.apall(input=filename, find="No", recenter="No", resize="No")

# Now identify the spectral lines in the arc lamps. Need to replace l1 l2 
# with the range of rows that have the spectral lines.

iraf.identify(filename, section="line 265 285")

# Tell the extracted spectrum what the wavelength solutions are.

iraf.hedit(images=[extracted_filename], fields=["REFSPEC1"], \
        value=[filename], add="Yes")

iraf.dispcor(input=extracted_filename, output=calibrated_filename)

# Plot the extracted spectrum?

iraf.splot(calibrated_filename)

