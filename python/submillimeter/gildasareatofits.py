# tested with Enthought Python 2.7.1 |EPD 7.0-2 (32-bit), os-x 10.7.
''' The purpose of this script is to take an .area file from a CLASS On-The-Fly map and generate a fits file with coordinate headers.

Inputs:
The file name, the center position of the map and the target name.

started at ATC2011, revised at ATC II, 2012,  douglase@bu.edu
'''

#file = raw_input('filename:')
#target_name=raw_input('target name:')
#RA=raw_input('RA, decimal:')
#DEC=raw_input('DEC, decimal:')
#s = raw_input('GILDAS CLASS Filename--> ')

#import packages
import numpy
import pyfits
import time

# example input criteria:
file = "sdd_fb12.aac_020_800x800.area"

target_name="M13, Eagle Nebula"
dir="/Users/edouglas/Documents/astronomycamp/atc2012/"
s=dir+file

RA =274.7  #target position, degrees RA,  18:19:32.9 -13:47:33.5
DEC = -13.7790277778

def areatofits(file,s,target_name,RA,DEC):
    area=numpy.genfromtxt(s, comments="!")

    numel = numpy.shape(area)[0]
    nx=0
    ny=0

    xvals=[area[0,1]]
    yvals=[area[0,2]]

#an inefficent loop through the values to find the dimensions of the array, ok b/c only 2 vectors:
    for i in range(numel-1):
        if area[i,1] not in xvals:
            xvals.append(area[i,1])        
        if area[i,2] not in yvals:
            yvals.append(area[i,2])



    vals=numpy.array(area[:,3])
    #test for rectangle:
    if numpy.size(xvals)*numpy.size(yvals) == numpy.size(vals):
        image=vals.reshape(numpy.size(xvals),numpy.size(yvals))
    else:
        print("rows or columns are not all the same size")

#rotate image as appropriate:
    image=numpy.rot90(image,k=3)
#image=numpy.fliplr(image)
#image=numpy.flipud(image)

#define header values:
    hdu = pyfits.PrimaryHDU(image)
    hdulist = pyfits.HDUList([hdu])

    CRPIX1  = numpy.size(xvals)/2.0 # 10.5 #n_el(RA)/2
    CRPIX2  =  numpy.size(yvals)/2.0 #10.5 #n_el(DEC)/2
    CRVAL1  = RA #274.7  #target position, degrees RA,  18:19:32.9 -13:47:33.5
    CRVAL2  = DEC #-13.7790277778  #62.0978889  #target position, deg. DEC
    CDELT1  = (numpy.max(xvals) - numpy.min(xvals))/(2*CRPIX1)/3600.0 #30/3600.0 <-WHERE DID THIS VALUE COME FROM?
        #resolution =(range vals)/nbins*(deg/arcsec)
    CDELT2  =  (numpy.max(yvals) - numpy.min(yvals))/(2*CRPIX1)/3600.0 #30/3600.0   #resolution
    CTYPE1  = 'RA---TAN' 
    CTYPE2  = 'DEC--TAN' 
    hdu.header.update('target', target_name, 'target name')
    hdu.header.update('CRPIX1', CRPIX1)
    hdu.header.update('CRPIX2', CRPIX2)
    hdu.header.update('CRVAL1', CRVAL1)
    hdu.header.update('CRVAL2', CRVAL2)
    hdu.header.update('CDELT1', CDELT1)
    hdu.header.update('CDELT2', CDELT2)
    hdu.header.update('CTYPE1', CTYPE1)
    hdu.header.update('CTYPE2', CTYPE2)

    hdulist.writeto(file+time.asctime()+'.fits')


areatofits(file,s,target_name,RA,DEC)
