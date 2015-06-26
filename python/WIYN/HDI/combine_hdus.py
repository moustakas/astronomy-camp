	
import astropy.io.fits as fits
import glob 

def combineHDUs(directory):
    '''takes a directory of HDI .fits files and 
    combines the imageHDU and PrimaryHDU headers
    into one new PrimaryHDU with the data set to
    the data from the imageHDU
    
    Parameters:
        directory:
            the directory that will be checked for .fits files 
            and populated with the original file name with _single_HDU.fits
            appended.
    '''
    filelist=glob.glob(directory+"/*.fits")

    for filename in filelist:
        print(filename)
        dataset=fits.open(filename)
        newheader=dataset[0].header+dataset[1].header
        newheader["HISTORY"]="2 element HDUList compressed into 1 primaryHDU"
        print(dataset[0].data)
        newfits=fits.HDUList([fits.PrimaryHDU(data=dataset[1].data,header=newheader)])
        newfits.writeto(filename[:-5]+"_single_HDU.fits",clobber=True)
        dataset.close()
