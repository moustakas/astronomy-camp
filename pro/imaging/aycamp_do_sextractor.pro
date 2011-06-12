pro aycamp_do_sextractor, imagelist
; generate SE catalogs

    datapath = file_dirname(imagelist[0])+'/'
    
    sexpath = getenv('AYCAMP_DIR')+'/data/'
    sexconfig = sexpath+'default.sex'
    sexparam = sexpath+'default.sex.param'
    sexconv = sexpath+'default.conv'
    sexnnw = sexpath+'default.nnw'
    
    weightlist = repstr(imagelist,'.fits','.weight.fits')
    catlist = repstr(imagelist,'.fits','.cat')
    
    nimage = n_elements(imagelist)
    config = init_sex_config(nimage)
    configfile = datapath+'sex.config'
    
    config.catalog_name = catlist
    config.weight_image = weightlist
    config.parameters_name = sexparam
    config.filter_name = sexconv
    config.starnnw_name = sexnnw

    config.catalog_type = 'FITS_LDAC'
    config.detect_thresh = 5.0
    config.analysis_thresh = 5.0
    config.weight_type = 'MAP_WEIGHT'
    config.weight_gain = 'N'
    config.interp_type = 'NONE'

    config.seeing_fwhm = 1.5
    config.mag_zeropoint = 28.0        ; arbitrary
    config.checkimage_type = 'NONE'    ; SEGMENTATION
;   config.checkimage_name = seglist

    mwrfits, config, configfile+'.fits', /create
    im_sex, imagelist, config, silent=silent ; do not pass CONFIGFILE

return
end    

