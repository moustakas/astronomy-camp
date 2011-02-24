pro aycamp_do_swarp, imagelist, object=object
; initialize the swarp configuration parameters

    nimage = n_elements(imagelist)
    datapath = file_dirname(imagelist[0])+'/'

    http = ''
;   http = 'http://sdss.physics.nyu.edu/ioannis/'+project+'/'
    xslswarp = http+'swarp.xsl'

    config = init_swarp_config()
    config.center_type = 'MANUAL' 
    config.pixelscale_type = 'MANUAL' 

    config.blank_badpixels = 'N'
    config.interpolate = 'N'
    config.write_fileinfo = 'Y'
    config.celestial_type = 'EQUATORIAL'
    config.copy_keywords = 'OBJECT,FILTER'
    config.xml_name = datapath+'swarp.xml'
    config.xsl_url = xslswarp

    config.write_xml = 'Y'
    config.subtract_back = 'Y'
    config.resampling_type = 'LANCZOS3'
    config.weight_type = 'MAP_WEIGHT'
    config.combine_type = 'AVERAGE'
;   config.combine_type = 'WEIGHTED'
    
; build a mosaic separately for each filter       
    info = aycamp_forage(imagelist)
    allfilt = strtrim(info.filter,2)
    filt = allfilt[uniq(allfilt,sort(allfilt))]
    nfilt = n_elements(filt)

; read the scamp headers to determine the central coordinates of the
; output mosaic, across all the filters
    headlist = repstr(imagelist,'.fits','.head')
    crval1 = dblarr(nimage)
    crval2 = dblarr(nimage)
    for ii = 0, nimage-1 do begin
       hdr = djs_readlines(headlist[ii])
       crval1[ii] = sxpar(hdr,'crval1')
       crval2[ii] = sxpar(hdr,'crval2')
    endfor
    cd2_2 = sxpar(hdr,'CD2_2')
    ra = im_dec2hms(djs_mean(crval1/15D),/colon)
    dec = im_dec2hms(djs_mean(crval2),/colon)

    for ifilt = 0, nfilt-1 do begin
; read the scamp headers to determine the central coordinates of the
; output mosaic, across all the filters
       these = where(filt[ifilt] eq allfilt,nthese)
       imlist = strtrim(imagelist[these],2)
       weightlist = repstr(imlist,'.fits','.weight.fits')

       config.center = ra+','+dec
       config.pixel_scale = string(fix(1000D*cd2_2*3600D)/1000D,format='(F5.3)')
       config.image_size = strtrim(info[0].naxis1,2)+','+$
         strtrim(info[0].naxis2,2)

       configfile = datapath+'swarp_'+filt[ifilt]+'.config'
       config.imageout_name = datapath+object+'_'+filt[ifilt]+'.fits'
       config.weightout_name = repstr(config.imageout_name,'.fits','.weight.fits')

       config.weight_image = strjoin(weightlist,',')

       splog, 'Building '+strtrim(config.imageout_name,2)
       mwrfits, config, configfile+'.fits', /create
       im_swarp, imlist, config, silent=silent
    endfor

return
end

