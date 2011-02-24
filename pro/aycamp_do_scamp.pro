pro aycamp_do_scamp, imagelist, suffix=suffix
; SCAMP

    if (n_elements(suffix) eq 0) then suffix = ''
    
    datapath = file_dirname(imagelist[0])+'/'
    catlist = repstr(imagelist,'.fits','.cat')
    nimage = n_elements(imagelist)

; initialize the scamp configuration parameters
    http = ''
;   http = 'http://sdss.physics.nyu.edu/ioannis/'+project+'/'
    configfile = datapath+'scamp'+suffix+'.config'
    xslscamp = http+'scamp.xsl'

    config = init_scamp_config()

; boost the weights
    config.astref_weight = '100.0'
    config.crossid_radius = 2.0 ; this is important

; the OBJECT names are different between 2003 & 2006, so use this
; header tag to force scamp to define a different *astrometric*
; instrument between the two epochs; use the FILTER keyword to define
; the three (BVR) photometric instruments
    config.astrinstru_key = 'FILTER'
    config.photinstru_key = 'FILTER'
    config.magzero_key = 'MAGZERO'
    config.save_refcatalog = 'N'
    config.refout_catpath = datapath
    config.mergedoutcat_type = 'NONE'

    config.checkplot_type = strjoin(['ASTR_CHI2','ASTR_INTERROR1D','ASTR_INTERROR2D',$
      'ASTR_REFERROR1D','ASTR_REFERROR2D','DISTORTION','FGROUPS','PHOT_ERROR','PHOT_ZPCORR',$
      'PHOT_ZPCORR3D'],',')
    config.checkplot_name = strjoin(datapath+['astr_chi2','astr_interror1d',$
      'astr_interror2d','astr_referror1d','referror2d','distort','fgroups',$
      'psphot_error','phot_zpcorr','phot_zpcorr3d']+suffix,',')

    config.xml_name = datapath+'scamp'+suffix+'.xml'
    config.xsl_url = xslscamp

    t0 = systime(1)
    maxiter = 1

; think carefully before increasing DEGREE; for example DEGREE>3 when
; making a single mosaic results in DOOM!       

    for iter = 0, maxiter do begin 
       case iter of
          0: begin
             config.distort_degrees = '1'
             config.mosaic_type = 'FIX_FOCALPLANE' ; 'UNCHANGED' ; 'LOOSE'
             config.pixscale_maxerr = '1.1'
             config.position_maxerr = '3.0'
             config.posangle_maxerr = '2.0'
             config.aheader_suffix = '.ahead'
          end
          1: begin
             config.distort_degrees = '3,3'
             config.mosaic_type = 'FIX_FOCALPLANE' ; 'UNCHANGED' ; 'LOOSE'
             config.pixscale_maxerr = '1.1'        ; '1.1'
             config.position_maxerr = '0.5'        ; '0.5'
             config.posangle_maxerr = '2.0'        ; '1.0'
             config.aheader_suffix = '.head'
          end
          else: begin
             config.distort_degrees = '5,5'
             config.mosaic_type = 'FIX_FOCALPLANE' ; 'UNCHANGED' ; 'LOOSE'
             config.pixscale_maxerr = '1.1'
             config.position_maxerr = '0.5'
             config.posangle_maxerr = '1.0'
             config.aheader_suffix = '.head'
          end
       endcase

       mwrfits, config, configfile+'.fits', /create
       im_scamp, catlist, config, silent=silent
    endfor 
    splog, 'Total time to run scamp = ', (systime(1)-t0)/60.0, ' minutes.'

return    
end 

