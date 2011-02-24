function remove_cosmics, image, weight
    reject_cr, image, weight, [0.496,0.246], rejects, $
      nrejects=nrejects, c2fudge=c2fudge, niter=10
    mask = image*0.0
    mask[rejects] = 1.0
    clean_image = djs_maskinterp(image,mask,iaxis=0,/const)
return, clean_image
end

pro pne_mosaics, sextractor=sextractor, scamp=scamp, $
  swarp=swarp, test=test, subtract=subtract, jpeg=jpeg
; jm10jun18ucsd - build mosaics!

    project = 'pne'
;   datapath = '/mount/moon1/ioannis/aycamp/'+project+'/'
;   datapath = getenv('HOME')+'/'+project+'/'
    datapath = getenv('AYCAMP_DATA')+'2010/kp09m/projects/'+project+'/'
    
; parse the file list
    allfiles = file_search(datapath+'*_????_*.fits')
    isweight = where(strmatch(allfiles,'*weight*'),comp=isimage)
    imagelist = allfiles[isimage]
    nimage = n_elements(imagelist)

; just do the 2048 images!    
    info = aycamp_forage(imagelist)
    keep = where(info.naxis1 eq 2048,nimage)
    imagelist = imagelist[keep]
    
; generate SE catalogs
    if keyword_set(sextractor) then aycamp_do_sextractor, imagelist

; astrometry; group by mosaic size
    if keyword_set(scamp) then begin
       info = aycamp_forage(imagelist)
       is512 = where(info.naxis1 eq 512,comp=is2048)
;      aycamp_do_scamp, imagelist[is512], suffix='_512'
       aycamp_do_scamp, imagelist[is2048], suffix='_2048'
    endif

; build the mosaics; work on each galaxy separately
    basefile = file_basename(imagelist)
    allobj = strarr(nimage)
    for ii = 0, nimage-1 do allobj[ii] = $
      strmid(basefile[ii],0,strpos(basefile[ii],'_'))
    obj = allobj[uniq(allobj,sort(allobj))]
    nobj = n_elements(obj)
    
    if keyword_set(swarp) then begin
       for iobj = 0, nobj-1 do begin
          splog, 'Working on '+obj[iobj]
          these = where(obj[iobj] eq allobj)
          aycamp_do_swarp, imagelist[these], object=obj[iobj]
       endfor
    endif
       
    if keyword_set(subtract) then begin
;      hb_on_weight = mrdfits(datapath+'Abell39_Hb.weight.fits')
;      hb_on = remove_cosmics(hb_on,hb_on_weight)


       hb_on = mrdfits(datapath+'Abell72_Hb.fits')
       ha_on = mrdfits(datapath+'Abell72_Ha6580.fits')
       ha_off = mrdfits(datapath+'Abell72_Ha6660.fits')
       
       ha = ha_on - 0.925*ha_off
       hb = hb_on - 0.925*ha_off
       nw_rgb_make, ha, (ha+hb)/2.0, hb, name='Abell72.jpeg', $
         scales=[1.0,1.0,1.1], nonlinearity=3.0, rebinfactor=2, $
         quality=90

       
       

stop       
       
       hb_on = mrdfits(datapath+'Abell39_Hb.fits')
       ha_on = mrdfits(datapath+'Abell39_Ha6580.fits')
       ha_off = mrdfits(datapath+'Abell39_Ha6660.fits')
       
       ha = ha_on - 0.925*ha_off
       hb = hb_on; - 0.925*ha_off
       nw_rgb_make, ha, (ha+hb)/2.0, hb, name='Abell39.jpeg', scales=[2.0,1.0,1.0], $
         nonlinearity=3.0, rebinfactor=2, quality=90

       
stop
    endif
       
;; JPEG
;    if keyword_set(jpeg) then begin
;       aycamp_do_jpeg, datapath+'M101_R.fits', datapath+'M101_B.fits', $
;         datapath+'M101_B.fits', jpegfile='~/public_html/tmp/M101_BR.jpeg', $
;         scales=[2.0,1.0,1.0]
;    endif
       
return
end
