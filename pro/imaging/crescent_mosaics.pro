pro crescent_mosaics, sextractor=sextractor, scamp=scamp, $
  swarp=swarp, test=test, jpeg=jpeg
; jm10jun18ucsd - build mosaics!

    project = 'crescent'
    datapath = '/mount/moon1/ioannis/aycamp/'+project+'/'
;   datapath = getenv('HOME')+'/'+project+'/'
;   datapath = getenv('AYCAMP_DATA')+'2010/kp09m/projects/'+project+'/'
    
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
       
;; JPEG
;    if keyword_set(jpeg) then begin
;       aycamp_do_jpeg, datapath+'M101_R.fits', datapath+'M101_B.fits', $
;         datapath+'M101_B.fits', jpegfile='~/public_html/tmp/M101_BR.jpeg', $
;         scales=[2.0,1.0,1.0]
;    endif
       
return
end
