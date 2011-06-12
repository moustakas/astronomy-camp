;+ 
; NAME:
;   AYCAMP_IMG_PROC
;
; PURPOSE:
;   Processes a set of images (overscan-subtract, bias-subtract, trim,
;   and flat-field).
;
; INPUTS:
;   esi -- dimg_strct defining the images of interest
;   img    -- intarr of images to process
;   flatnm -- Flat root to use (e.g. 'Flats/SkyFltN')
;
; RETURNS:
;
; OUTPUTS:
;   Processed image (e.g. Final/f_ccd001.fits)
;
; OPTIONAL KEYWORDS:
;  OUTPTH =  Output directory (default = 'Final/')
;  INTER  =  Interactive OV fitting
;  DELOV  =  Delete ov files when through
;  MSK    =  Default name for mask file
;  NOGAIN =  Do not apply the gain
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_imgproc, dimg, STDS, 'Flats/SkyFltN'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   30-Jul-2002 Written by JXP
;-

pro aycamp_img_proc, struct, outroot=outroot, nozip=nozip, $
  clobber=clobber
    
    if (n_elements(struct) eq 0) then begin 
       doc_library, 'aycamp_img_proc'
       return
    endif 
    
    img = where((struct.flg_anly ne 0) and $
      (strtrim(struct.flat_fil,2) ne ''),nimg)
    if (nimg eq 0) then begin
       splog, 'No good images!'
       return
    endif

    if (n_elements(outroot) eq 0) then outroot = 'Final/f_'
    imfiles = strtrim(outroot+struct[img].img_root,2)
    ivarfiles = repstr(imfiles,'.fits','.weight.fits')

; sort by flat-field name
    allflatfiles = strtrim(struct[img].flat_fil,2)
    flatfiles = allflatfiles[uniq(allflatfiles,sort(allflatfiles))]
    nflat = n_elements(flatfiles)

    for ii = 0, nflat-1 do begin
       these = where(flatfiles[ii] eq allflatfiles,nthese)
       dflat = mrdfits(flatfiles[ii],/silent)
       mask = dflat gt 0.75 ; bad-pixel mask [0=bad]

; loop on each individual image so that we don't have to write
; out each overscan-subtracted image
       for jj = 0, nthese-1 do begin
          if file_test(imfiles[these[jj]]+'.gz') and $
            (keyword_set(clobber) eq 0) then begin
             splog, 'Final file '+imfiles[these[jj]]+' exists...skipping'
          endif else begin
             aycamp_img_overscan, struct[img[these[jj]]], rawsub=rawimage, $
               rawivar=rawivar, hdr=hdr, /nowrite
             biasfile = strtrim(struct[img[these[jj]]].bias_fil,2)
             bias = mrdfits(biasfile,/silent)
             image = (rawimage - bias)/(dflat+(dflat eq 0))*(dflat ne 0)*mask
             ivar = rawivar*dflat^2*mask
;            test = image*0.0
;            sz = size(image,/dim)
;            medwidth = 301L
;            for irow = 0L, sz[1]-1L do for icol = 0L, sz[0]-1L do test[icol,irow] = image[icol,irow] - $
;              median(image[(icol-medwidth/2L)>0L:(icol+medwidth/2L)<(sz[0]-1L),irow])
;; reject cosmic rays and interpolate over bad pixels for cosmetic
;; purposes
;            sky, image[where(image gt 0.0)], skymode, skysig, /silent
;            reject_cr, image-skymode, ivar, [0.496,0.246], rejects, $
;              nrejects=nrejects, c2fudge=c2fudge, niter=10
;            ivar[rejects] = 0
             image = djs_maskinterp(image,ivar eq 0,iaxis=0,/const)
; update the header, making sure the filter name appears in the
; header; also add the basic astrometric header
             sxaddpar, hdr, 'BIASCOR', file_basename(biasfile), $
               ' bias frame', before='HISTORY'
             sxaddpar, hdr, 'FLATCOR', file_basename(flatfiles[ii]), $
               ' dome flat-field', before='HISTORY'
             sxaddhist, "'Bias-subtracted "+hogg_iso_date()+"'", hdr
             sxaddhist, "'Flat-fielded "+hogg_iso_date()+"'", hdr

             filt = sxpar(hdr,'FILTER',count=fcheck)
             if (fcheck eq 0) then sxaddpar, hdr, 'FILTER', $
               strtrim(struct[img[these[jj]]].filter,2), ' filter name', $
               before='HISTORY'

             sz = size(image,/dim)
             pixscale = struct[img[these[jj]]].pixscale/3600D
             astr = hogg_make_astr(struct[img[these[jj]]].ra,$
               struct[img[these[jj]]].dec,sz[0]*pixscale,$
               sz[1]*pixscale,pixscale=pixscale,orientation=orientation)
             putast, hdr, astr
             hdr = hdr[where((strmatch(strcompress(hdr,/remove),'*HISTORY*PUTAST*',/fold) eq 0))]
             sxaddhist, "'Preliminary WCS astrometry added "+hogg_iso_date()+"'", hdr

; write out
             splog, 'Writing '+imfiles[these[jj]]
             mwrfits, float(image), imfiles[these[jj]], hdr, /silent, /create
             if (keyword_set(nozip) eq 0) then spawn, 'gzip -f '+imfiles[these[jj]], /sh

             splog, 'Writing '+ivarfiles[these[jj]]
             mwrfits, float(ivar), ivarfiles[these[jj]], hdr, /silent, /create
             if (keyword_set(nozip) eq 0) then spawn, 'gzip -f '+ivarfiles[these[jj]], /sh
          endelse
       endfor 
    endfor

return
end
