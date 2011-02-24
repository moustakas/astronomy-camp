;+ 
; NAME:
;   AYCAMP_IMG_OVERSCAN
;
; PURPOSE:
;    Overscan subtracts a list of images of a structure. 
;
; CALLING SEQUENCE:
;   
;   aycamp_img_overscan, struct, imgs, ORDR=, /ERASE, /INTER, /NOSV
;
; INPUTS:
;   struct -- dimg_strct defining the images of interest
;   imgs   -- integer array defining the images to process
;
; RETURNS:
;
; OUTPUTS:
;   ovimg - fits files in the dir OV
;
; OPTIONAL KEYWORDS:
;   /erase - Erase pre-existing ov files
;   /inter - Interactive OV fitting (uses x1dfit)
;   /nosv  -- Dont save the OV image!
;   ORDR=  -- Order to fit the overscan
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   aycamp_img_overscan, nght1_strct, ovimg
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   22-July-2001 Written by JXP
;-

pro rcos20_sbig_overscan, filename, rawsub, rawivar, hdr=hdr, $
  gain=gain, rnoise=rnoise, image=image, verbose=verbose, $
  bin=bin

;  no overscan for rows...
    if NOT keyword_set(image) then begin
;  look for filename
       filelist = lookforgzip(strtrim(filename,2))
       if filelist[0] EQ '' then begin
          print, 'Could not find file named ', filename
          return
       endif
       image = xmrdfits(filelist[0],0,hdr,/fscale,/silent)
       hdr = hdr[where(strcompress(hdr,/rem) ne '')]
    endif

    dims = size(image, /dim)
    ncol = dims[0]
    nrow = dims[1]
    ccdsum = strcompress(sxpar(hdr,'CCDSUM'), /rem)
    bin_col = long(strmid(ccdsum, 0, 1))
    bin_row = long(strmid(ccdsum, 1, 1))
    bin = [bin_col, bin_row]

    gain = double(sxpar(hdr, 'EGAIN'))
    rnoise = double(sxpar(hdr, 'RDNOISE'))

    temp_image = image*gain  ; [electrons]
    rawsub = temp_image
    rawivar = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) + rnoise^2) ; [electrons]

return
end

pro kp09m_s2kb_overscan, filename, rawsub, rawivar, hdr=hdr, $
  gain=gain, rnoise=rnoise, image=image, verbose=verbose, $
  bin=bin

;  no overscan for rows...
    if NOT keyword_set(image) then begin
;  look for filename
       filelist = lookforgzip(strtrim(filename,2))
       if filelist[0] EQ '' then begin
          print, 'Could not find file named ', filename
          return
       endif
       image = xmrdfits(filelist[0],0,hdr,/fscale,/silent)
       hdr = hdr[where(strcompress(hdr,/rem) ne '')]
    endif

    dims = size(image, /dim)
    ncol = dims[0]
    nrow = dims[1]
    ccdsum = strcompress(sxpar(hdr,'CCDSUM'), /rem)
    bin_col = long(strmid(ccdsum, 0, 1))
    bin_row = long(strmid(ccdsum, 1, 1))
    bin = [bin_col, bin_row]

    gain = double(sxpar(hdr, 'GAIN'))
    rnoise = double(sxpar(hdr, 'RDNOISE'))
    detector = sxpar(hdr, 'DETECTOR')
    datasec = strcompress(sxpar(hdr, 'DATASEC'), /rem)
    biassec = strcompress(sxpar(hdr, 'BIASSEC'), /rem)
    data_arr = long(strsplit(datasec, '[*:*,*:*]', /extract))
    bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))
    nbias = bias_arr[1]-bias_arr[0] + 1L
    
; these are small buffers to avoid using the overscan region too close 
; to data 
    oscan_buffer = 3L
    imagecol = data_arr[1]
    imagerow = data_arr[3]

    rawsub = fltarr(imagecol,imagerow)
    rawivar = fltarr(imagecol,imagerow)

    biascols = lindgen(nbias-oscan_buffer) + imagecol + oscan_buffer
    oscan1 = djs_avsigclip(image[biascols, 0:imagerow-1L], 1)
    sig_res = 7
    nhalf =  long(sig_res)*4L
    xkern = dindgen(2*nhalf+1)-nhalf
    kernel = gauss1(xkern, [0.0, sig_res, 1.0])
    oscan = convol(oscan1, kernel, /edge_truncate)
    osub = replicate(1, imagecol) # oscan
    temp_image = (image[0:imagecol-1L, 0:imagerow-1L] - osub)*gain  ; [electrons]
    rawsub[0:imagecol-1L, *] = temp_image
    rawivar[0:imagecol-1L, *] = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) + rnoise^2) ; [electrons]

; mark bad columns
    
; update the header
    sxaddpar, hdr, 'OVERSCAN', '['+$
      string(min(biascols)+1,format='(I0)')+':'+$
      string(max(biascols)+1,format='(I0)')+','+$
      string(1,format='(I0)')+':'+$
      string(imagerow,format='(I0)')+']', ' overscan region (median='+$
      strtrim(string(median(oscan),format='(F12.2)'),2)+')', before='HISTORY'
    sxaddhist, "'Overscan region subtracted "+hogg_iso_date()+"'", hdr

return
end

pro aycamp_img_overscan, struct, noivar=noivar, nowrite=nowrite, $
  ovroot=ovroot, rawsub=rawsub, rawivar=rawivar, hdr=hdr, $
  clobber=clobber

    if (n_elements(struct) eq 0) then begin 
       doc_library, 'aycamp_img_overscan'
       return
    endif
    nimage = n_elements(struct)
    
; input/output file names
    if (n_elements(ovroot) eq 0) then ovroot = 'OV/ov_'
    rawfiles = strtrim(struct.rootpth+struct.img_root,2)
    ovfiles = strtrim(ovroot+struct.img_root,2)

; which telescope/detector?
    tel = strtrim(struct.tel,2)
    ccd = strtrim(struct.ccd,2)

    for ii = 0, nimage-1 do begin
       if file_test(ovfiles[ii]+'.gz') and (keyword_set(clobber) eq 0) then begin
          splog, 'Overscan-subtracted file '+ovfiles[ii]+' exists...skipping'
       endif else begin
          case ccd[ii] of
             'S2KB': kp09m_s2kb_overscan, rawfiles[ii], $ ; KPNO/0.9-meter
               rawsub, rawivar, hdr=hdr 
             'SBIG': rcos20_sbig_overscan, rawfiles[ii], $ ; 20inch
               rawsub, rawivar, hdr=hdr 
             else: splog, 'Unknown instrument/telescope!'
          endcase
; write out
          if (keyword_set(nowrite) eq 0) then begin
             splog, 'Writing '+ovfiles[ii]
             mwrfits, float(rawsub), ovfiles[ii], hdr, /silent, /create
             if (keyword_set(noivar) eq 0) then $ ; do not write out ivar map
               mwrfits, float(rawivar), ovfiles[ii], hdr, /silent
             spawn, 'gzip -f '+ovfiles[ii], /sh
          endif
       endelse
    endfor

return
end
