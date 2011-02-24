;+
; NAME:
;   AYCAMP_FITS_CUBE()
;
; PURPOSE:
;   Read many FITS images into a structure cube.
;
; INPUTS:
;   fitslist - string array of FITS file names
;   ext - extension number to read
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   cube  - structure cube with an image and header field 
;
; COMMENTS:
;   Very slow for a large number of images!  Also the maximum header
;   size is hard-coded to 1000 elements.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2001 Aug 12, U of A, written
;   jm02nov1uofa - added additional error checking
;   jm02nov27uofa - generalized to read in one-dimensional FITS
;                   files
;   jm04nov05uofa - updated the HEADER string array from 500 to
;                   1000 elements
;   jm10jun01ucsd - ported to AYCAMP/IMG
;-

function aycamp_fits_cube, fitslist, ext=ext

    nfits = n_elements(fitslist)
    if nfits eq 0L then begin
       doc_library, 'aycamp_fits_cube'
       return, -1
    endif

    if (size(fitslist[0],/type) ne 7) or $
      (strmatch(fitslist[0],'*.fits*') eq 0) then begin
       print, 'File list must be type string FITS files'
       return, -1
    endif

    if (n_elements(ext) eq 0) then ext = 0
    
; read in the first image
    if (file_test(fitslist[0]) eq 0) then begin
       splog, 'FITS file '+fitslist[0]+' not found'
       return, -1
    endif

    image = mrdfits(fitslist[0],ext,head,/fscale,/silent)
    ndim = size(image,/n_dim)
    imsz = size(image,/dim)
    if (ndim ne 1) and (ndim ne 2) then begin
       splog, 'FITS file must be one- or two-dimensional!'
       return, -1
    endif
    
; create the data cube, assuming all the images are the same size;
; store the headers in an array of pointers
    if (ndim eq 1) then $
      template = {fname: '', image: make_array(imsz[0],/float), header: strarr(1000)} else $
      template = {fname: '', image: make_array(imsz[0],imsz[1],/float), header: strarr(1000)}
    cube = replicate(template,nfits)

    cube.fname = fitslist
    cube[0].image = image
    cube[0].header[0:n_elements(head)-1L] = head
    for ii = 1L, nfits-1L do begin
       if (file_test(fitslist[ii]) eq 0) then begin
          splog, 'FITS file '+fitslist[i]+' not found'
          return, -1
       endif

       image = mrdfits(fitslist[ii],ext,head,/fscale,/silent)
       szt = size(image,/dim)
       sdim = size(image,/n_dim)
       if (sdim ne ndim) then begin
          splog, 'FITS files have incompatible dimensions!'
          return, -1
       endif
       if (sdim eq 1) then begin
          if (imsz[0] ne szt[0]) then begin
             splog, 'FITS files are not the same size!'
             return, -1
          endif
       endif else begin
          if ((imsz[0] ne szt[0]) or (imsz[1] ne szt[1])) then begin
             splog, 'FITS files are not the same size!'
             return, -1
          endif
       endelse 
       cube[ii].image = image
       cube[ii].header[0:n_elements(head)-1L] = head
    endfor

return, cube
end
