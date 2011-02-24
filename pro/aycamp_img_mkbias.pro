;+
; NAME:
;   AYCAMP_IMG_MKBIAS
;
; PURPOSE:
;   Build the coadded master bias frame(s).
;
; INPUTS: 
;   struct - AYCAMP/IMG inspection structure
;   mmem - maximum memory to use in combining dome flats 
;
; KEYWORD PARAMETERS: 
;   clobber - overwrite existing files
;
; OUTPUTS: 
;   Overscan-subtracted files and the master bias frames.
;
; COMMENTS:
;   Separate coadded bias frames are built for each different
;   binning. 
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 May 31, UCSD - based on XIDL code written by
;     J. Prochaska & collaborators; ported into the 
;-

pro aycamp_img_mkbias, struct, ovroot=ovroot, mmem=mmem, clobber=clobber

    if (n_elements(struct) eq 0) then begin 
       doc_library, 'aycamp_img_mkbias'
       return
    endif 
    if (n_elements(ovroot) eq 0) then ovroot = 'OV/ov_'

; find all binnings
    gd = where(struct.flg_anly NE 0)
    bintot = struct[gd].cbin + 10*struct[gd].rbin
    bintyp = bintot[uniq(bintot, sort(bintot))]
    nbin = n_elements(bintyp)

    for qq = 0, nbin-1 do begin
       splog, 'Building bias frame for binning '+strtrim(bintyp[qq],2)
       row = bintyp[qq]/10
       col = bintyp[qq] - row*10
       outfil = 'Bias/Bias'+strtrim(col,2)+'x'+strtrim(row,2)+'.fits'
       if file_test(outfil) and (keyword_set(clobber) eq 0) then begin
          splog, 'Output file '+outfil+' exists; use /CLOBBER' 
       endif else begin
; grab all the bias frames
          these = where((struct.type eq 'ZRO') and (struct.rbin eq row) and $
            (struct.cbin eq col) and (struct.flg_anly ne 0),nthese)
          if (nthese eq 0) then begin
             splog, 'No bias frames found!'
          endif else begin
; overscan-subtract and then combine
             aycamp_img_overscan, struct[these], /noivar, $
               ovroot=ovroot, clobber=clobber
             if (nthese gt 1) then begin
                xcombine, ovroot+struct[these].img_root, $
                  img_zro, head, fcomb=0, mmem=mmem
             endif else begin
                img_zro = xmrdfits(ovroot+struct[these].img_root, /fscale, /silent)
             endelse
; write out
             mkhdr, hdr, float(img_zro)
             sxdelpar, hdr, 'COMMENT'
             sxdelpar, hdr, 'DATE'
             sxaddpar, hdr, 'RBIN', fix(row), ' row bin factor'
             sxaddpar, hdr, 'CBIN', fix(col), ' column bin factor'
             sxaddpar, hdr, 'NCOMBINE', fix(nthese), ' number of images combined'
             sxaddpar, hdr, 'JDMEAN', djs_mean(struct[these].date), $
               ' mean Julian date'
             sxaddhist, "'Master bias frame created "+hogg_iso_date()+"'", hdr
             splog, 'Writing '+outfil
             mwrfits, float(img_zro), outfil, hdr, /create, /silent
          endelse
       endelse
    endfor 

return
end
