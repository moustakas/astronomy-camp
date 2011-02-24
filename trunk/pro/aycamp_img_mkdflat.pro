;+
; NAME:
;   AYCAMP_IMG_MKDFLAT
;
; PURPOSE:
;   Build the coadded master dome flat(s).
;
; INPUTS: 
;   struct - AYCAMP/IMG inspection structure
;   mmem - maximum memory to use in combining zero flats 
;
; KEYWORD PARAMETERS: 
;   unity - make a flat of all ones
;   clobber - overwrite existing files
;
; OUTPUTS: 
;   Overscan-subtracted files and the master dome flats.
;
; COMMENTS:
;   Separate coadded bias frames are built for every combination of
;   binning and filter.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 May 31, UCSD - based on XIDL code written by
;     J. Prochaska & collaborators; ported into the 
;-

pro aycamp_img_mkdflat, struct, ovroot=ovroot, outroot=outroot, $
  mmem=mmem, unity=unity, clobber=clobber 

    if (n_elements(struct) eq 0) then begin 
       doc_library, 'aycamp_img_mkdflat'
       return
    endif 

    if (n_elements(ovroot) eq 0) then ovroot = 'OV/ov_'
    if (n_elements(outroot) eq 0) then outroot = 'Flats/DFlat_'

; group by filter name and binning
    gd = where((struct.flg_anly NE 0) and $
      (strtrim(struct.type,2) ne 'ZRO') and $
      (strtrim(struct.type,2) ne 'DRK'),ngd)
    if (ngd eq 0) then message, 'Fix me'
      
    allfilt = strtrim(struct[gd].filter,2)
    filt = allfilt[uniq(allfilt,sort(allfilt))]
    nfilt = n_elements(filt)

    bintot = struct[gd].cbin + 10*struct[gd].rbin
    bintyp = bintot[uniq(bintot,sort(bintot))]
    nbin = n_elements(bintyp)

    for ii = 0, nfilt-1 do begin
       for jj = 0, nbin-1 do begin
          splog, 'Building master dome flat for filter '+$
            filt[ii]+', binning '+strtrim(bintyp[jj],2)
          row = bintyp[jj]/10
          col = bintyp[jj] - row*10
          outfil = outroot+filt[ii]+strtrim(col,2)+'x'+strtrim(row,2)+'.fits'
          if file_test(outfil) and (keyword_set(clobber) eq 0) then begin
             splog, 'Output file '+outfil+' exists; use /CLOBBER' 
          endif else begin
; grab the relevant dome flats
             these = where((struct.type EQ 'DFT') AND (struct.rbin EQ row) and $
               (struct.cbin eq col) and(struct.flg_anly ne 0) and $
               (strtrim(struct.filter,2) eq filt[ii]),nthese)
             if (nthese eq 0) then begin
                splog, 'No dome flats found!'
             endif else begin
; overscan-subtract and then combine, unless /UNITY
                xhalf = struct[these[0]].naxis1/2
                yhalf = struct[these[0]].naxis2/2
                sstatsec = '['+$
                  string(xhalf*0.5,format='(I0)')+':'+$
                  string(xhalf*1.5,format='(I0)')+','+$
                  string(yhalf*0.5,format='(I0)')+':'+$
                  string(yhalf*1.5,format='(I0)')+']'
                statsec = xregtovec(sstatsec)
                if keyword_set(unity) then begin
                   img_dflat = 1.0+0.0*xmrdfits(ovroot+$
                     struct[these[0]].img_root,/fscale,/silent)
                endif else begin
                   aycamp_img_overscan, struct[these], /noivar, $
                     ovroot=ovroot, clobber=clobber
                   if (nthese gt 1) then begin
                      xcombine, ovroot+struct[these].img_root, $
                        img_dflat, head, fcomb=0, mmem=mmem, $
                        scale='MED', statsec=sstatsec
                   endif else begin
                      img_dflat = xmrdfits(ovroot+struct[these].img_root,/fscale,/silent)
                   endelse
                endelse
; normalize
                norm = djs_median(img_dflat[statsec[0]:statsec[1],$
                  statsec[2]:statsec[3]])
                if (norm le 0.0) then message, 'Problem here'
                img_dflat = temporary(img_dflat)/norm
; write out
                mkhdr, hdr, float(img_dflat)
                sxdelpar, hdr, 'COMMENT'
                sxdelpar, hdr, 'DATE'
                sxaddpar, hdr, 'RBIN', fix(row), ' row bin factor'
                sxaddpar, hdr, 'CBIN', fix(col), ' column bin factor'
                sxaddpar, hdr, 'FILTER', filt[ii], ' filter name'
                sxaddpar, hdr, 'NCOMBINE', fix(nthese), ' number of images combined'
                sxaddpar, hdr, 'JDMEAN', djs_mean(struct[these].date), $
                  ' mean Julian date'
                sxaddhist, "'Master dome flat created "+hogg_iso_date()+"'", hdr
                splog, 'Writing '+outfil
                mwrfits, float(img_dflat), outfil, hdr, /create, /silent
             endelse
          endelse
       endfor 
    endfor

return
end
