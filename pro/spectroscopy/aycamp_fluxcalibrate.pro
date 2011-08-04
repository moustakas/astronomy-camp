;+
; NAME:
;   AYCAMP_FLUXCALIBRATE
;
; PURPOSE:
;   Flux-calibrate, trim crap pixels from the extremes of the
;   wavelength range, and convert to standard FITS format a spectrum
;   that has been processed by LONG_COADD or LONG_COMBSPEC.
;
; INPUTS: 
;   infile - 
;
; OPTIONAL INPUTS: 
;
;
; KEYWORD PARAMETERS: 
;
;
; OUTPUTS: 
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Jul 18, UCSD
;-

pro aycamp_fluxcalibrate, infile, outfile=outfile, ntrim=ntrim, $
  sensfuncfile=sensfuncfile, writetxt=writetxt, clobber=clobber, $
  gzip=gzip

    nspec = n_elements(infile)
    if (nspec eq 0) then begin
       doc_library, 'aycamp_fluxcalibrate'
       return
    endif

    nout = n_elements(outfile)
    if (nout eq 0) then begin
       outfile = repstr(repstr(infile,'.fits','_f.fits'),'.gz','')
    endif else begin
       if (nout ne nspec) then message, 'Dimensions of INFILE and OUTFILE must match!'
    endelse

    fluxcal = 0
    if (n_elements(sensfuncfile) ne 0) then begin
       if file_test(sensfuncfile) then fluxcal = 1 else $
         splog, 'Sensitivity function '+sensfuncfile+' not found!'
    endif
    
    if (n_elements(ntrim) eq 0) then ntrim = 10 ; [pixels]

    for iobj = 0, nspec-1 do begin
       if aycamp_file_test(outfile,clobber=clobber) then continue

       if fluxcal then begin
          long_fluxcal, infile[iobj], sensfunc=sensfuncfile, $
            outfil=outfile[iobj]
       endif else spawn, 'cp -f '+infile[iobj]+' '+outfile[iobj], /sh

; trim, etc.
       flux = mrdfits(outfile[iobj],0,hdr,/silent)
       ferr = mrdfits(outfile[iobj],1,/silent)
       wave = mrdfits(outfile[iobj],2,/silent)
       npix = n_elements(wave)
       flux = 1E-17*flux[ntrim:npix-ntrim-1]
       ferr = 1E-17*ferr[ntrim:npix-ntrim-1]
       wave = wave[ntrim:npix-ntrim-1]
       npix = n_elements(wave)

; interpolate over pixels affected by strong sky lines
;      mask = ((wave gt 5545.0) and (wave lt 5595.0)) or $
;        ((wave gt 6280.0) and (wave lt 6307.0)) or $
;        ((wave gt 6345.0) and (wave lt 6368.0))
;      flux = djs_maskinterp(flux,mask,wave,/const)
; rebin linearly in wavelength
       dwave = ceil(100D*(max(wave)-min(wave))/(npix-1))/100D
       newwave = dindgen((max(wave)-min(wave))/dwave+1)*dwave+min(wave)
       newflux = rebin_spectrum(flux,wave,newwave)
       newvar = rebin_spectrum(ferr^2,wave,newwave)
       newferr = sqrt(newvar*(newvar gt 0.0)) ; enforce positivity

; build the final header
       mkhdr, newhdr, float(newflux), /exten
       sxdelpar, newhdr, 'DATE'
       sxdelpar, newhdr, 'COMMENT'
       sxdelpar, newhdr, 'END'
       newhdr = [newhdr,[$
         hdr[where(strmatch(hdr,'OBJECT*'))],$
         hdr[where(strmatch(hdr,'RA *'))],$
         hdr[where(strmatch(hdr,'DEC *'))],$
         hdr[where(strmatch(hdr,'EQUINOX*'))],$
         hdr[where(strmatch(hdr,'FILENAME*'))],$
         hdr[where(strmatch(hdr,'DATE-OBS*'))],$
         hdr[where(strmatch(hdr,'TIME-OBS*'))],$
         hdr[where(strmatch(hdr,'UT *'))],$
         hdr[where(strmatch(hdr,'AIRMASS*'))],$
         hdr[where(strmatch(hdr,'NEXP *'))],$
         hdr[where(strmatch(hdr,'EXPT_AVG=*'))],$
         hdr[where(strmatch(hdr,'ROTANGLE*'))],$
         hdr[where(strmatch(hdr,'TELESCOP*'))],$
         hdr[where(strmatch(hdr,'OBSERVAT*'))],$
         hdr[where(strmatch(hdr,'OBSERVER*'))],$
         hdr[where(strmatch(hdr,'DISPERSE*'))],$
         hdr[where(strmatch(hdr,'APERTURE*'))]]]
       
       sxaddpar, newhdr, 'CRVAL1', min(newwave), ' wavelength at CRPIX1'
       sxaddpar, newhdr, 'CRPIX1', 1D, ' reference pixel number'
       sxaddpar, newhdr, 'CD1_1', dwave, ' dispersion [Angstrom/pixel]'
       sxaddpar, newhdr, 'CDELT1', dwave, ' dispersion [Angstrom/pixel]'
       sxaddpar, newhdr, 'CTYPE1', 'LINEAR', ' projection type'
       mwrfits, float(newflux), outfile[iobj], newhdr, /create
       mwrfits, float(newferr), outfile[iobj], newhdr
       if keyword_set(gzip) then spawn, 'gzip -f '+outfile[iobj], /sh

; optionally write out a text spectrum
       if keyword_set(writetxt) then begin
          txthdr = newhdr[where((strcompress(newhdr,/remove) ne '') and $
            (strmatch(newhdr,'*XTENSION*') eq 0) and $
            (strmatch(newhdr,'*SIMPLE*') eq 0) and $
            (strmatch(newhdr,'*BITPIX*') eq 0) and $
            (strmatch(newhdr,'*NAXIS *') eq 0) and $
            (strmatch(newhdr,'*EXTEND*') eq 0) and $
            (strmatch(newhdr,'*PCOUNT*') eq 0) and $
            (strmatch(newhdr,'*GCOUNT*') eq 0) and $
            (strmatch(newhdr,'*END *') eq 0))]
          txtfile = repstr(outfile[iobj],'.fits','.txt')
          openw, lun, txtfile, /get_lun
          for ii = 0, n_elements(txthdr)-1 do printf, lun, '# '+txthdr[ii]
          printf, lun, '# '
          aycamp_niceprintf, lun, wave, flux, ferr
          free_lun, lun
       endif
    endfor 

return
end

