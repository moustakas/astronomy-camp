pro aycamp_plotspec, file, nsmooth=nsmooth, overplot=overplot, _extra=extra, $
  xrange=xrange, flux=flux, wave=wave, scale=scale1, baseline=baseline, $
  pos=pos, postscript=postscript, pdf=pdf, objname=objname
; jm11jun29ucsd - plot a 1D extracted spectrum
    
    if (n_elements(file) eq 0) then begin
       doc_library, 'aycamp_plotspec'
       return
    endif

    if (file_test(file) eq 0) then begin
       splog, 'Spectral file '+file+' not found!'
       return
    endif

    if (n_elements(scale1) eq 0) then scale = 1.0 else scale = scale1
    
    flux = mrdfits(file,0,hdr)
    wave = aycamp_make_wave(hdr)

    if keyword_set(baseline) then begin
       get_element, wave, xrange, xx
       flux = flux-djs_median(flux[xx[0]:xx[1]])
    endif

    if keyword_set(postscript) then begin
       psfile = repstr(file,'.fits','.ps')
       aycamp_plotconfig, 0, pos, psfile=psfile, height=4.5
    endif
    
    if keyword_set(overplot) then begin
       djs_oplot, wave, scale*flux, psym=10, _extra=extra, xrange=xrange
    endif else begin
       if (n_elements(scale1) ne 0) then yprefix = '10^{-'+strtrim(fix(alog10(scale)),2)+'} ' else $
         yprefix = ''
       ytitle = 'Flux ('+yprefix+'erg s^{-1} cm^{-2} \AA^{-1})'
       djs_plot, wave, scale*flux, xsty=3, ysty=3, psym=10, xrange=xrange, $
         xtitle='Observed Wavelength (\AA)', ytitle=ytitle, _extra=extra
       if (n_elements(objname) ne 0) then legend, objname, /right, /top, box=0
    endelse

    if keyword_set(postscript) then aycamp_plotconfig, /psclose, $
      psfile=psfile, pdf=pdf
    
    
return
end
    
