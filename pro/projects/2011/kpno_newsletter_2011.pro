pro kpno_newsletter_2011
; jm11jul20ucsd - build figures for the KPNO Newsletter

; ---------------------------------------------------------------------------
; make a QAplot of WISE1810
    datapath = getenv('AYCAMP_DATA')+'2011/bok/projects/wise/'    

    flux = mrdfits(datapath+'WISE1810+69_f.fits',0,hh) & wave = make_wave(hh)
    flux = djs_maskinterp(flux,(wave gt 5565) and (wave lt 5585))

    zobj = 0.00769
    wave = wave/(1+zobj)
    flux = flux/interpol(flux,wave,5500.0)
    
    psfile = getenv('AYCAMP_DATA')+'2011/bok/projects/newsletter/wise.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0

    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, position=pos, $
      xrange=[3580,6960], yrange=[0.5,6.5], xtitle='Rest Wavelength (\AA)', $
      ytitle='Relative Flux'
    djs_oplot, wave, flux, psym=10
    legend, 'WISE-1', /left, /top, box=0

;   wline = [6732.0,6716.0,$
;     6584.0,6548.0,6562.8]
;   wtext = textoidl(['[SII]','[SII]','[NII]','[NII]','H\alpha'])
;   im_lineid_plot, wave, flux, wline, wtext, psym=10, xsty=1, ysty=1, position=pos, $
;     xrange=[3580,6990], yrange=[0.5,6.5], xtitle='Observed Wavelength (\AA)', $
;     ytitle='Relative Flux', lcharthick=!p.charthick, ymargin=[
    
    im_plotconfig, psfile=psfile, /psclose;, /pdf

stop    
    
; ---------------------------------------------------------------------------
; make a single plot with the three supernovae we typed, plus 2011dn (M51)
    datapath = getenv('AYCAMP_DATA')+'2011/bok/projects/sne/'    
    dnflux = mrdfits(datapath+'sn2011dn_2011_06_23_f.fits',0,hh) & dnwave = make_wave(hh)
    dvflux = mrdfits(datapath+'sn2011dv_f.fits',0,hh) & dvwave = make_wave(hh)
    dwflux = mrdfits(datapath+'sn2011dw_f.fits',0,hh) & dwwave = make_wave(hh)
    dhflux = mrdfits(datapath+'sn2011dh_f.fits',0,hh) & dhwave = make_wave(hh)

; clean up sky-subtraction residuals
    dvflux = djs_maskinterp(dvflux,(dvwave gt 5565) and (dvwave lt 5584))
    dwflux = djs_maskinterp(dwflux,((dwwave gt 5565) and (dwwave lt 5584)) or $
      ((dwwave gt 6293) and (dwwave lt 6307)))

; normalize
    dnflux = dnflux/interpol(dnflux,dnwave,5500.0)
    dvflux = dvflux/interpol(dvflux,dvwave,5500.0)
    dwflux = dwflux/interpol(dwflux,dwwave,5500.0)
    dhflux = dhflux/interpol(dhflux,dhwave,5500.0)
    
    psfile = getenv('AYCAMP_DATA')+'2011/bok/projects/newsletter/sne.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0

    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, position=pos, $
      xrange=[3580,6990], yrange=[0.5,6.5], xtitle='Observed Wavelength (\AA)', $
      ytitle='Flux (arbitrary units)'
    dnflux = dnflux+4.0 & djs_oplot, dnwave, dnflux, psym=10
    dwflux = dwflux+3.0 & djs_oplot, dwwave, dwflux, psym=10
    dvflux = dvflux+1.8 & djs_oplot, dvwave, dvflux, psym=10
    dhflux = dhflux+0.7 & djs_oplot, dhwave, dhflux, psym=10

    xyouts, 6400.0, interpol(dnflux,dnwave,6400.0)+0.3, 'SN 2011dn', align=0.5, charsize=1.4
    xyouts, 6400.0, interpol(dwflux,dwwave,6400.0)+0.3, 'SN 2011dw', align=0.5, charsize=1.4
    xyouts, 6400.0, interpol(dvflux,dvwave,6400.0)+0.3, 'SN 2011dv', align=0.5, charsize=1.4
    xyouts, 6400.0, interpol(dhflux,dhwave,6400.0)+0.6, 'SN 2011dh', align=0.5, charsize=1.4
    
    im_plotconfig, psfile=psfile, /psclose;, /pdf
    
return
end
    
