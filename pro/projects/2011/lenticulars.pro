function readspec, galaxy, wave=wave

    datapath = getenv('AYCAMP_DATA')+'2011/bok/'
    case galaxy of
       'M102': begin
          date = '28jun11'
          file = 'M102_28jun11.0306_aper1.fits'
       end
       'NGC5485': begin
          date = '28jun11'
          file = 'NGC5485_28jun11.0302_aper1.fits'
       end
       'NGC6548': begin
          date = '28jun11'
          file = 'NGC6548_28jun11.0304_aper1.fits'
       end
       'NGC5354': begin
          date = '28jun11'
          file = 'NGC5354_28jun11.0298_aper1.fits'          
       end
       'NGC5864': begin
          date = '28jun11'
          file = 'NGC5864_28jun11.0300_aper1.fits'         
       end
       'NGC6654': begin
          date = '28jun11'
          file = 'NGC6654_28jun11.0308_aper1.fits'
       end
       'NGC5195': begin
          date = '27jun11'
          file = 'NGC5195_27jun11.0229_aper1.fits'
       end
       'NGC5273': begin
          date = '27jun11'
          file = 'NGC5273_27jun11.0227_aper1.fits'
       end
       'NGC5353': begin
          date = '27jun11'
          file = 'NGC5353_27jun11.0231_aper1.fits'
       end
       'NGC5710': begin
          date = '27jun11'
          file = 'NGC5710_27jun11.0233_aper1.fits'
       end
    endcase

    flux = mrdfits(datapath+date+'/spec1d/'+file,0,hdr)
    wave = make_wave(hdr)
       
return, flux
end

pro lenticulars

    datapath = getenv('AYCAMP_DATA')+'2011/projects/lenticulars/'

; --------------------------------------------------    
; no bar    
    flux1_nobar = aycamp_readspec(datapath+'NGC5353_27jun11.0231_aper1.fits',wave=wave1_nobar)
    flux2_nobar = aycamp_readspec(datapath+'NGC5485_28jun11.0302_aper1.fits',wave=wave2_nobar)
    flux3_nobar = aycamp_readspec(datapath+'NGC5354_28jun11.0298_aper1.fits',wave=wave3_nobar)

    wave1_nobar = wave1_nobar/(1+0.007755)
    wave2_nobar = wave2_nobar/(1+0.006671)
    wave3_nobar = wave3_nobar/(1+0.008603)
    
    flux2_int_nobar = interpol(flux2_nobar,wave2_nobar,wave1_nobar)
    flux3_int_nobar = interpol(flux3_nobar,wave3_nobar,wave1_nobar)

    flux1_scale_nobar = flux1_nobar/interpol(flux1_nobar,wave1_nobar,5500.0)
    flux2_scale_nobar = flux2_int_nobar/interpol(flux2_int_nobar,wave1_nobar,5500.0)
    flux3_scale_nobar = flux3_int_nobar/interpol(flux3_int_nobar,wave1_nobar,5500.0)

    bigflux_nobar = [[flux1_scale_nobar],[flux2_scale_nobar],[flux3_scale_nobar]]
    avgflux_nobar = total(bigflux_nobar,2)/3.0

; --------------------------------------------------    
; bar    
    flux1_bar = aycamp_readspec(datapath+'NGC5864_28jun11.0300_aper1.fits',wave=wave1_bar)
    flux2_bar = aycamp_readspec(datapath+'NGC5195_27jun11.0229_aper1.fits',wave=wave2_bar)
    flux3_bar = aycamp_readspec(datapath+'NGC6548_28jun11.0304_aper1.fits',wave=wave3_bar)

;   wave1_bar = wave1_bar/(1+0.001551)
;   wave2_bar = wave2_bar/(1+0.006288)
;   wave3_bar = wave3_bar/(1+0.022362)
    
    flux2_int_bar = interpol(flux2_bar,wave2_bar,wave1_bar)
    flux3_int_bar = interpol(flux3_bar,wave3_bar,wave1_bar)

    flux1_scale_bar = flux1_bar/interpol(flux1_bar,wave1_bar,5500.0)
    flux2_scale_bar = flux2_int_bar/interpol(flux2_int_bar,wave1_bar,5500.0)
    flux3_scale_bar = flux3_int_bar/interpol(flux3_int_bar,wave1_bar,5500.0)

    bigflux_bar = [[flux1_scale_bar],[flux2_scale_bar],[flux3_scale_bar]]
    avgflux_bar = total(bigflux_bar,2)/3.0
    
    psfile = datapath+'lenticulars_avg.ps'
    aycamp_plotconfig, 6, pos, psfile=psfile, xmargin=[1.2,0.4], $
      width=6.9

    xrange = [3540.0,6890]

    djs_plot, wave1_nobar, flux1_scale_nobar, psym=10, xsty=3, ysty=3, $
      ytitle='Relative Flux', $
      position=pos[*,0], color='blue', xrange=xrange, $
      xtickname=replicate(' ',10), thick=5
    djs_oplot, wave1_nobar, flux2_scale_nobar, color='red', psym=10, thick=5
    djs_oplot, wave1_nobar, flux3_scale_nobar, color='green', psym=10, thick=5
    djs_oplot, wave1_nobar, avgflux_nobar, color='black', thick=5, psym=10
    legend, 'S0 Average', /right, /bottom, box=0

    djs_plot, wave1_bar, flux1_scale_bar, psym=10, xsty=3, ysty=3, $
      xtitle='Rest Wavelength (\AA)', ytitle='Relative Flux', $
      position=pos[*,1], color='blue', xrange=xrange, /noerase, $
      yrange=[0.1,1.3], thick=5
    legend, 'SB0 Average', /right, /bottom, box=0
      
    djs_oplot, wave1_bar, flux2_scale_bar, color='red', psym=10, thick=5
    djs_oplot, wave1_bar, flux3_scale_bar, color='green', psym=10, thick=5
    djs_oplot, wave1_bar, avgflux_bar, color='black', thick=5, psym=10

; compare both

    aycamp_plotconfig, 0, pos, xmargin=[1.2,0.4], height=1.5, width=6.9

    djs_plot, wave1_nobar, avgflux_nobar, psym=10, xsty=3, ysty=3, $
      ytitle='Relative Flux', xtitle='Rest Wavelength (\AA)', $
      position=pos, color='red', thick=5
    djs_oplot, wave1_bar, avgflux_bar, color='blue', thick=5, psym=10

    im_legend, ['S0 Average','SB0 Average'], /right, /bottom, box=0, $
      textcolor=['red','blue']

    aycamp_plotconfig, psfile=psfile, /psclose, /pdf
;   spawn, 'convert '+repstr(psfile,'.ps','.png'), /sh

stop    
    
    psfile = 'lenticulars.ps'
    aycamp_plotconfig, 6, pos, psfile=psfile, yspace=1.0
    
; pair 1    
    flux = readspec('M102',wave=wave)
    djs_plot, wave, smooth(1E15*flux,3), position=pos[*,0], psym=10, xsty=3, $
      ysty=3, xtitle='Observed Wavelength (\AA)', $
      ytitle='Flux (10^{-15} erg s^{-1} cm^{-2} \AA^{-1})'
    legend, ['M102','S0'], /left, /top, box=0
      
    flux = readspec('NGC5195',wave=wave)
    djs_plot, wave, 1E15*flux, position=pos[*,1], /noerase, $
      psym=10, xsty=3, ysty=3, xtitle='Observed Wavelength (\AA)', $
      ytitle='Flux (10^{-15} erg s^{-1} cm^{-2} \AA^{-1})'
    legend, ['NGC 5195','SB0'], /left, /top, box=0
      
; pair 2
    flux = readspec('NGC5354',wave=wave)
    djs_plot, wave, 1E15*flux, position=pos[*,0], psym=10, xsty=3, $
      ysty=3, xtitle='Observed Wavelength (\AA)', $
      ytitle='Flux (10^{-15} erg s^{-1} cm^{-2} \AA^{-1})'
    legend, ['NGC 5354','S0'], /left, /top, box=0
      
    flux = readspec('NGC5710',wave=wave)
    djs_plot, wave, 1E15*flux, position=pos[*,1], /noerase, $
      psym=10, xsty=3, ysty=3, xtitle='Observed Wavelength (\AA)', $
      ytitle='Flux (10^{-15} erg s^{-1} cm^{-2} \AA^{-1})'
    legend, ['NGC 5710','SB0'], /left, /top, box=0
      
; pair 3
    flux = readspec('NGC5485',wave=wave)
    djs_plot, wave, 1E15*flux, position=pos[*,0], psym=10, xsty=3, $
      ysty=3, xtitle='Observed Wavelength (\AA)', $
      ytitle='Flux (10^{-15} erg s^{-1} cm^{-2} \AA^{-1})'
    legend, ['NGC 5485','S0'], /left, /top, box=0
      
    flux = readspec('NGC6654',wave=wave)
    djs_plot, wave, 1E15*flux, position=pos[*,1], /noerase, $
      psym=10, xsty=3, ysty=3, xtitle='Observed Wavelength (\AA)', $
      ytitle='Flux (10^{-15} erg s^{-1} cm^{-2} \AA^{-1})'
    legend, ['NGC 6654','SB0'], /left, /top, box=0
      
; pair 4
    flux = readspec('NGC5353',wave=wave)
    djs_plot, wave, 1E15*flux, position=pos[*,0], psym=10, xsty=3, $
      ysty=3, xtitle='Observed Wavelength (\AA)', $
      ytitle='Flux (10^{-15} erg s^{-1} cm^{-2} \AA^{-1})'
    legend, ['NGC 5353','S0'], /left, /top, box=0
      
    flux = readspec('NGC5864',wave=wave)
    djs_plot, wave, 1E15*flux, position=pos[*,1], /noerase, $
      psym=10, xsty=3, ysty=3, xtitle='Observed Wavelength (\AA)', $
      ytitle='Flux (10^{-15} erg s^{-1} cm^{-2} \AA^{-1})'
    legend, ['NGC 5864','SB0'], /left, /top, box=0
      
; pair 5
    flux = readspec('NGC5273',wave=wave)
    djs_plot, wave, 1E15*flux, position=pos[*,0], psym=10, xsty=3, $
      ysty=3, xtitle='Observed Wavelength (\AA)', $
      ytitle='Flux (10^{-15} erg s^{-1} cm^{-2} \AA^{-1})'
    legend, ['NGC 5273','S0'], /left, /top, box=0
      
    flux = readspec('NGC6548',wave=wave)
    djs_plot, wave, 1E15*flux, position=pos[*,1], /noerase, $
      psym=10, xsty=3, ysty=3, xtitle='Observed Wavelength (\AA)', $
      ytitle='Flux (10^{-15} erg s^{-1} cm^{-2} \AA^{-1})'
    legend, ['NGC 6548','SB0'], /left, /top, box=0
      
    aycamp_plotconfig, psfile=psfile, /psclose, /pdf

    
    
stop    


return
end
    
