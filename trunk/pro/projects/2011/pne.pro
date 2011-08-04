pro pne

    datapath = getenv('AYCAMP_DATA')+'2011/bok/24jun11/spec1d/'
    psfile = 'pne.ps'
;   xr = [6500,6700]
;   yr = [-0.5,12]
    
    aycamp_plotconfig, 0, pos, psfile=psfile, height=5.0, $
      xmargin=[1.5,0.4], width=6.6

    aycamp_plotspec, datapath+'CatsEyeNebula_24jun11.0170_aper1.fits', $
      xsty=1, ysty=3, xrange=xr, yrange=yr, psym=10, scale=1D11, $
      xtitle='Observed Wavelength (\AA)', $
      ytitle='Flux (10^{-11} erg s^{-1} cm^{-2} \AA^{-1})', $
      position=pos, thick=5, title="Cat's Eye"

    aycamp_plotspec, datapath+'NGC6058_24jun11.0166_aper1.fits', $
      xsty=1, ysty=3, xrange=xr, yrange=yr, psym=10, scale=1D13, $
      xtitle='Observed Wavelength (\AA)', $
      ytitle='Flux (10^{-13} erg s^{-1} cm^{-2} \AA^{-1})', $
      position=pos, thick=5, title='NGC 6058'

    aycamp_plotspec, datapath+'RingNebula_24jun11.0173_aper3.fits', $
      xsty=1, ysty=3, xrange=xr, yrange=yr, psym=10, scale=1D12, $
      xtitle='Observed Wavelength (\AA)', $
      ytitle='Flux (10^{-12} erg s^{-1} cm^{-2} \AA^{-1})', $
      position=pos, thick=5, title='Ring Nebula'

    aycamp_plotspec, datapath+'dumbellnebulaonWD_24jun11.0176_aper1.fits', $
      xsty=1, ysty=3, xrange=xr, yrange=yr, psym=10, scale=1D13, $
      xtitle='Observed Wavelength (\AA)', $
      ytitle='Flux (10^{-13} erg s^{-1} cm^{-2} \AA^{-1})', $
      position=pos, thick=5, title='Dumbell Nebula - WD'

    aycamp_plotspec, datapath+'dumbellnebulaongas_24jun11.0177_aper2.fits', $
      xsty=1, ysty=3, xrange=xr, yrange=yr, psym=10, scale=1D13, $
      xtitle='Observed Wavelength (\AA)', $
      ytitle='Flux (10^{-13} erg s^{-1} cm^{-2} \AA^{-1})', $
      position=pos, thick=5, title='Dumbell Nebula - Gas'

    
    aycamp_plotconfig, psfile=psfile, /psclose, /pdf


return
end
    
