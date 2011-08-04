pro ngc4559
; jm11aug04ucsd - make some pretty plots for the NGC4559 rotation
; curve project 

    datapath = getenv('AYCAMP_DATA')+'2011/bok/projects/ngc4559/'

    psfile = datapath+'ngc4559_halpha.ps'
    scale = 1D15
    xr = [6520,6650]
    yr = [-0.5,5.0]
    
    aycamp_plotconfig, 0, pos, psfile=psfile, height=5.0, $
      xmargin=[1.1,0.6], width=6.8

    aycamp_plotspec, datapath+'NGC4559_aper1_f.fits', $
      xsty=1, ysty=1, xrange=xr, yrange=yr, psym=10, scale=scale, $
      xtitle='Observed Wavelength (\AA)', $
      ytitle='Flux (10^{-16} erg s^{-1} cm^{-2} \AA^{-1})', $
      position=pos, /baseline, thick=7, title='NGC 4559'
    
    aycamp_plotspec, datapath+'NGC4559_aper2_f.fits', $
      position=pos, /noerase, color='red', scale=scale, psym=10, $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      xrange=xr, yrange=yr, /baseline, thick=7, /overplot

    aycamp_plotspec, datapath+'NGC4559_aper3_f.fits', $
      position=pos, /noerase, color='blue', scale=scale, psym=10, $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      xrange=xr, yrange=yr, /baseline, thick=7, /overplot
    
    aycamp_plotspec, datapath+'NGC4559_aper4_f.fits', $
      position=pos, /noerase, color='dark green', scale=scale, psym=10, $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      xrange=xr, yrange=yr, /baseline, thick=7, /overplot

    aycamp_plotspec, datapath+'NGC4559_aper5_f.fits', $
      position=pos, /noerase, color='orange', scale=scale, psym=10, $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      xrange=xr, yrange=yr, /baseline, thick=7, /overplot

    legend, ['R = -4.12 kpc','R = -2.18 kpc','R = -0.61 kpc','R = Galaxy Center','R = +1.89 kpc'], $
      /left, /top, box=0, textcolor=djs_icolor(['black','red','blue','dark green','orange']), $
      margin=0, charsize=1.8
    
    aycamp_plotconfig, psfile=psfile, /psclose, /pdf

; rotation curve
    dist = 8.859E3 ; distance from NED [kpc]
    rad_pixel = [13,40,62,70.5,97] ; [pixels]
    rad_arcsec = (rad_pixel-70.5)*1.6666   ; [arcsec]
    rad_kpc = rad_arcsec/3600.0*!pi/180.0*dist ; [kpc]
    lam = [6580.9,6582.2,6583.0,6582.7,6584.3] ; H-alpha centroid [Angstrom]
    vel = 2.99E5*(lam-6562.85)/lam ; [km/s]
    niceprint, rad_pixel, rad_arcsec, rad_kpc, lam, vel
    
    psfile = 'ngc4559_rotationcurve.ps'
    aycamp_plotconfig, 0, pos, psfile=psfile, height=5.0
    djs_plot, rad_kpc, vel-900, psym=8, xrange=[-6,+6], yrange=150*[-1,1], $
      xsty=3, ysty=3, symsize=3, ytitle='Velocity (km s^{-1})', $
      xtitle='Distance from Center (kpc)', charsize=2
    djs_oplot, !x.crange, [0,0], line=5, thick=5
    djs_oplot, [0,0], !y.crange, line=5, thick=5
    aycamp_plotconfig, psfile=psfile, /psclose, /pdf

return
end
    
