pro rotationcurve

    datapath = getenv('AYCAMP_DATA')+'2011/bok/28jun11/spec1d/'
    outpath = getenv('AYCAMP_DATA')+'2011/projects/rotationcurve/'

    psfile = outpath+'rotationcurve.ps'
    scale = 1D15
    xr = [6500,6700]
    yr = [-0.5,12]
    
    aycamp_plotconfig, 0, pos, psfile=psfile, height=5.0, $
      xmargin=[1.1,0.6], width=6.8

    aycamp_plotspec, datapath+'NGC3448_28jun11.0292_aper1.fits', $
      xsty=1, ysty=1, xrange=xr, yrange=yr, psym=10, scale=scale, $
      xtitle='Observed Wavelength (\AA)', $
      ytitle='Flux (10^{-16} erg s^{-1} cm^{-2} \AA^{-1})', $
      position=pos, /baseline, thick=7, title='NGC 3448'
    
    aycamp_plotspec, datapath+'NGC3448_28jun11.0292_aper2.fits', $
      position=pos, /noerase, color='red', scale=scale, psym=10, $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      xrange=xr, yrange=yr, /baseline, thick=7

    aycamp_plotspec, datapath+'NGC3448_28jun11.0292_aper5.fits', $
      position=pos, /noerase, color='blue', scale=scale, psym=10, $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      xrange=xr, yrange=yr, /baseline, thick=7
    
    aycamp_plotspec, datapath+'NGC3448_28jun11.0292_aper6.fits', $
      position=pos, /noerase, color='dark green', scale=scale, psym=10, $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      xrange=xr, yrange=yr, /baseline, thick=7

    aycamp_plotspec, datapath+'NGC3448_28jun11.0292_aper7.fits', $
      position=pos, /noerase, color='orange', scale=scale, psym=10, $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      xrange=xr, yrange=yr, /baseline, thick=7

    im_legend, ['R = -1400 lyr','R = -700 lyr','R = +760 lyr','R = +1330 lyr','R = +2400 lyr'], $
      /left, /top, box=0, textcolor=['black','red','blue','dark green','orange']
    
    aycamp_plotconfig, psfile=psfile, /psclose, /pdf


; rotation curve
    radius = [-1394.8,-696.2,760.8,1331.5,2409.3]
    vel = [1494.6,1494.2,1340.2,1332.97,1292.6]
    
    psfile = 'scurve.ps'
    aycamp_plotconfig, 0, pos, psfile=psfile, height=5.0
    djs_plot, radius, vel-1350.0, psym=8, xrange=[-2500,2500], yrange=[1200,1550]-1350, $
      xsty=3, ysty=3, symsize=3, ytitle='Velocity (km s^{-1})', $
      xtitle='Distance from Center (lyr)', charsize=2
    djs_oplot, !x.crange, [0,0], line=5, thick=5
    djs_oplot, [0,0], !y.crange, line=5, thick=5
    aycamp_plotconfig, psfile=psfile, /psclose, /pdf

return
end
    
