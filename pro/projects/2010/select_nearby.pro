pro select_nearby
; jm10jun15 - select candidate nearby galaxies for camp/2010

    allcat = read_rc3()
    ra = 15*im_hms2dec(allcat.ra)
    dec = im_hms2dec(allcat.dec)

    keep = where((dec gt 5.0) and (ra gt 180.0) and $
      (ra lt 200.0) and (allcat.bmag gt 0.0) and $
      (allcat.r25 gt 0.0),nobj)
    cat = allcat[keep]
    ra = ra[keep]
    dec = dec[keep]

    djs_plot, ra, dec, psym=6, xsty=3, ysty=3
    
    im_plothist, cat.per*24, bin=0.5 

    airmass_plots, '2010-06-18', cat._raj2000, cat._dej2000, $
      obj=cat.gcvs, psname='select_cvs.ps', /bigpost
    
    
stop    

return
end
    
