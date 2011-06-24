pro select_cvs_2011
; jm11jun15ucsd - select some candidate CVs for camp/2011

    allcat = mrdfits(getenv('CATALOGS_DIR')+'/cataclysmics/catalog_of_cvs.fits.gz',1)
    keep = where((allcat.orb_per gt 0.0) and (allcat.orb_per*24 lt 10.0) and $
      (allcat._dej2000 gt 5.0) and (allcat._raj2000 gt 180.0) and $
      (allcat._raj2000 lt 300.0) and (allcat.mag1 gt 0.0) and $
      (allcat.mag1 lt 15),nobj)
    cat = allcat[keep]

    djs_plot, cat._raj2000, cat._dej2000, psym=6, xsty=3, ysty=3, sym=3
    im_plothist, cat.orb_per*24, bin=0.5 

    airmass_plots, '2011-06-25', cat._raj2000, cat._dej2000, $
      obj=cat.name, psname='select_cvs.ps', /bigpost, /pdf
    
    
stop    

return
end
    
