pro wise_observing_prep, getdss=getdss, write_regions=write_regions
; jm11may24ucsd - prepare for megacam observing run of the DLS field

    wisepath = getenv('AYCAMP_DIR')+'/projects/wise/'
    
    cat = rsex(wisepath+'wise_candidates.sex')
    airmass_plots, '2011-06-25', 15D*im_hms2dec(cat.ra), im_hms2dec(cat.dec), $
      object=cat.name, obsname='kpno', psname='wise_airmass.ps', /post, /big

return
end
    
