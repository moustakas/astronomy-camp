pro select_sao
; jm10jun15 - select SAO stars

    jj = read_sao()

    keep = where((jj._raj2000/15.0 gt 11.0) and (jj._raj2000/15.0 lt 15.0) and $
      (jj._dej2000 gt 20.0) and (jj._dej2000 lt 70.0) and $
;     (jj.vmag gt 8.0) and (strtrim(jj.sptype,2) ne ''))
    cat = jj[keep]
    cat = cat[sort(cat._raj2000)]

    this = where(strmatch(cat.sptype,'*A0*'))
    struct_print, struct_trimtags(cat[this],select=[$
      '_raj2000','_dej2000','vmag','sptype',$
      'ra2000','de2000'])
    niceprint, cat[this].dm, im_dec2hms(cat[this]._raj2000/15.0,/col), $
      im_dec2hms(cat[this]._dej2000,/col), cat[this].sptype, cat[this].vmag
    
    
    
    

    
    djs_plot, cat._raj2000, cat._dej2000, psym=6, xsty=3, ysty=3
    
    
    im_plothist, cat.per*24, bin=0.5 

    airmass_plots, '2010-06-18', cat._raj2000, cat._dej2000, $
      obj=cat.gcvs, psname='select_cvs.ps', /bigpost
    
    
stop    

return
end
    
