pro select_quasars
; jm10jun16 - select SDSS quasars

    allsdss = read_10schneider()

    cut1 = where((allsdss._raj2000/15.0 gt 14.0) and $
      (allsdss._raj2000/15.0 lt 22.0) and $
      (allsdss._dej2000 gt 20.0) and $
      (allsdss._dej2000 lt 70.0) and (allsdss.gmag gt 0.0) and $
      (allsdss.z gt 1.5 and allsdss.z lt 3.3))
    sdss = allsdss[cut1]

    plot, sdss.z, sdss.gmag, psym=6, ysty=3, symsize=0.2
    
    cut2 = where(sdss.gmag lt 17.5)
    final = sdss[cut2]
    djs_oplot, final.z, final.gmag, psym=6, color='red', symsize=0.2

    struct_print, struct_trimtags(final,select=['sdss',$
      'raj2000','dej2000','z','gmag','rmag']) 
    
stop    

return
end
    
