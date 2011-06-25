function aycamp_find_calspec, oldplan, radius=radius
; jm09dec17ucsd - given a PLAN, return the indices of CALSPEC
; standards 
; jm10jun03ucsd - ported into AYCAMP repository

    if (n_elements(radius) eq 0) then radius = 1.0 ; [degree]

    nobj = n_elements(oldplan)
    starinfo = replicate({starfile: '...', starname: '...'},nobj)
    newplan = struct_addtags(oldplan,starinfo)

; read the calspec database info file
    allstd = aycamp_rsex(getenv('AYCAMP_DIR')+'/data/calspec_info.txt')

; grab the coordinates of the standards from the header
    sci = where(strtrim(newplan.flavor,2) eq 'science',nsci)
    if (nsci eq 0) then return, newplan
    info = aycamp_forage('Raw/'+strtrim(newplan[sci].filename,2))
;   struct_print, struct_trimtags(info,sel=['object','ra','dec'])
    
; spherematch generously
    delvarx, m1, m2
    for ii = 0, n_elements(info)-1 do begin
       srcor, 15D*hms2dec(allstd.ra), hms2dec(allstd.dec), $
         15D*hms2dec(info[ii].ra), hms2dec(info[ii].dec), radius*3600.0, $
         m11, m21, /silent, spherical=2
       if (m11[0] ne -1) then begin
          if (n_elements(m1) eq 0) then begin
             m1 = m11
             m2 = ii
          endif else begin
             m1 = [m1,m11]
             m2 = [m2,ii]
          endelse
       endif 
    endfor
    
;   spherematch, 15.0*hms2dec(allstd.ra), hms2dec(allstd.dec), $
;     15.0*hms2dec(info.ra), hms2dec(info.dec), radius, $
;     m1, m2, max=0
    if (m2[0] ne -1) then begin
       newplan[sci[m2]].flavor = 'std'
       newplan[sci[m2]].starfile = repstr(allstd[m1].file,'.fits.gz','')
       newplan[sci[m2]].starname = allstd[m1].star
    endif
;   struct_print, newplan

return, newplan
end
    
