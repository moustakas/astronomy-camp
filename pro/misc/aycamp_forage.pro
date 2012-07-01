;+
; NAME:
;   AYCAMP_FORAGE()
;
; PURPOSE:
;   Retrieve useful FITS header information. 
;
; INPUTS:
;   fitslist - list of FITS images [NFITS]
;   ext - extension number to read (default 0)
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   forage - data structure with lots of goodies [NFITS] 
;
; EXAMPLE:
;   IDL> info = aycamp_forage('*.fits')
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 May 30, UCSD
;-

function aycamp_forage, fitslist1, ext=ext

    nfits1 = n_elements(fitslist1)
    if (nfits1 eq 0L) then begin
       doc_library, 'kp09m_forage'
       return, -1
    endif

    fitslist = file_search(fitslist1,count=nfits)

    for jj = 0L, nfits-1L do begin
       if (file_test(fitslist[jj]) eq 0) then $
         splog, 'File '+fitslist[jj]+' not found'
       hh = headfits(fitslist[jj],ext=ext)
       str = aycamp_hdr2struct(hh)
       forage1 = struct_trimtags(str,except=$
         ['*HISTORY*','*COMMENT*','_END'])
; add tags here
       add = {file: fitslist[jj]}
       forage1 = struct_addtags(add,temporary(forage1))
       
       if (n_elements(forage) eq 0L) then $
         forage = aycamp_empty_structure(forage1,ncopies=nfits)
       forage[jj] = aycamp_struct_assign(forage1,forage[jj])
    endfor

return, reform(forage)
end    
