;+
; NAME:
;   AYCAMP_EMPTY_STRUCTURE()
;
; PURPOSE:
;   Empty a data structure semi-intelligently.
;
; INPUTS: 
;   oldstruct - old data structure
;
; OPTIONAL INPUTS: 
;   ncopies - number of times to copy the output structure 
;   empty_value - replace non-string fields with this value
;   empty_string - replace string fields with this string 
;
; KEYWORD PARAMETERS: 
;   extra - extra keywords for STRUCT_TRIMTAGS()
;
; OUTPUTS: 
;   newstruct - empties structure
;
; OPTIONAL OUTPUTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2005 Apr 07, U of A
;-

function aycamp_empty_structure, oldstruct, ncopies=ncopies, $
  empty_value=empty_value, empty_string=empty_string, $
  _extra=extra

    newstruct = oldstruct[0]
    struct_assign, {junk: 0}, newstruct

    if (n_elements(empty_value) ne 0) then $
      for itag = 0, n_tags(newstruct)-1 do $
        if size(newstruct.(itag),/tname) ne 'STRING' then $
          newstruct.(itag) = empty_value
    if (n_elements(empty_string) ne 0) then $
      for itag = 0, n_tags(newstruct)-1 do $
        if size(newstruct.(itag),/tname) eq 'STRING' then $
          newstruct.(itag) = empty_string
    
    newstruct = im_struct_trimtags(newstruct,_extra=extra)
    if (n_elements(ncopies) ne 0) then begin
       if (ncopies gt 1) then newstruct = replicate(newstruct,ncopies)
    endif

return, newstruct
end    
