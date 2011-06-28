;+
; NAME:
;   AYCAMP_SENSFUNC
;
; PURPOSE:
;  Simple wrapper to build the sensitivity function.
;
; INPUTS: 
;
;
; OPTIONAL INPUTS: 
;
;
; KEYWORD PARAMETERS: 
;
;
; OUTPUTS: 
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Jun 27, UCSD 
;
;-

pro aycamp_sensfunc, stdfiles, stdnames, sensfuncfile=sensfuncfile, $
  nogrey=nogrey
          
; do the fit, masking the bluest and reddest pixels
    ncol = 1200 & nmask1 = 5 & nmask2 = 5
    inmask = intarr(ncol)+1
    inmask[0:nmask1-1] = 0
    inmask[ncol-nmask2-1:ncol-1] = 0

    sens = long_sensfunc(stdfiles,sensfuncfile,sensfit=sensfit,$
      std_name=stdnames,wave=wave,flux=flux,nogrey=nogrey,inmask=inmask,$
      /msk_balm)
    
return
end
