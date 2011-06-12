;+
; NAME:
;   REDUCE_KP09M_2009
;
; PURPOSE:
;   Reduce all the KPNO/0.9-meter data from AYCamp/2009.
;;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;   inspect - 
;   mkbias - 
;   mkdflat - 
;   assign_calib - 
;   proc - 
;   clobber - 
;
; OUTPUTS: 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 May 31, UCSD
;-

pro reduce_kp09m_2009, datapath, night=night, inspect=inspect, $
  mkbias=mkbias, mkdflat=mkdflat, assign_calib=assign_calib, $
  proc=proc, unpack_projects=unpack_projects, clobber=clobber

    if (n_elements(datapath) eq 0) then datapath = getenv('HOME')+'/'+$
      'home/research/data/aycamp/2009/0.9m/iredux/'
    splog, 'Data path '+datapath
    tel = 'kp09m'
    ccd = 'S2KB'

    if (n_elements(night) eq 0) then night = $
      ['21jun2009','22jun2009','23jun2009','24jun2009']
    nnight = n_elements(night)
    structfile = night+'_imgstrct.fits'

; #########################
; preliminary reductions

; inspect all the files, removing problematic or irrelevant exposures 
    if keyword_set(inspect) then begin 
       for inight = 0, nnight-1 do begin
          pushd, datapath+night[inight]
          aycamp_img_strct, allstruct, ccd=ccd, tel=tel, $
            outfil='all_'+structfile[inight], clobber=clobber
          case night[inight] of
             '21jun2009': keep = lindgen(n_elements(allstruct)) ; all OK!
             '22jun2009': keep = where($
               (strmatch(allstruct.img_root,'*.0001.*',/fold) eq 0) and $ ; autoguider problem
               (strmatch(allstruct.img_root,'*.0016.*',/fold) eq 0) and $ ; dome problem
               (strmatch(allstruct.obj,'*service*',/fold) eq 0))          ; service observing
             '23jun2009': keep = where($
               (strmatch(allstruct.obj,'*service*',/fold) eq 0)) ; service observing
             '24jun2009': keep = where($
               (strmatch(allstruct.obj,'*service*',/fold) eq 0)) ; service observing
             else: message, 'Update me!'
          endcase
          struct = allstruct[keep]
          struct_print, struct_trimtags(struct,select=['img_root',$
            'obj','type','exp','am','filter','naxis1','naxis2']), $
            file=repstr(structfile[inight],'.fits','.txt')
          mwrfits, struct, structfile[inight], /create
          popd
       endfor
    endif

; build the master bias frames, including the 512x512 mini-biases 
    if keyword_set(mkbias) then begin
       for inight = 0, nnight-1 do begin
          pushd, datapath+night[inight]
          struct = mrdfits(structfile[inight],1)
          aycamp_img_mkbias, struct, clobber=clobber
          allbias = file_search('Bias/Bias?x?.fits',count=nbias)
          for ibias = 0, nbias-1 do begin
             info = file_info(allbias[ibias])
             if (info.symlink eq 0) then begin
                splog, 'Reading '+allbias[ibias]
                bias2048 = mrdfits(allbias[ibias],0,hdr,/silent)
                bias512 = bias2048[768:1279,768:1279] ; see the observing log
                outfile = repstr(allbias[ibias],'.fits','_512.fits')
                splog, 'Writing '+outfile
                mwrfits, bias512, outfile, hdr, /create
             endif 
          endfor 
          popd
       endfor
    endif 

; build the master dome flat-fields, including the 512x512
; mini-flat-fields 
    if keyword_set(mkdflat) then begin
       for inight = 0, nnight-1 do begin
          pushd, datapath+night[inight]
          struct = mrdfits(structfile[inight],1)
          aycamp_img_mkdflat, struct, clobber=clobber
          alldflat = file_search('Flats/DFlat_*?x?.fits',count=ndflat)
          for idflat = 0, ndflat-1 do begin
             info = file_info(alldflat[idflat])
             if (info.symlink eq 0) then begin
                splog, 'Reading '+alldflat[idflat]
                dflat2048 = mrdfits(alldflat[idflat],0,hdr,/silent)
                dflat512 = dflat2048[768:1279,768:1279] ; see the observing log
                outfile = repstr(alldflat[idflat],'.fits','_512.fits')
                splog, 'Writing '+outfile
                mwrfits, dflat512, outfile, hdr, /create
             endif
          endfor 
          popd
       endfor 
    endif

; #########################
; assign each science image a bias frame and dome flat-field
    if keyword_set(assign_calib) then begin
       pushd, datapath
       for inight = 0, nnight-1 do begin
          biasfiles = file_search(night[inight]+'/Bias/Bias*.fits',count=nbias)
          if (nbias gt 0) then begin
             biasinfo1 = aycamp_forage(biasfiles)
             if (n_elements(biasinfo) eq 0) then $
               biasinfo = biasinfo1 else $
                 biasinfo = [biasinfo,biasinfo1]
          endif
          dflatfiles = file_search(night[inight]+'/Flats/DFlat*.fits',count=ndflat)
          if (ndflat gt 0) then begin
             dflatinfo1 = aycamp_forage(dflatfiles)
             if (n_elements(dflatinfo) eq 0) then $
               dflatinfo = dflatinfo1 else $
                 dflatinfo = [dflatinfo,dflatinfo1]
          endif 
       endfor 

       for inight = 0, nnight-1 do begin
          thisstructfile = night[inight]+'/'+night[inight]+'_imgstrct.fits'
          struct = mrdfits(thisstructfile,1,/silent)
          obj = where(((struct.type eq 'OBJ') or (struct.type eq 'STD')) and $
            (struct.flg_anly ne 0),nobj)
          for iobj = 0, nobj-1 do begin
             bin = strtrim(struct[obj[iobj]].cbin,2)+'x'+$
               strtrim(struct[obj[iobj]].rbin,2)
             filt = strtrim(struct[obj[iobj]].filter,2)
             naxis2 = strtrim(struct[obj[iobj]].naxis2,2)
; do the necessary bias frames and dome flats exist?  if not, then
; create a soft link pointing to the appropriate file from another
; night (if it exists)
             if (naxis2 eq 512) then suffix = '_512' else suffix = ''
             biasfile = 'Bias/Bias'+bin+suffix+'.fits'
             dflatfile = 'Flats/DFlat_'+filt+bin+suffix+'.fits'
             if (file_test(night[inight]+'/'+biasfile) eq 0) then begin
                if (n_elements(biasinfo) ne 0) then begin
                   these = where((biasinfo.cbin eq struct[obj[iobj]].cbin) and $
                     (biasinfo.rbin eq struct[obj[iobj]].rbin) and $
                     (biasinfo.naxis2 eq struct[obj[iobj]].naxis2),nthese)
                   if (nthese eq 0) then begin
                      splog, 'No master bias frame found for binning '+$
                        bin+' and axis size '+strtrim(naxis2,2)+'x'+strtrim(naxis2,2)+'!'
                   endif else begin
                      mindiff = min(abs(biasinfo[these].jdmean-struct[obj[iobj]].date),mindx)
                      pushd, night[inight]
                      spawn, 'ln -sfv ../../'+strtrim(biasinfo[these[mindx]].file,2)+$
                        ' '+file_basename(biasfile), /sh
                      popd
                      struct[obj[iobj]].bias_fil = biasfile ; assigned bias frame
                   endelse
                endif
             endif else struct[obj[iobj]].bias_fil = biasfile ; assigned bias frame

; find flats with the right binning and filter name
             if (file_test(night[inight]+'/'+dflatfile) eq 0) then begin
                if (n_elements(dflatinfo) ne 0) then begin
                   these = where((dflatinfo.cbin eq struct[obj[iobj]].cbin) and $
                     (dflatinfo.rbin eq struct[obj[iobj]].rbin) and $
                     (strtrim(dflatinfo.filter,2) eq filt) and $
                     (dflatinfo.naxis2 eq struct[obj[iobj]].naxis2),nthese)
                   if (nthese eq 0) then begin
                      splog, 'No master dome-flat found for binning '+$
                        bin+', axis size '+strtrim(naxis2,2)+'x'+$
                        strtrim(naxis2,2)+' and filter '+filt+'!'
                   endif else begin
                      mindiff = min(abs(dflatinfo[these].jdmean-struct[obj[iobj]].date),mindx)
                      pushd, night[inight]
                      spawn, 'ln -sfv ../../'+strtrim(dflatinfo[these[mindx]].file,2)+$
                        ' '+file_basename(dflatfile), /sh
                      popd
                      struct[obj[iobj]].flat_fil = dflatfile ; assigned dome-flat field
                   endelse
                endif
             endif else struct[obj[iobj]].flat_fil = dflatfile ; assigned dome-flat field
          endfor ; close object loop
          splog, '### Night '+night[inight]
          splog, 'Updating '+thisstructfile
          mwrfits, struct, thisstructfile, /create
          struct_print, struct_trimtags(struct,select=$
            ['img_root','obj','type','exp',$
            'am','filter','naxis1','naxis2']), $
            file=repstr(thisstructfile,'.fits','.txt')
; test info for the screen
          struct_print, struct_trimtags(struct[obj],select=$
            ['obj','filter','naxis2','cbin','rbin','bias_fil','flat_fil'])
          print
       endfor    ; close night loop
       popd
    endif    

; #########################
; overscan-subtract, bias-subtract, and flat-field all the science
; images
    if keyword_set(proc) then begin
       for inight = 0, nnight-1 do begin
          pushd, datapath+night[inight]
          struct = mrdfits(structfile[inight],1)
          obj = where(((struct.type eq 'OBJ') or (struct.type eq 'STD')) and $
            (struct.flg_anly ne 0),nobj)
;         obj = obj[0] & nobj = 1
; check for missing bias frames and dome flats
          miss = where(strtrim(struct[obj].bias_fil eq '') or $
            strtrim(struct[obj].flat_fil eq ''),nmiss)
          if (nmiss ne 0) then begin
             splog, 'The following objects are missing bias and/or dome flat-fields!'
             stop
          endif
          aycamp_img_proc, struct[obj], outroot=outroot, clobber=clobber
          popd
       endfor 
    endif 

; #########################
; unpack the various imaging 
    if keyword_set(unpack_projects) then begin
       if (file_test(datapath+'baboquivari',/dir) eq 0) then $
         spawn, 'mkdir -p '+datapath+'baboquivari', /sh
       if (file_test(datapath+'landolt',/dir) eq 0) then $
         spawn, 'mkdir -p '+datapath+'landolt', /sh
;      for inight = 1, nnight-1 do begin
       for inight = 0, nnight-1 do begin
          struct = mrdfits(datapath+night[inight]+'/'+night[inight]+$
            '_imgstrct.fits',1,/silent)
          obj = where(((struct.type eq 'OBJ') or (struct.type eq 'STD')) and $
            (struct.flg_anly ne 0),nobj)
          struct = struct[obj]
          struct_print, struct_trimtags(struct,select=['img_root','obj','filter'])
; standard-star fields (don't need the weight maps)
          std = where(strmatch(struct.obj,'*landolt*',/fold) or $
            strmatch(struct.obj,'*standard*',/fold),nstd)
          if (nstd ne 0) then begin
             for ii = 0, nstd-1 do begin
                file = strtrim(struct[std[ii]].img_root,2)
                spawn, 'rsync -auv '+datapath+night[inight]+$
                  '/Final/f_'+file+'.gz '+datapath+'landolt/', /sh
;               spawn, 'rsync -auv '+datapath+night[inight]+$
;                 '/Final/f_'+repstr(file,'.fits','.weight.fits.gz')+' '+$
;                 datapath+'landolt/', /sh
             endfor
          endif
; Baboquivari
          babo = where(strmatch(struct.obj,'*baboquivari*',/fold),nbabo)
          if (nbabo ne 0) then begin
             for ii = 0, nbabo-1 do begin
                file = strtrim(struct[babo[ii]].img_root,2)
                spawn, 'rsync -auv '+datapath+night[inight]+$
                  '/Final/f_'+file+'.gz '+datapath+'baboquivari/', /sh
                spawn, 'rsync -auv '+datapath+night[inight]+$
                  '/Final/f_'+repstr(file,'.fits','.weight.fits.gz')+' '+$
                  datapath+'baboquivari/', /sh
             endfor
          endif 
       endfor
    endif

return
end
    
