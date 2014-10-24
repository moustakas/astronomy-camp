;+
; NAME:
;   REDUCE_BOK_2011
; PURPOSE:
;   Process all the Bok B&C spectroscopy obtained during ATC2011. 
; INPUTS: 
;   night - optional input to process a specific night of data
; KEYWORD PARAMETERS: 
;   preproc - pre-process all the data: make the relevant directories;
;     repair bad pixels; fix headers; and write the FITS files to the
;     'Raw' directory
;   plan - build the planfile for each night
;   calib - process all the calibration data (biases, domeflats,
;     skyflats, and wavelength maps)
;   standards - process all the standard stars
;   science - process all the science frames
;   sensfunc - 
; OUTPUTS: 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Jun 15, UCSD
;-

; long_setflex, 'Science/sci-23jun11.0092.fits.gz', 'bok_sky_bcs_400.sav'

pro reduce_bok_2011, night, preproc=preproc, plan=plan, calib=calib, $
  standards=standards, science=science, sensfunc=sensfunc, chk=chk, $
  unpack_wise=unpack_wise, unpack_sne=unpack_sne, unpack_vyaqr=unpack_vyaqr, $
  unpack_ngc4559=unpack_ngc4559, unpack_lenticulars=unpack_lenticulars, $
  unpack_pne=unpack_pne, unpack_rotationcurve=unpack_rotationcurve, $
  clobber=clobber

    unpack_something = (keyword_set(unpack_wise) or keyword_set(unpack_sne) or $
        keyword_set(unpack_vyaqr) or keyword_set(unpack_ngc4559) or $
        keyword_set(unpack_lenticulars) or keyword_set(unpack_pne) or $
        keyword_set(unpack_rotationcurve)) ? 1 : 0

    datapath = getenv('AYCAMP_DATA')+'2011/bok/'
    projectpath = datapath+'projects/'
    if (file_test(projectpath,/dir) eq 0) then $
      spawn, 'mkdir -p '+projectpath, /sh
    
    badpixfile = getenv('AYCAMP_DIR')+'/data/bok_badpix.dat'
    if (file_test(badpixfile) eq 0) then message, $
      'Bad pixel file '+badpixfile+' not found'
;    sensfuncfile = datapath+'sensfunc_2011.fits'

    if (n_elements(night) eq 0) then night = ['22jun11','23jun11',$
      '24jun11','27jun11','28jun11','29jun11']
    nnight = n_elements(night)

; ##################################################
; pre-processing: read the data in the "rawdata" directory, repair bad
; pixels, remove cosmic rays, clean up headers, and move the spectra
; to the "Raw" subdirectory for further processing 
    if keyword_set(preproc) then begin
       for inight = 0, nnight-1 do begin
          if (file_test(datapath+night[inight],/dir) eq 0) then $
            spawn, 'mkdir '+night[inight], /sh
          pushd, datapath+night[inight]
; make some directories 
          if (file_test('spec1d',/dir) eq 0) then $
            spawn, 'mkdir '+'spec1d', /sh
          if (file_test('Raw',/dir) eq 0) then $
            spawn, 'mkdir '+'Raw', /sh
          if (file_test('Fspec',/dir) eq 0) then $
            spawn, 'mkdir '+'Fspec', /sh
          if (file_test('Science',/dir) eq 0) then $
            spawn, 'mkdir '+'Science', /sh
          
;         sensfuncfile = datapath+'sensfunc_'+night[inight]+'.fits'
          planfile = 'plan_'+night[inight]+'.par'
          calibplanfile = 'plan_calib_'+night[inight]+'.par'
          
          splog, 'Reading '+badpixfile
          readcol, badpixfile, x1, x2, y1, y2, comment='#', $
            format='L,L,L,L', /silent
          
          allfiles = file_search(datapath+'rawdata/'+night[inight]+'/*.fits*', $
              count=nspec)
          if (nspec eq 0) then message, 'No files found in '+datapath+'rawdata/'

          for iobj = 0, nspec-1 do begin
             outfile = repstr('Raw/'+file_basename(allfiles[iobj]),'.gz','')
             if file_test(outfile+'.gz') and (keyword_set(clobber) eq 0) $
                    then begin
                splog, 'Output file '+outfile+' exists; use /CLOBBER'
             endif else begin
                image = mrdfits(allfiles[iobj],0,hdr,/fscale,/silent)

                sxaddpar, hdr, 'INSTRUME', 'bcspeclamps', ' instrument name'
                sxaddpar, hdr, 'DISPERSE', '400/4889', ' disperser'
                sxaddpar, hdr, 'APERTURE', '2.5', ' slit width'

                if strmatch(allfiles[iobj],'*22jun11.004[5-6].*',/fold) then $
                    sxaddpar, hdr, 'APERTURE', '4.5'
                if strmatch(allfiles[iobj],'*22jun11.0052.*',/fold) then $
                    sxaddpar, hdr, 'APERTURE', '4.5'

                if strmatch(allfiles[iobj],'*23jun11.008[3-4].*',/fold) then $
                    sxaddpar, hdr, 'APERTURE', '4.5'
                if strmatch(allfiles[iobj],'*23jun11.0095.*',/fold) then $
                    sxaddpar, hdr, 'APERTURE', '4.5'
                if strmatch(allfiles[iobj],'*23jun11.0109.*',/fold) then $
                    sxaddpar, hdr, 'APERTURE', '4.5'

                if strmatch(allfiles[iobj],'*24jun11.0159*',/fold) then begin
                   sxaddpar, hdr, 'OBJECT', 'HZ44 4.5 slit'
                   sxaddpar, hdr, 'APERTURE', '4.5'
                endif
                if strmatch(allfiles[iobj],'*24jun11.0160*',/fold) then $
                    sxaddpar, hdr, 'OBJECT', 'HZ44 2.5 slit'
                if strmatch(allfiles[iobj],'*24jun11.0180.*',/fold) then $
                    sxaddpar, hdr, 'APERTURE', '4.5'
                if strmatch(allfiles[iobj],'*24jun11.0181.*',/fold) then $
                    sxaddpar, hdr, 'APERTURE', '4.5'

                if strmatch(allfiles[iobj],'*27jun11.0223.*',/fold) then $
                    sxaddpar, hdr, 'APERTURE', '4.5'
                if strmatch(allfiles[iobj],'*27jun11.0244.*',/fold) then $
                    sxaddpar, hdr, 'APERTURE', '4.5'
                if strmatch(allfiles[iobj],'*27jun11.025[2-3].*',/fold) then $
                    sxaddpar, hdr, 'APERTURE', '4.5'
                
                if strmatch(allfiles[iobj],'*28jun11.029[0-1].*',/fold) then $
                    sxaddpar, hdr, 'APERTURE', '4.5'
                if strmatch(allfiles[iobj],'*28jun11.0310.*',/fold) then $
                    sxaddpar, hdr, 'APERTURE', '4.5'

                type = sxpar(hdr,'imagetyp')
                if (strlowcase(strtrim(type,2)) eq 'object') then begin
                   dims = size(image,/dim)
                   badpixmask = image*0.0
                   for ipix = 0, n_elements(x1)-1 do $
                     badpixmask[(x1[ipix]-1)>0:(x2[ipix]-1)<(dims[0]-1),$
                     (y1[ipix]-1)>0:(y2[ipix]-1)<(dims[1]-1)] = 1
                   image = djs_maskinterp(image,badpixmask,iaxis=0,/const)
                endif
;               exp = sxpar(hdr,'exptime')
;               if (strlowcase(strtrim(type,2)) eq 'object') and (exp gt 300.0) then begin
;                  splog, 'Identifying cosmic rays '
;                  aycamp_la_cosmic, image, gain=-1.0, readn=0.0, $
;                    outlist=newim, sigclip=5.0
;                  image = newim
;               endif

; the area of the detector that was read out changed slightly over the
; course of the camp, so standardize that here
                ccdsize = strcompress(sxpar(hdr,'CCDSIZE'),/rem)
                datasec = strcompress(sxpar(hdr,'DATASEC'),/rem)
                biassec = strcompress(sxpar(hdr,'BIASSEC'),/rem)
                ccdsum = strcompress(sxpar(hdr, 'CCDSUM'), /rem)
                rowbin = fix(strmid(ccdsum,1,1))

                data_arr = long(strsplit(datasec,'[*:*,*:*]',/extract))
                ccd_arr = long(strsplit(ccdsize,'[*:*,*:*]',/extract))
                data_arr = [1,1200,1,124]
                
                bias_arr = long(strsplit(biassec,'[*:*,*:*]',/extract))
                new_data_arr = '['+strtrim(data_arr[0],2)+':'+$
                  strtrim(data_arr[1],2)+','+strtrim(data_arr[2],2)+':'+$
                  strtrim(data_arr[3],2)+']'
                new_bias_arr = '['+strtrim(bias_arr[0],2)+':'+$
                  strtrim(bias_arr[1],2)+','+strtrim(bias_arr[2],2)+':'+$
                  strtrim(bias_arr[3],2)+']'
                
                sxaddpar, hdr, 'DATASEC', new_data_arr
                sxaddpar, hdr, 'TRIMSEC', new_data_arr
                sxaddpar, hdr, 'AMPSEC', new_data_arr
                sxaddpar, hdr, 'BIASSEC', new_bias_arr
                
                aycamp_mwrfits, image[*,0:data_arr[3]-1], outfile, hdr, /clobber
             endelse
          endfor 
          popd
       endfor 
    endif 
       
; ##################################################
; make the plan files
    if keyword_set(plan) then begin
       for inight = 0, nnight-1 do begin
          planfile = 'plan_'+night[inight]+'.par'
          calibplanfile = 'plan_calib_'+night[inight]+'.par'

          pushd, datapath+night[inight]
          allfiles = file_search('Raw/*.fits*',count=nspec)
          if (nspec eq 0) then begin
             splog, 'No files found in '+'Raw/'
             continue
          endif
          long_plan, '*.fits.gz', 'Raw/', planfile=planfile
          old = yanny_readone(planfile,hdr=hdr,/anony)
          new = aycamp_find_calspec(old,radius=radius)
          twi = where(strmatch(new.target,'*skyflat*'),ntwi)
          if (ntwi ne 0) then new[twi].flavor = 'twiflat'
; remove crap exposures
          case night[inight] of
             '22jun11': keep = where($
               (strmatch(new.filename,'*test*') eq 0) and $ ; Feige34/crap
               (strmatch(new.filename,'*.000[1-6].*') eq 0) and $ ; focus
               (strmatch(new.filename,'*.0045.*') eq 0)) ; arc w/ wrong slit
             '23jun11': keep = where($
               (strmatch(new.filename,'*test*') eq 0) and $ ; Feige34/crap
               (strmatch(new.filename,'*.0112.*') eq 0) and $ ; wrong object!
               (strmatch(new.filename,'*.005[4-7].*') eq 0)) ; focus exposures
             '24jun11': keep = where($
               (strmatch(new.filename,'*test*') eq 0) and $
               (strmatch(new.filename,'*.0168.*') eq 0) and $ ; Cat's Eye saturated!
               (strmatch(new.filename,'*.0181.*') eq 0)) ; wrong slit
             '27jun11': keep = where($
               (strmatch(new.filename,'*test*') eq 0))
             '28jun11': keep = where($
               (strmatch(new.filename,'*test*') eq 0))
             '29jun11': keep = where($
               (strmatch(new.filename,'*test*') eq 0))
             else: message, 'Code me up!'
          endcase          
          new = new[keep]
          struct_print, new

; parameter defaults
          hdr = hdr[where(strcompress(hdr,/remove) ne '')]
          maxobj = where(strmatch(hdr,'*maxobj*'))
          hdr[maxobj] = 'maxobj 3'
          hdr = [hdr,'nozap 1','nohelio 1','noflex 0','novac 1','skytrace 0']

          yanny_write, planfile, ptr_new(new), hdr=hdr, /align
          write_calib_planfile, planfile, calibplanfile
          popd
       endfor
    endif 

; ##################################################
; reduce the calibration data
    if keyword_set(calib) then begin
       for inight = 0, nnight-1 do begin
          calibplanfile = 'plan_calib_'+night[inight]+'.par'
          pushd, datapath+night[inight]
          long_reduce, calibplanfile, /justcalib, calibclobber=clobber
          popd
       endfor
    endif
       
; ##################################################
; reduce and extract the standards, if any
    if keyword_set(standards) then begin
       for inight = 0, nnight-1 do begin
          calibplanfile = 'plan_calib_'+night[inight]+'.par'
          pushd, datapath+night[inight]
          long_reduce, calibplanfile, /juststd, sciclobber=clobber, chk=chk
          popd
       endfor
    endif

; ##################################################
; build the sensitivity function from all the nights 
;   push this to its own routine!
    if keyword_set(sensfunc) then begin
       delvarx, stdplan
       for inight = 0, nnight-1 do begin
          planfile = datapath+night[inight]+'/plan_'+night[inight]+'.par'
          stdplan1 = yanny_readone(planfile)
          std = where((strmatch(stdplan1.starname,'*...*') eq 0),nstd)
          if (nstd ne 0) then begin
             stdplan1 = stdplan1[std]
             stdplan1.filename = datapath+night[inight]+'/Science/std-'+$
               strtrim(stdplan1.filename,2)
             if (n_elements(stdplan) eq 0) then stdplan = stdplan1 else $
               stdplan = [stdplan,stdplan1]
          endif else splog, 'No standard stars observed on night '+ $
              night[inight]+'!'
       endfor
       if (n_elements(stdplan) eq 0) then $
          splog, 'No standard stars observed!' $
       else begin
          temp = (stdplan.maskname[sort(stdplan.maskname)])[$
              uniq(stdplan.maskname[sort(stdplan.maskname)])]
          sensfunc_maskname = (temp[sort(temp)])[uniq(temp[sort(temp)])]

          temp = (stdplan.wave[sort(stdplan.wave)])[$
              uniq(stdplan.wave[sort(stdplan.wave)])]
          sensfunc_wave = (temp[sort(temp)])[uniq(temp[sort(temp)])]

          temp = (stdplan.grating[sort(stdplan.grating)])[$
              uniq(stdplan.grating[sort(stdplan.grating)])]
          sensfunc_grating = (temp[sort(temp)])[uniq(temp[sort(temp)])]

          for imaskname = 0, n_elements(sensfunc_maskname)-1 do begin
            for iwave = 0, n_elements(sensfunc_wave)-1 do begin
              for igrating = 0, n_elements(sensfunc_grating)-1 do begin
    ;         keep = lindgen(n_elements(stdplan))
                keep = where($
                  (stdplan.maskname eq sensfunc_maskname[imaskname]) and $
                  (stdplan.wave eq sensfunc_wave[iwave]) and $
                  (stdplan.grating eq sensfunc_grating[igrating]) and $
                  (strmatch(stdplan.filename,'*22jun11.0052*',/fold) eq 0) and $
                  (strmatch(stdplan.filename,'*23jun11.0095*',/fold) eq 0) and $
                  (strmatch(stdplan.filename,'*27jun11.0244*',/fold) eq 0) and $
                  (strmatch(stdplan.filename,'*27jun11.025[2-3]*', $
                      /fold) eq 0) and $
                  (strmatch(stdplan.filename,'*23jun11.0109*',/fold) eq 0), $
                  nkeep)
                if (nkeep gt 0) then begin
                  stdplan_temp = stdplan[keep]
                  struct_print, stdplan_temp
                  nstd = n_elements(stdplan_temp)

                  sensfuncfile = "sensfunc_2011_"+strjoin(strsplit($
                      sensfunc_grating[igrating],'/',/extract),'.')+ $
                      "grating_"+strcompress(sensfunc_maskname[imaskname], $
                      /remove_all)+"slit_"+strcompress(sensfunc_wave[iwave], $
                      /remove_all)+"tilt.fits"

                  splog, 'Building '+sensfuncfile+' from '+ $
                      string(nstd,format='(I0)')+' standards'
              
                  aycamp_sensfunc, strtrim(stdplan_temp.filename,2), $
                      strtrim(stdplan_temp.starfile,2), nogrey=0, $
                      sensfuncfile=sensfuncfile
                endif
              endfor
            endfor
          endfor
       endelse
    endif 

; ##################################################
; reduce the objects by reading the extraction file
    if keyword_set(science) then begin
       for inight = 0, nnight-1 do begin
          pushd, datapath+night[inight]
          planfile = datapath+night[inight]+'/plan_'+night[inight]+'.par'
          handfile = datapath+night[inight]+'/handextract_'+night[inight]+'.par'
          if (file_test(handfile) eq 0) then begin
             splog, 'Handfile '+handfile+' not found!'
             continue
          endif
          hand = yanny_readone(handfile,/anonymous)

          delvarx, hand_fwhm, hand_x, hand_y
;          for ii = 2, 3 do begin
          for ii = 0, n_elements(hand)-1 do begin
             scifile = repstr(strtrim(hand[ii].filename,2),'.gz','')
             
             gdhand = where(hand[ii].hand_fwhm gt 0.0,ngdhand)
             if (ngdhand gt 0) then begin
                hand_fwhm = hand[ii].hand_fwhm[gdhand]
                hand_x = hand[ii].hand_x[gdhand]
                hand_y = hand[ii].hand_y[gdhand]
             endif else begin
                 hand_fwhm = 0
                 hand_x = 0
                 hand_y = 0
             endelse

             long_reduce, planfile, sciclobber=clobber, /justsci, $
               onlysci=scifile, hand_fwhm=hand_fwhm, hand_x=hand_x, $
               hand_y=hand_y, box_rad=hand[ii].box_rad, $
               maxobj=hand[ii].maxobj, nolocal=hand[ii].nolocal, $
               chk=chk, /skymask, /verbose, peak_smth=1.0
          endfor
          popd
       endfor 
    endif 

; ##################################################
; unpack each project in turn
    if unpack_something eq 1 then begin
       first = 0
       for inight = 0, nnight-1 do begin
          scifiles1 = file_search(datapath+night[inight]+$
              '/Science/sci-*.fits*', count=nscifiles)

          if (nscifiles gt 0) then begin
            if (first eq 0) then begin 
              scifiles = scifiles1 
              first = 1
            endif else  begin
              scifiles = [scifiles,scifiles1]
            endelse
          endif
       endfor
       allinfo = aycamp_forage(scifiles)
    endif

; -------------------------
; WISE    
    if keyword_set(unpack_wise) then begin
       outpath = projectpath+'wise/'
       if (file_test(outpath,/dir) eq 0) then spawn, 'mkdir -p '+outpath, /sh
       
       info = allinfo[where(strmatch(allinfo.object,'*wise*',/fold))]
       obj = strcompress(info.object,/remove)

       allgrp = spheregroup(15D*hms2dec(info.ra),hms2dec(info.dec),15D/3600.0)
       grp = allgrp[uniq(allgrp,sort(allgrp))]

       for ig = 0, n_elements(grp)-1 do begin
          these = where(grp[ig] eq allgrp,nthese)

          aperture = strcompress(info[these[0]].aperture,/remove)
          tilt = strcompress(info[these[0]].tiltpos,/remove)
          grating = strjoin(strsplit(strcompress(info[these[0]].disperse, $
              /remove),'/',/extract),'.')
          sensfuncfile = "sensfunc_2011_"+grating+"grating_"+aperture+"slit_"+ $
              tilt+"tilt.fits"

          coadd_outfile = outpath+obj[these[0]]+'.fits'
          aycamp_niceprint, info[these].file, obj[these]
          long_coadd, info[these].file, 1, outfil=coadd_outfile, /medscale, $
            box=0, check=0, /norej, /nosharp
; flux calibrate and write out the final 1D FITS and ASCII spectra
          outfile = repstr(coadd_outfile,'.fits','_f.fits')
          aycamp_fluxcalibrate, coadd_outfile, outfile=outfile, $
            sensfuncfile=sensfuncfile, /clobber, /writetxt
          aycamp_plotspec, outfile, /postscript, scale=1D16, objname=obj[these[0]]
       endfor
    endif

; -------------------------
; supernovae
    if keyword_set(unpack_sne) then begin
       outpath = projectpath+'sne/'
       if (file_test(outpath,/dir) eq 0) then spawn, 'mkdir -p '+outpath, /sh
       
       info = allinfo[where($
         strmatch(allinfo.object,'*M51*',/fold) or $
         strmatch(allinfo.object,'*psn*',/fold) or $
         strmatch(allinfo.object,'*2011dn*',/fold) or $
         strmatch(allinfo.object,'*ptf*',/fold))]

; give everything a nice name
       obj = strcompress(info.object,/remove)
       ww = where(strmatch(obj,'*ugc*',/fold) or strmatch(obj,'*2011dn*',/fold))
       obj[ww] = 'sn2011dn_'+repstr(info[ww].date,'-','_')
       obj[where(strmatch(obj,'*M51*',/fold))] = 'sn2011dh'
       obj[where(strmatch(obj,'*1631*',/fold))] = 'sn2011dw'
       obj[where(strmatch(obj,'*1612*',/fold))] = 'sn2011dv'
       obj[where(strmatch(obj,'*1759*',/fold))] = 'PSNJ1759'
       aycamp_niceprint, info.file, obj

       allgrp = spheregroup(15D*hms2dec(info.ra),hms2dec(info.dec),15D/3600.0)
       grp = allgrp[uniq(allgrp,sort(allgrp))]

       for ig = 0, n_elements(grp)-1 do begin
          these = where(grp[ig] eq allgrp,nthese)

          aperture = strcompress(info[these[0]].aperture,/remove)
          tilt = strcompress(info[these[0]].tiltpos,/remove)
          grating = strjoin(strsplit(strcompress(info[these[0]].disperse, $
              /remove),'/',/extract),'.')
          sensfuncfile = "sensfunc_2011_"+grating+"grating_"+aperture+"slit_"+ $
              tilt+"tilt.fits"

          coadd_outfile = outpath+obj[these[0]]+'.fits'
          aycamp_niceprint, info[these].file, obj[these]
; set the objid (need to figure these out individually using
; long_look) 
          objid = replicate(1,nthese)
          other = where($
            strmatch(info[these].file,'*22jun11.004[8-9]*',/fold) or $ ; =sn2011dn_2011_06_23
            strmatch(info[these].file,'*23jun11.010[6-7]*',/fold) or $ ; =sn2011dn_2011_06_24
            strmatch(info[these].file,'*27jun11.024[7-9]*',/fold))     ; =PTF11gjq
          if (other[0] ne -1) then objid[other] = 2
          long_coadd, info[these].file, objid, outfil=coadd_outfile, /medscale, $
            box=0, check=0, /norej, /nosharp
; flux calibrate and write out the final 1D FITS and ASCII spectra
          outfile = repstr(coadd_outfile,'.fits','_f.fits')
          aycamp_fluxcalibrate, coadd_outfile, outfile=outfile, $
            sensfuncfile=sensfuncfile, /clobber, /writetxt
          aycamp_plotspec, outfile, /postscript, scale=1D16, objname=obj[these[0]]
       endfor
    endif
    
; -------------------------
; NGC4559 rotation curve project
    if keyword_set(unpack_ngc4559) then begin
       outpath = projectpath+'ngc4559/'
       if (file_test(outpath,/dir) eq 0) then spawn, 'mkdir -p '+outpath, /sh

; one set of back-to-back observations of this object
       info = allinfo[where(strmatch(strcompress(allinfo.object,/remove),'*4559*'))]
       obj = 'NGC4559'

; we extracted five apertures for this object centered on H-alpha;
; coadd each spectrum/aperture individually
       nap = 5
       for jj = 0, nap-1 do begin
          aperture = strcompress(info[0].aperture,/remove)
          tilt = strcompress(info[0].tiltpos,/remove)
          grating = strjoin(strsplit(strcompress(info[0].disperse, $
              /remove),'/',/extract),'.')
          sensfuncfile = "sensfunc_2011_"+grating+"grating_"+aperture+"slit_"+ $
              tilt+"tilt.fits"

          coadd_outfile = outpath+obj+'_aper'+string(jj+1,format='(I0)')+'.fits'
          long_coadd, info.file, jj+1, outfil=coadd_outfile, /medscale, $
            box=0, check=0, /norej;, /nosharp
; flux calibrate and write out the final 1D FITS and ASCII spectra
          outfile = repstr(coadd_outfile,'.fits','_f.fits')
          aycamp_fluxcalibrate, coadd_outfile, outfile=outfile, $
            sensfuncfile=sensfuncfile, /clobber, /writetxt
          aycamp_plotspec, outfile, /postscript, scale=1D16, $
            objname=obj+'_aper'+string(jj+1,format='(I0)')
       endfor
    endif

; -------------------------
; Vy Aqr cataclysmic variable project
    if keyword_set(unpack_vyaqr) then begin
       outpath = projectpath+'vyaqr/'
       if (file_test(outpath,/dir) eq 0) then spawn, 'mkdir -p '+outpath, /sh

; one set of back-to-back observations of this object
       info = allinfo[where(strmatch(strcompress(allinfo.object,/remove),'*vyaqr*',/fold),nobj)]
       obj = 'VyAqr'
       
       for jj = 0, nobj-1 do begin
          aperture = strcompress(info[0].aperture,/remove)
          tilt = strcompress(info[0].tiltpos,/remove)
          grating = strjoin(strsplit(strcompress(info[0].disperse, $
              /remove),'/',/extract),'.')
          sensfuncfile = "sensfunc_2011_"+grating+"grating_"+aperture+"slit_"+ $
              tilt+"tilt.fits"

          coadd_outfile = outpath+obj+'_spec'+string(jj+1,format='(I2.2)')+'.fits'
          long_coadd, info.file, 1, outfil=coadd_outfile, /medscale, $
            box=0, check=0, /norej;, /nosharp
; flux calibrate and write out the final 1D FITS and ASCII spectra
          outfile = repstr(coadd_outfile,'.fits','_f.fits')
          aycamp_fluxcalibrate, coadd_outfile, outfile=outfile, $
            sensfuncfile=sensfuncfile, /clobber, /writetxt
          aycamp_plotspec, outfile, /postscript, scale=1D16, $
            objname=obj+'_spec'+string(jj+1,format='(I0)')
       endfor
    endif
    
; -------------------------
; lenticular bars vs no bars project
    if keyword_set(unpack_lenticulars) then begin
       outpath = projectpath+'lenticulars/'
       if (file_test(outpath,/dir) eq 0) then spawn, 'mkdir -p '+outpath, /sh

       info = allinfo[where($
         strmatch(allinfo.object,'*5273*',/fold) or $
         strmatch(allinfo.object,'*5195*',/fold) or $
         strmatch(allinfo.object,'*5353*',/fold) or $
         strmatch(allinfo.object,'*5710*',/fold) or $
         strmatch(allinfo.object,'*5354*',/fold) or $
         strmatch(allinfo.object,'*5864*',/fold) or $
         strmatch(allinfo.object,'*102*',/fold) or $
         strmatch(allinfo.object,'*5485*',/fold) or $
         strmatch(allinfo.object,'*6548*',/fold) or $
         strmatch(allinfo.object,'*6654*',/fold))]
       obj = strcompress(info.object,/remove)
       aycamp_niceprint, info.file, obj

       allgrp = spheregroup(15D*hms2dec(info.ra),hms2dec(info.dec),15D/3600.0)
       grp = allgrp[uniq(allgrp,sort(allgrp))]

       for ig = 0, n_elements(grp)-1 do begin
          these = where(grp[ig] eq allgrp,nthese)

          aperture = strcompress(info[these[0]].aperture,/remove)
          tilt = strcompress(info[these[0]].tiltpos,/remove)
          grating = strjoin(strsplit(strcompress(info[these[0]].disperse, $
              /remove),'/',/extract),'.')
          sensfuncfile = "sensfunc_2011_"+grating+"grating_"+aperture+"slit_"+ $
              tilt+"tilt.fits"

          coadd_outfile = outpath+obj[these[0]]+'.fits'
          aycamp_niceprint, info[these].file, obj[these]
          long_coadd, info[these].file, 1, outfil=coadd_outfile, /medscale, $
            box=0, check=0, /norej, /nosharp
; flux calibrate and write out the final 1D FITS and ASCII spectra
          outfile = repstr(coadd_outfile,'.fits','_f.fits')
          aycamp_fluxcalibrate, coadd_outfile, outfile=outfile, $
            sensfuncfile=sensfuncfile, /clobber, /writetxt
          aycamp_plotspec, outfile, /postscript, scale=1D16, objname=obj[these[0]]
       endfor
    endif

; -------------------------
; planetary nebulae
    if keyword_set(unpack_pne) then begin
       outpath = projectpath+'pne/'
       if (file_test(outpath,/dir) eq 0) then spawn, 'mkdir -p '+outpath, /sh

       info = allinfo[where($
         strmatch(allinfo.object,'*egg*',/fold) or $
         strmatch(allinfo.object,'*6058*',/fold) or $
         strmatch(allinfo.object,'*ring*',/fold) or $
         strmatch(allinfo.object,'*dumb*',/fold) or $
         strmatch(allinfo.object,'*omega*',/fold) or $
         strmatch(allinfo.object,'*cat*',/fold))]
       obj = strcompress(info.object,/remove)
;      aycamp_niceprint, info.file, obj, info.ra, info.dec

       allgrp = spheregroup(15D*hms2dec(info.ra),hms2dec(info.dec),30D/3600.0)
       grp = allgrp[uniq(allgrp,sort(allgrp))]

       for ig = 0, n_elements(grp)-1 do begin
          these = where(grp[ig] eq allgrp,nthese)

          aperture = strcompress(info[these[0]].aperture,/remove)
          tilt = strcompress(info[these[0]].tiltpos,/remove)
          grating = strjoin(strsplit(strcompress(info[these[0]].disperse, $
              /remove),'/',/extract),'.')
          sensfuncfile = "sensfunc_2011_"+grating+"grating_"+aperture+"slit_"+ $
              tilt+"tilt.fits"

          coadd_outfile = outpath+obj[these[0]]+'.fits'
          aycamp_niceprint, info[these].file, obj[these]
          long_coadd, info[these].file, 1, outfil=coadd_outfile, /medscale, $
            box=1, check=0, /norej, /nosharp
; flux calibrate and write out the final 1D FITS and ASCII spectra
          outfile = repstr(coadd_outfile,'.fits','_f.fits')
          aycamp_fluxcalibrate, coadd_outfile, outfile=outfile, $
            sensfuncfile=sensfuncfile, /clobber, /writetxt
          aycamp_plotspec, outfile, /postscript, scale=1D16, objname=obj[these[0]]
       endfor
    endif

; -------------------------
; rotation curve project with NGC3448 - incomplete reductions!!
    if keyword_set(unpack_rotationcurve) then begin
       outpath = projectpath+'rotationcurve/'
       if (file_test(outpath,/dir) eq 0) then spawn, 'mkdir -p '+outpath, /sh

       info = allinfo[where($
         strmatch(allinfo.object,'*3448*',/fold) or $
         strmatch(allinfo.object,'*4100*',/fold) or $
         strmatch(allinfo.object,'*4565*',/fold))]
       obj = ['NGC4565','NGC4565','NGC4565_dustlane','NGC3448',$
         'NGC3448','NGC4100','NGC4100']
;      aycamp_niceprint, info.file, info.object, obj

       allgrp = spheregroup(15D*hms2dec(info.ra),hms2dec(info.dec),10D/3600.0)
       grp = allgrp[uniq(allgrp,sort(allgrp))]

       for ig = 0, n_elements(grp)-1 do begin
          these = where(grp[ig] eq allgrp,nthese)

          aperture = strcompress(info[these[0]].aperture,/remove)
          tilt = strcompress(info[these[0]].tiltpos,/remove)
          grating = strjoin(strsplit(strcompress(info[these[0]].disperse, $
              /remove),'/',/extract),'.')
          sensfuncfile = "sensfunc_2011_"+grating+"grating_"+aperture+"slit_"+ $
              tilt+"tilt.fits"

          coadd_outfile = outpath+obj[these[0]]+'.fits'
          aycamp_niceprint, info[these].file, obj[these]
          long_coadd, info[these].file, 1, outfil=coadd_outfile, /medscale, $
            box=1, check=0, /norej, /nosharp
; flux calibrate and write out the final 1D FITS and ASCII spectra
          outfile = repstr(coadd_outfile,'.fits','_f.fits')
          aycamp_fluxcalibrate, coadd_outfile, outfile=outfile, $
            sensfuncfile=sensfuncfile, /clobber, /writetxt
          aycamp_plotspec, outfile, /postscript, scale=1D16, objname=obj[these[0]]
       endfor
    endif

return
end
