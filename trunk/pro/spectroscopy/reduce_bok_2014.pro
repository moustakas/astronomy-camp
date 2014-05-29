;+
; NAME:
;   REDUCE_BOK_2014
; PURPOSE:
;   Process all the Bok B&C spectroscopy obtained during ATC2014. 
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
;   J. Moustakas, 2013 Jun 15, UCSD
;-

; long_setflex, 'Science/sci-23jun11.0092.fits.gz', 'bok_sky_bcs_400.sav'

pro reduce_bok_2014, night, preproc=preproc, plan=plan, calib=calib, $
  standards=standards, science=science, sensfunc=sensfunc, chk=chk, $
  unpack_kepler=unpack_kepler, unpack_ngc5273=unpack_ngc5273, $
  unpack_sne=unpack_sne, unpack_eclipse=unpack_eclipse, $
  unpack_transit=unpack_transit, unpack_standards=unpack_standards, $
  clobber=clobber

    ; !!!!!!!UPDATE - need to update this to reflect the projects done during
    ;                 ATC2014.
    unpack_something = (keyword_set(unpack_kepler) or $
        keyword_set(unpack_ngc5273) or keyword_set(unpack_sne) or $
        keyword_set(unpack_eclipse) or keyword_set(unpack_transit) or $
        keyword_set(unpack_standards)) ? 1 : 0

    datapath = getenv('AYCAMP_DATA')+'2014/bok/'
    projectpath = datapath+'projects/'
    if (file_test(projectpath,/dir) eq 0) then $
      spawn, 'mkdir -p '+projectpath, /sh
    
    badpixfile = getenv('AYCAMP_DIR')+'/data/bok_badpix.dat'
    if (file_test(badpixfile) eq 0) then message, $
      'Bad pixel file '+badpixfile+' not found'
;    sensfuncfile = datapath+'sensfunc_2013.fits'

    ;!!!!!!UPDATE - update this to reflect the nights we observed during ATC2014
    if (n_elements(night) eq 0) then night = ['25jun13','26jun13','27jun13']
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

                ; Overwrite any incorrect header info here.
                ;if strmatch(allfiles[iobj],'*28jun11.0310.*',/fold) then $
                ;    sxaddpar, hdr, 'APERTURE', '4.5'

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
             ;!!!!!!!UPDATE - Update to reflect errors in data during ATC2014
             '25jun13': keep = where($
               (strmatch(new.filename,'*test*') eq 0) and $
               (strmatch(new.filename,'*.0079.*') eq 0) and $
               (strmatch(new.filename,'*.008[0-2].*') eq 0))
             '26jun13' : keep = where($
               (strmatch(new.filename,'*test*') eq 0))
             '27jun13' : keep = where($
               (strmatch(new.filename,'*test*') eq 0) and $
               (strmatch(new.filename,'*.002[5-6].*') eq 0) and $
               (strmatch(new.filename,'*.004[6-7].*') eq 0) and $
               (strmatch(new.filename,'*.0078.*') eq 0))
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
                  ; Remove any bad standards here.
                  ;!!!!!!!UPDATE - remove any ATC2014 bad standards
                  (strmatch(stdplan.filename,'*22jun11.0052*',/fold) eq 0), $
                  nkeep)
                if (nkeep gt 0) then begin
                  stdplan_temp = stdplan[keep]
                  struct_print, stdplan_temp
                  nstd = n_elements(stdplan_temp)

                  sensfuncfile = "sensfunc_2014_"+strjoin(strsplit($
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

;!!!!!!!!!!!!UPDATE - need to add/remove/modify projects below to reflect the
;                     projects done during ATC2014

; -------------------------
; Kepler   
    if keyword_set(unpack_kepler) then begin
       outpath = projectpath+'kepler/'
       if (file_test(outpath,/dir) eq 0) then spawn, 'mkdir -p '+outpath, /sh
       
       info = allinfo[where(strmatch(allinfo.object,'*KID*',/fold))]
       obj = strcompress(info.object,/remove)

       allgrp = spheregroup(15D*hms2dec(info.ra),hms2dec(info.dec),15D/3600.0)
       grp = allgrp[uniq(allgrp,sort(allgrp))]

       for ig = 0, n_elements(grp)-1 do begin
          these = where(grp[ig] eq allgrp,nthese)

          aperture = strcompress(info[these[0]].aperture,/remove)
          tilt = strcompress(info[these[0]].tiltpos,/remove)
          grating = strjoin(strsplit(strcompress(info[these[0]].disperse, $
              /remove),'/',/extract),'.')
          sensfuncfile = "sensfunc_2014_"+grating+"grating_"+aperture+"slit_"+ $
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
; NGC5273   
    if keyword_set(unpack_ngc5273) then begin
       outpath = projectpath+'ngc5273/'
       if (file_test(outpath,/dir) eq 0) then spawn, 'mkdir -p '+outpath, /sh
       
       info = allinfo[where(strmatch(allinfo.object,'*5273*',/fold))]
       obj = strcompress(info.object,/remove)

       allgrp = spheregroup(15D*hms2dec(info.ra),hms2dec(info.dec),15D/3600.0)
       grp = allgrp[uniq(allgrp,sort(allgrp))]

       for ig = 0, n_elements(grp)-1 do begin
          these = where(grp[ig] eq allgrp,nthese)

          aperture = strcompress(info[these[0]].aperture,/remove)
          tilt = strcompress(info[these[0]].tiltpos,/remove)
          grating = strjoin(strsplit(strcompress(info[these[0]].disperse, $
              /remove),'/',/extract),'.')
          sensfuncfile = "sensfunc_2014_"+grating+"grating_"+aperture+"slit_"+ $
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
; Supernovae   
    if keyword_set(unpack_sne) then begin
       outpath = projectpath+'sne/'
       if (file_test(outpath,/dir) eq 0) then spawn, 'mkdir -p '+outpath, /sh
       
       info = allinfo[where(strmatch(allinfo.object,'*PSN*',/fold) or $
           strmatch(allinfo.object,'*UGC*',/fold))]
       obj = strcompress(info.object,/remove)

       allgrp = spheregroup(15D*hms2dec(info.ra),hms2dec(info.dec),15D/3600.0)
       grp = allgrp[uniq(allgrp,sort(allgrp))]

       for ig = 0, n_elements(grp)-1 do begin
          these = where(grp[ig] eq allgrp,nthese)

          aperture = strcompress(info[these[0]].aperture,/remove)
          tilt = strcompress(info[these[0]].tiltpos,/remove)
          grating = strjoin(strsplit(strcompress(info[these[0]].disperse, $
              /remove),'/',/extract),'.')
          sensfuncfile = "sensfunc_2014_"+grating+"grating_"+aperture+"slit_"+ $
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
; Eclipsing binary   
    if keyword_set(unpack_eclipse) then begin
       outpath = projectpath+'eclipse/'
       if (file_test(outpath,/dir) eq 0) then spawn, 'mkdir -p '+outpath, /sh
       
       info = allinfo[where(strmatch(allinfo.object,'*TCas0*',/fold))]
       obj = strcompress(info.object,/remove)

       allgrp = spheregroup(15D*hms2dec(info.ra),hms2dec(info.dec),15D/3600.0)
       grp = allgrp[uniq(allgrp,sort(allgrp))]

       for ig = 0, n_elements(grp)-1 do begin
          these = where(grp[ig] eq allgrp,nthese)

          aperture = strcompress(info[these[0]].aperture,/remove)
          tilt = strcompress(info[these[0]].tiltpos,/remove)
          grating = strjoin(strsplit(strcompress(info[these[0]].disperse, $
              /remove),'/',/extract),'.')
          sensfuncfile = "sensfunc_2014_"+grating+"grating_"+aperture+"slit_"+ $
              tilt+"tilt.fits"

          coadd_outfile = outpath+strjoin(strsplit(obj[these[0]],'/',/extract),$
              '-')+'.fits'
          aycamp_niceprint, info[these].file, obj[these]
          long_coadd, info[these].file, 2, outfil=coadd_outfile, /medscale, $
            box=0, check=0, /norej, /nosharp
; flux calibrate and write out the final 1D FITS and ASCII spectra
          outfile = repstr(coadd_outfile,'.fits','_f.fits')
          aycamp_fluxcalibrate, coadd_outfile, outfile=outfile, $
            sensfuncfile=sensfuncfile, /clobber, /writetxt
          aycamp_plotspec, outfile, /postscript, scale=1D16, objname=obj[these[0]]
       endfor
    endif

; -------------------------
; Transit   
    if keyword_set(unpack_transit) then begin
       outpath = projectpath+'transit/'
       if (file_test(outpath,/dir) eq 0) then spawn, 'mkdir -p '+outpath, /sh
       
       info = allinfo[where(strmatch(allinfo.object,'*CoRoT*',/fold))]
       obj = strcompress(info.object,/remove)

       allgrp = spheregroup(15D*hms2dec(info.ra),hms2dec(info.dec),15D/3600.0)
       grp = allgrp[uniq(allgrp,sort(allgrp))]

       for ig = 0, n_elements(grp)-1 do begin
          these = where(grp[ig] eq allgrp,nthese)

          aperture = strcompress(info[these[0]].aperture,/remove)
          tilt = strcompress(info[these[0]].tiltpos,/remove)
          grating = strjoin(strsplit(strcompress(info[these[0]].disperse, $
              /remove),'/',/extract),'.')
          sensfuncfile = "sensfunc_2014_"+grating+"grating_"+aperture+"slit_"+ $
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
; Standards  
    if keyword_set(unpack_standards) then begin
       outpath = projectpath+'standards/'
       if (file_test(outpath,/dir) eq 0) then spawn, 'mkdir -p '+outpath, /sh
       
       info = allinfo[where(strmatch(allinfo.object,'*Vega*',/fold) or $
           strmatch(allinfo.object,'*Cyg*',/fold))]
       obj = strcompress(info.object,/remove)

       allgrp = spheregroup(15D*hms2dec(info.ra),hms2dec(info.dec),15D/3600.0)
       grp = allgrp[uniq(allgrp,sort(allgrp))]

       for ig = 0, n_elements(grp)-1 do begin
          these = where(grp[ig] eq allgrp,nthese)

          aperture = strcompress(info[these[0]].aperture,/remove)
          tilt = strcompress(info[these[0]].tiltpos,/remove)
          grating = strjoin(strsplit(strcompress(info[these[0]].disperse, $
              /remove),'/',/extract),'.')
          sensfuncfile = "sensfunc_2014_"+grating+"grating_"+aperture+"slit_"+ $
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

return
end
