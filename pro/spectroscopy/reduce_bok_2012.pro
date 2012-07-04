;+
; NAME:
;   REDUCE_BOK_2012
; PURPOSE:
;   Process all the Bok B&C spectroscopy obtained during ATC2012. 
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
;   J. Moustakas, 2011 Jun 15, UCSD (reduce_bok_2012.pro, starting point)
;   B. Swift, 2012 Jun 20-28, Steward Obs.
;-

; long_setflex, 'Science/sci-23jun11.0092.fits.gz', 'bok_sky_bcs_400.sav'

pro reduce_bok_2012, night, preproc=preproc, plan=plan, calib=calib, $
        standards=standards, sensfunc=sensfunc, science=science, $
        clobber=clobber, chk=chk, camp=camp, $
        unpack_blackhole=unpack_blackhole, unpack_test=unpack_test
;perhaps change the unpack keywords to "unpack=unpack, _EXTRA=extra"?
; where it would be called with "... /unpack, /sne", for example.
; this way, the code doesn't have to be tied to an individual year's
; projects. and then break out the "unpack" code into another
; procedure anyway...
  
  unpack_something = (keyword_set(unpack_blackhole) or keyword_set(unpack_test))
  
  ;check if any of the step keywords are defined?  if not, complain, exit.

;this needs to read the rawdata/ subdirectories.  if 'night' not
;defined, then read that listing and present as a menu selection,
;including "ALL OF THE ABOVE" as a valid option.  maybe also give the
;option to select multiple nights in some easy fashion.

  if not keyword_set(camp) then camp='2012'
  
  if camp eq '2012a' then $
     night0 = ['21jun12','22jun12'] $ ;'23jun12'];,'26jun12','27jun12']
  else if camp eq '2012' then $
     night0 = ['21jun12','22jun12','23jun12','26jun12','01jul12'] 

  if (n_elements(night) eq 0) then night = night0
  nnight = n_elements(night)
  
  ;this line is ALSO tied to a particular camp 
  datapath = getenv('AYCAMP_DATA')+camp+'/bok/'
  projectpath = datapath+'projects/'
  if (file_test(projectpath,/dir) eq 0) then $
     spawn, 'mkdir -p '+projectpath, /sh
  
  ;
  ;make the campers build a bad pixel mask?  could be interesting.
  ;
  badpixfile = getenv('AYCAMP_DIR')+'/data/bok_badpix.dat'
  if (file_test(badpixfile) eq 0) then message, $
     'Bad pixel file '+badpixfile+' not found'

    ;AGAIN: NEED TO BREAK THIS OUT TO BE AGNOSTIC TO CAMP DATASET
    ;ALSO, NEED IT TO BE SPECIFIC TO SETUPS!
    sensfuncfile = datapath+'sensfunc_2012.fits'
    
    if ((night eq '21jun12') or (night eq '22jun12') or (night eq '23jun12') $
        (array_equal(night,night0)) or $
        (array_equal(night,['21jun12','22jun12'])) or $
        (array_equal(night,['21jun12','23jun12'])) or $
        (array_equal(night,['22jun12','23jun12'])) or $
        (array_equal(night,['22jun12','21jun12'])) or $
        (array_equal(night,['23jun12','21jun12'])) or $
        (array_equal(night,['23jun12','22jun12']))) then $
        sensfuncfile = datapath+$
            'sensfunc_2012_400line4889_6.2tilt_2.5slit.fits' $
    else if (night eq '26jun12') then $
        sensfuncfile = datapath+$
            'sensfunc_2012_1200line7847_25.0tilt_2.5slit.fits'$
    else if (night eq '01jul12') then $
        sensfuncfile = datapath+$
            'sensfunc_2012_400line4889_6.58tilt_2.5slit.fits' $
    else begin 
        print, 'ERROR - CODE THIS UP: NIGHT not included in '+$
            'REDUCE_BOK_2012.PRO near line 90-ish.'
        return
    endelse


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
          
          allfiles = file_search(datapath+'rawdata/'+night[inight]+'/*.fits*',$
              count=nspec)
          if (nspec eq 0) then message, 'No files found in '+datapath+'rawdata/'

          for iobj = 0, nspec-1 do begin
             outfile = repstr('Raw/'+file_basename(allfiles[iobj]),'.gz','')
             if file_test(outfile+'.gz') and (keyword_set(clobber) eq 0) then $
                 begin
                splog, 'Output file '+outfile+' exists; use /CLOBBER'
             endif else begin
                image = mrdfits(allfiles[iobj],0,hdr,/fscale,/silent)

                sxaddpar, hdr, 'INSTRUME', 'bcspeclamps', ' instrument name'
                sxaddpar, hdr, 'DISPERSE', '400/4889', ' disperser'
                sxaddpar, hdr, 'APERTURE', '2.5', ' slit width'

                if strmatch(allfiles[iobj],'*21jun12*',/fold) then sxaddpar, $
                    hdr, 'TILTPOS', '6.2'
                if strmatch(allfiles[iobj],'*21jun12*',/fold) then sxaddpar, $
                    hdr, 'OBSERVER', 'Camp'

                if strmatch(allfiles[iobj],'*22jun12*',/fold) then sxaddpar, $
                    hdr, 'INSFILTE', 'UV-36+Schott-8612' 
;we didn't think we had a filter in, so we put 'none'.  turns out we did.  
;Betsey probably thinks we're huge idiots.

                if strmatch(allfiles[iobj],'*26jun12*',/fold) then sxaddpar, $
                    hdr, 'DISPERSE', '1200/7847', ' disperser' 
                if strmatch(allfiles[iobj],'*26jun12*',/fold) then sxaddpar, $
                    hdr, 'INSFOCUS', '006' 
                if strmatch(allfiles[iobj],'*26jun12*',/fold) then sxaddpar, $
                    hdr, 'TILTPOS', '20.0' 
                if strmatch(allfiles[iobj],'*26jun12_00[0-2]*',/fold) then $
                    sxaddpar, hdr, 'TILTPOS', '25.0' 

                type = sxpar(hdr,'imagetyp')
                if (strlowcase(strtrim(type,2)) eq 'object') then begin
                   dims = size(image,/dim)
                   badpixmask = image*0.0
                   for ipix = 0, n_elements(x1)-1 do $
                     badpixmask[(x1[ipix]-1)>0:(x2[ipix]-1)<(dims[0]-1),$
                     (y1[ipix]-1)>0:(y2[ipix]-1)<(dims[1]-1)] = 1
                   image = djs_maskinterp(image,badpixmask,iaxis=0,/const)
                endif
                exp = sxpar(hdr,'exptime')
                if (strlowcase(strtrim(type,2)) eq 'object') and $
                    (exp gt 300.0) then begin
                   splog, 'Identifying cosmic rays '
                   aycamp_la_cosmic, image, gain=-1.0, readn=0.0, $
                     outlist=newim, sigclip=5.0
                   image = newim
                endif

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
                
                aycamp_mwrfits, image[*,0:data_arr[3]-1], outfile, hdr, $
                    clobber=clobber
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
;THIS NEEDS TO BE BROKEN OUT OF HERE!!!
          case night[inight] of
             ;; CAMP ATC2012-A
             '21jun12': keep = where( $
                (strmatch(new.filename,'*test*') eq 0) and $  
                (strmatch(new.filename,'*.000[1-9].*') eq 0) and $ ; focus
                (strmatch(new.filename,'*.001[0-3].*') eq 0) and $ ; focus, 13=saturated continuum
                (strmatch(new.filename,'*.0014.*') eq 0) and $ ; continuum lamp, won't use it
                (strmatch(new.filename,'*.001[5-6].*') eq 0) and $ ; standard, not at parallactic
                (strmatch(new.filename,'*.003[1-2].*') eq 0) and $     ; wierd bias?
                (strmatch(new.filename,'*.0019.*') eq 0) and $     ; wierd bias?
                (strmatch(new.filename,'*.0040.*') eq 0) and $     ; object, but lamps were on
                (strmatch(new.filename,'*.0044.*') eq 0) and $     ; standard, not at parallactic
                (strmatch(new.filename,'*.004[6-9].*') eq 0) and $ ; bad flats
                (strmatch(new.filename,'*.005[0-2].*') eq 0))     ; bad flats 

             '22jun12': keep = where( $
                (strmatch(new.filename,'*test*') eq 0) and $ 
                (strmatch(new.filename,'*.0009.*') eq 0) and $ ; continuum, won't use it
                (strmatch(new.filename,'*.0010.*') eq 0)) ; continuum, won't use it
             
             '23jun12': keep = where( $
                (strmatch(new.filename,'*test*') eq 0))

             '26jun12': keep = where( $
                (strmatch(new.filename,'*test*') eq 0))

             '27jun12': keep = where( $
                (strmatch(new.filename,'*test*') eq 0)) 
             
             '01jul12': keep = where( $
                (strmatch(new.filename,'*test*') eq 0))

             '03jul12': keep = where( $
                (strmatch(new.filename,'*test*') eq 0))

             '04jul12': keep = where( $
                (strmatch(new.filename,'*test*') eq 0))

             '05jul12': keep = where( $
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
          endif else splog, 'No standard stars observed on night '+$
              night[inight]+'!'
       endfor
       if (n_elements(stdplan) eq 0) then splog, $
           'No standard stars observed!' $
        else begin
          keep = where( $
            ; List any bad standards here.
            (strmatch(stdplan.filename,'*23jun12_0118*',/fold) eq 0) and $
            (strmatch(stdplan.filename,'*23jun12_0121*',/fold) eq 0) and $
            (strmatch(stdplan.filename,'*01jul12.0061*',/fold) eq 0))

          stdplan = stdplan[keep]
          struct_print, stdplan
          nstd = n_elements(stdplan)
          splog, 'Building '+sensfuncfile+' from '+string(nstd,format='(I0)')+$
              ' standards'
          
          aycamp_sensfunc, strtrim(stdplan.filename,2), $
              strtrim(stdplan.starfile,2), nogrey=0, sensfuncfile=sensfuncfile
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
          for ii = 0, n_elements(hand)-1 do begin
             scifile = repstr(strtrim(hand[ii].filename,2),'.gz','')
             
             gdhand = where(hand[ii].hand_fwhm gt 0.0,ngdhand)
             if (ngdhand gt 0) then begin
                hand_fwhm = hand[ii].hand_fwhm[gdhand]
                hand_x    = hand[ii].hand_x[gdhand]
                hand_y    = hand[ii].hand_y[gdhand]
             endif

             long_reduce, planfile, sciclobber=clobber, /justsci, $
               onlysci=scifile, hand_fwhm=hand_fwhm, hand_x=hand_x, $
               hand_y=hand_y, box_rad=hand[ii].box_rad, $
               maxobj=hand[ii].maxobj, nolocal=hand[ii].nolocal, $
               chk=chk, /skymask, /verbose, peak_smth=1.0
          endfor
          popd
       endfor 
    endif 

;########################################################
; Unpack the various projects.

    if unpack_something eq 1 then begin
        for inight = 0, nnight-1 do begin
            scifiles1 = file_search(datapath+night[inight]+$
                '/Science/sci-*.fits*', count=niscifiles)
            if niscifiles gt 0 then begin
                if (inight eq 0) then scifiles = scifiles1 $
                else scifiles = [scifiles,scifiles1]
            endif
        endfor
        allinfo = aycamp_forage(scifiles)
    endif

; ------------------------------------------------
; The test data from the second Advanced Camp.

    if keyword_set(unpack_test) then begin
       outpath = projectpath+'test/'
       if (file_test(outpath,/dir) eq 0) then spawn, 'mkdir -p '+outpath, /sh

       info = allinfo[where($
         strmatch(allinfo.object,'*RandomStar*',/fold))]
       obj = strcompress(info.object,/remove)
;      aycamp_niceprint, info.file, obj, info.ra, info.dec

       allgrp = spheregroup(15D*hms2dec(info.ra),hms2dec(info.dec),30D/3600.0)
       grp = allgrp[uniq(allgrp,sort(allgrp))]

       for ig = 0, n_elements(grp)-1 do begin
          these = where(grp[ig] eq allgrp,nthese)
          coadd_outfile = outpath+obj[these[0]]+'.fits'
          aycamp_niceprint, info[these].file, obj[these]
          long_coadd, info[these].file, 1, outfil=coadd_outfile, /medscale, $
            box=1, check=0, /norej, /nosharp
; flux calibrate and write out the final 1D FITS and ASCII spectra
          outfile = repstr(coadd_outfile,'.fits','_f.fits')
          aycamp_fluxcalibrate, coadd_outfile, outfile=outfile, $
            sensfuncfile=sensfuncfile, /clobber, /writetxt
          aycamp_plotspec, outfile, /postscript, scale=1D16, $
            objname=obj[these[0]]
       endfor
    endif

; ----------------------------------------------------------------
; The binary black holes

    if keyword_set(unpack_blackhole) then begin
       outpath = projectpath+'BlackHole/'
       if (file_test(outpath,/dir) eq 0) then spawn, 'mkdir -p '+outpath, /sh

       info = allinfo[where($
         strmatch(allinfo.object,'*SS433*',/fold) or $
         strmatch(allinfo.object,'*Cygnus*',/fold))]
       obj = strcompress(info.object,/remove)
;      aycamp_niceprint, info.file, obj, info.ra, info.dec

       allgrp = spheregroup(15D*hms2dec(info.ra),hms2dec(info.dec),30D/3600.0)
       grp = allgrp[uniq(allgrp,sort(allgrp))]

       for ig = 0, n_elements(grp)-1 do begin
          these = where(grp[ig] eq allgrp,nthese)
          coadd_outfile = outpath+obj[these[0]]+'.fits'
          aycamp_niceprint, info[these].file, obj[these]
          long_coadd, info[these].file, 1, outfil=coadd_outfile, /medscale, $
            box=1, check=0, /norej, /nosharp
; flux calibrate and write out the final 1D FITS and ASCII spectra
          outfile = repstr(coadd_outfile,'.fits','_f.fits')
          aycamp_fluxcalibrate, coadd_outfile, outfile=outfile, $
            sensfuncfile=sensfuncfile, /clobber, /writetxt
          aycamp_plotspec, outfile, /postscript, scale=1D16, $
            objname=obj[these[0]]
       endfor
    endif

return
end
