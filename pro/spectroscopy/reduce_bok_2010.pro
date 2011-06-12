;+
; NAME:
;   REDUCE_BOK_2010
;
; PURPOSE:
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Jun 15, UCSD
;-

pro reduce_bok_2010, night, fixbadpix=fixbadpix, plan=plan, calib=calib, $
  standards=standards, science=science, sensfunc=sensfunc, $
  unpack_projects=unpack_projects, fluxcal=fluxcal, clobber=clobber, $
  chk=chk, qaplot=qaplot

    datapath = getenv('AYCAMP_DATA')+'2010/bok/'
    splog, 'Data path '+datapath

    badpixfile = getenv('AYCAMP_DIR')+'/data/bok_badpix.dat'
    if (file_test(badpixfile) eq 0) then message, $
      'Bad pixel file '+badpixfile+' not found'
    sensfuncfile = datapath+'sensfunc_2010.fits'

    trim_sedge = 5 ; trim region
    trim_wave = trim_sedge*2
    
    if (n_elements(night) eq 0) then night = $
      ['16jun10','17jun10','18jun10','20jun10','21jun10']
    nnight = n_elements(night)

    for inight = 0, nnight-1 do begin
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

;      sensfuncfile = 'sensfunc_'+night[inight]+'.fits'
       planfile = 'plan_'+night[inight]+'.par'
       calibplanfile = 'plan_calib_'+night[inight]+'.par'

; ##################################################
; reading the data in the "rawdata" directory, repair bad pixels,
; remove cosmic rays, clean up headers, and move the spectra to the
; "Raw" subdirectory for further processing 
       if keyword_set(fixbadpix) then begin
          splog, 'Reading '+badpixfile
          readcol, badpixfile, x1, x2, y1, y2, comment='#', $
            format='L,L,L,L', /silent
          allfiles = file_search('rawdata/*.fits*',count=nspec)
          if (nspec eq 0) then begin
             splog, 'No files found in '+'rawdata/'
             continue
          endif
          for iobj = 0, nspec-1 do begin
             outfile = repstr('Raw/'+file_basename(allfiles[iobj]),'.gz','')
             if file_test(outfile+'.gz') and (keyword_set(clobber) eq 0) then begin
                splog, 'Output file '+outfile+' exists; use /CLOBBER'
             endif else begin
                image = mrdfits(allfiles[iobj],0,hdr,/fscale,/silent)
                sxaddpar, hdr, 'INSTRUME', 'bcspeclamps', ' instrument name'
                sxaddpar, hdr, 'DISPERSE', '400/4889', ' disperser'
                sxaddpar, hdr, 'APERTURE', '2.5', ' slit width'
; fix headers and bad pixels
                if strmatch(allfiles[iobj],'*18Jun10_0050*',/fold) then $
                  sxaddpar, hdr, 'OBJECT', 'VeilNebula No.3'
                type = sxpar(hdr,'imagetyp')
                if (strlowcase(strtrim(type,2)) eq 'object') then begin
                   dims = size(image,/dim)
                   badpixmask = image*0.0
                   for ipix = 0, n_elements(x1)-1 do $
                     badpixmask[(x1[ipix]-1)>0:(x2[ipix]-1)<(dims[0]-1),$
                     (y1[ipix]-1)>0:(y2[ipix]-1)<(dims[1]-1)] = 1
                   image = djs_maskinterp(image,badpixmask,iaxis=0,/const)
                endif
                im_mwrfits, image, outfile, hdr, /clobber

;; trim the top and bottom NTRIM pixels; need to update the header 
;                ntrim = 5
;
;                datasec = strcompress(sxpar(hdr,'DATASEC'),/rem)
;                biassec = strcompress(sxpar(hdr,'BIASSEC'),/rem)
;                data_arr = long(strsplit(datasec,'[*:*,*:*]',/extract))
;                bias_arr = long(strsplit(biassec,'[*:*,*:*]',/extract))
;
;                sxaddpar, hdr, 'DATASEC', '['+strtrim(data_arr[0],2)+':'+$
;                  strtrim(data_arr[1],2)+','+strtrim(data_arr[2],2)+':'+$
;                  strtrim(data_arr[3]-2*ntrim,2)+']'
;                sxaddpar, hdr, 'TRIMSEC', '['+strtrim(data_arr[0],2)+':'+$
;                  strtrim(data_arr[1],2)+','+strtrim(data_arr[2],2)+':'+$
;                  strtrim(data_arr[3]-2*ntrim,2)+']'
;                sxaddpar, hdr, 'BIASSEC', '['+strtrim(bias_arr[0],2)+':'+$
;                  strtrim(bias_arr[1],2)+','+strtrim(bias_arr[2],2)+':'+$
;                  strtrim(bias_arr[3]-2*ntrim,2)+']'
;                im_mwrfits, image[*,ntrim:data_arr[3]-ntrim-1], outfile, hdr, /clobber
             endelse
          endfor
       endif 
       
; ##################################################
; make the plan files
       if keyword_set(plan) then begin
          allfiles = file_search('Raw/*.fits*',count=nspec)
          if (nspec eq 0) then begin
             splog, 'No files found in '+'Raw/'
             continue
          endif
          long_plan, '*.fits.gz', 'Raw/', planfile=planfile
; make an edit script here
          old = yanny_readone(planfile,hdr=hdr)
          new = aycamp_find_calspec(old,radius=radius)
          twi = where(strmatch(new.target,'*twi*'),ntwi)
          if (ntwi ne 0) then new[twi].flavor = 'twiflat'
; remove crap exposures
          case night[inight] of
             '16jun10': keep = lindgen(n_elements(new))
             '17jun10': keep = where($
               (strmatch(new.filename,'*0016*') eq 0) and $ ; Feige34/crap
               (strmatch(new.filename,'*0017*') eq 0) and $ ; Feige34/crap
               (strmatch(new.filename,'*0032*') eq 0) and $ ; CygX1
               (strmatch(new.filename,'*0033*') eq 0) and $ ; CygX1
               (strmatch(new.filename,'*0034*') eq 0) and $ ; CygX1
               (strmatch(new.filename,'*0023*') eq 0) and $
               (strmatch(new.filename,'*0024*') eq 0) and $
               (strmatch(new.filename,'*0025*') eq 0) and $
               (strmatch(new.filename,'*0026*') eq 0))
             '18jun10': keep = where($
               (strmatch(new.filename,'*0036*') eq 0) and $ ; SN2010ei 
               (strmatch(new.filename,'*0037*') eq 0))      ; SN2010ei
             '20jun10': begin
                keep = where($
                  (strmatch(new.filename,'*0048*') eq 0) and $
                  (strmatch(new.filename,'*0064*') eq 0))
             end
             '21jun10': keep = lindgen(n_elements(new))
             else: message, 'Fix me'
          endcase          
          new = new[keep]
          struct_print, new

; parameter defaults
          hdr = hdr[where(strcompress(hdr,/remove) ne '')]
          maxobj = where(strmatch(hdr,'*maxobj*'))
          hdr[maxobj] = 'maxobj 1'
          hdr = [hdr,'nozap 1','nohelio 1','noflex 1','novac 1','skytrace 0']

          yanny_write, planfile, ptr_new(new), hdr=hdr, /align
          write_calib_planfile, planfile, calibplanfile
       endif 

; ##################################################
; reduce the calibration data
       if keyword_set(calib) then begin
          long_reduce, calibplanfile, /justcalib, calibclobber=clobber
       endif

; ##################################################
; reduce and extract the standards, if any
       if keyword_set(standards) then begin
          long_reduce, calibplanfile, /juststd, sciclobber=clobber
       endif

; ##################################################
; reduce the objects
       if keyword_set(science) then begin
          long_reduce, planfile, sciclobber=clobber, /justsci, chk=chk, /nolocal
       endif
    endfor 
       
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
          endif else splog, 'No standard stars observed on night '+night[inight]+'!'
       endfor
       if (n_elements(stdplan) ne 0) then begin
          struct_print, stdplan
          keep = where($
            (strmatch(stdplan.filename,'*21jun10.0034*',/fold) eq 0) and $
            (strmatch(stdplan.filename,'*18jun10_0033*',/fold) eq 0) and $
            (strmatch(stdplan.filename,'*18jun10_0034*',/fold) eq 0) and $
            (strmatch(stdplan.filename,'*20jun10.0078*',/fold) eq 0) and $
            (strmatch(stdplan.filename,'*17jun10_0018*',/fold) eq 0))
          stdplan = stdplan[keep]
          nstd = n_elements(stdplan)
          splog, 'Building '+sensfuncfile+' from '+string(nstd,format='(I0)')+' standards'
          stdfiles = strtrim(stdplan.filename,2)
          std_names = strtrim(stdplan.starfile,2)

; do the fit, masking the bluest and reddest pixels
          ncol = 1200 & nmask = 10
          inmask = intarr(ncol)+1
          inmask[0:nmask-1] = 0
          inmask[ncol-nmask-1:ncol-1] = 0
          sens = long_sensfunc(stdfiles,sensfuncfile,sensfit=sensfit,$
            std_name=std_names,wave=wave,flux=flux,nogrey=0,inmask=inmask,$
;            /msk_balm)
       endif else splog, 'No standard stars observed!'
    endif

; ##################################################
;   push this to its own routine!
    if keyword_set(fluxcal) then begin
       for inight = 0, nnight-1 do begin
          pushd, datapath+night[inight]
          infiles = file_search('Science/sci-*.fits*',count=nspec)
          info = aycamp_forage(infiles)
          ra = 15D*im_hms2dec(info.ra)
          dec = im_hms2dec(info.dec)
          obj = strcompress(info.object,/remove)
          allgrp = spheregroup(ra,dec,15/3600.0)
          grp = allgrp[uniq(allgrp,sort(allgrp))]

pick up here - coadd observations across
unpack projects!
          
          
          for ig = 0, n_elements(grp)-1 do begin
             these = where(grp[ig] eq allgrp,nthese)
             outfile = 'Fspec/'+obj[these[0]]+'.fits'
             skyfile = 'Fspec/'+obj[these[0]]+'.sky.fits'
             niceprint, infiles[these]
             long_coadd, infiles[these], 1, wave=wave, flux=flux, $
               ivar=ivar, outfil=outfile, skyfil=skyfile, /box
          endfor

; flux-calibrate, trim crap pixels from each end and convert to
; standard FITS format
          infiles = file_search('Fspec/*.fits*',count=nspec)
          nosky = where(strmatch(infiles,'*sky*',/fold) eq 0,nspec)
          infiles = infiles[nosky]
          outfiles = 'spec1d/'+file_basename(infiles)
          for iobj = 0, nspec-1 do begin
             long_fluxcal, infiles[iobj], sensfunc=sensfuncfile, $
               outfil=outfiles[iobj]
; trim, etc.
             flux = mrdfits(outfiles[iobj],0,hdr,/silent)
             ferr = mrdfits(outfiles[iobj],1,/silent)
             wave = mrdfits(outfiles[iobj],2,/silent)
             npix = n_elements(wave)
             flux = 1E-17*flux[trim_wave:npix-trim_wave-1]
             ferr = 1E-17*ferr[trim_wave:npix-trim_wave-1]
             wave = wave[trim_wave:npix-trim_wave-1]
             npix = n_elements(wave)
;; interpolate over pixels affected by strong sky lines
;             mask = ((wave gt 5545.0) and (wave lt 5595.0)) or $
;               ((wave gt 6280.0) and (wave lt 6307.0)) or $
;               ((wave gt 6345.0) and (wave lt 6368.0))
;             flux = djs_maskinterp(flux,mask,wave,/const)
; rebin linearly in wavelength
             dwave = ceil(100D*(max(wave)-min(wave))/(npix-1))/100D
             newwave = dindgen((max(wave)-min(wave))/dwave+1)*dwave+min(wave)
             newflux = rebin_spectrum(flux,wave,newwave)
             newvar = rebin_spectrum(ferr^2,wave,newwave)
             newferr = sqrt(newvar*(newvar gt 0.0)) ; enforce positivity
; build the final header; also grab the redshift from NED
             sxaddpar, hdr, 'CRVAL1', min(wave), ' wavelength at CRPIX1'
             sxaddpar, hdr, 'CRPIX1', 1D, ' reference pixel number'
             sxaddpar, hdr, 'CD1_1', dwave, ' dispersion [Angstrom/pixel]'
             sxaddpar, hdr, 'CDELT1', dwave, ' dispersion [Angstrom/pixel]'
             sxaddpar, hdr, 'CTYPE1', 'LINEAR', ' projection type'
; write out                   
             openw, lun, repstr(outfiles[iobj],'.fits','.txt'), /get_lun
             niceprintf, lun, wave, flux, ferr
             free_lun, lun

             mwrfits, float(newflux), outfiles[iobj], hdr, /create
             mwrfits, float(newferr), outfiles[iobj], hdr
             spawn, 'gzip -f '+outfiles[iobj], /sh
          endfor

; build a QAplot for all the spectra from this night
          psfile = 'spec1d/qa_'+night[inight]+'.ps'
          im_plotconfig, 8, pos, psfile=psfile
          for iobj = 0, nspec-1 do begin
             flux = mrdfits(outfiles[iobj]+'.gz',0,hdr,/silent)
             ferr = mrdfits(outfiles[iobj]+'.gz',1,/silent)
             wave = make_wave(hdr)
             obj = sxpar(hdr,'object')
             djs_plot, wave, 1D16*flux, xsty=3, ysty=3, ps=10, $
               position=pos, xtitle='Wavelength (\AA)', $
               ytitle='Flux (10^{-16} '+flam_units()+')'
             legend, [obj], /right, /top, box=0
          endfor
          im_plotconfig, psfile=psfile, /psclose, /gzip
       endfor
    endif


    
    
stop    
    
    
; ##################################################
; coadd multiple exposures of the same object and flux-calibrate 
;   push this to its own routine!
    if keyword_set(fluxcal) then begin
       for inight = 0, nnight-1 do begin
          planfile = datapath+night[inight]+'/plan_'+night[inight]+'.par'
          plan1 = yanny_readone(planfile)
          sci = where(strtrim(plan1.flavor,2) eq 'science')
          plan1 = plan1[sci]
          plan1.filename = datapath+night[inight]+'/Science/sci-'+strtrim(plan1.filename,2)
          if (inight eq 0) then plan = plan1 else $
            plan = [plan,plan1]
       endfor

       infiles = strtrim(plan.filename,2)
       nspec = n_elements(infiles)
       for ii = 0, nspec-1 do begin
          suffix = strmid(repstr(file_basename(infiles[ii]),'.fits.gz',''),11,/reverse)
          obj = strcompress(sxpar(headfits(infiles[ii]),'object'),/remove)
          outfile = datapath+'spec1d/'+obj+'_'+suffix+'.fits'
          long_coadd, infiles[ii], 1, outfil=outfile
       endfor

; now read the files back in, flux-calibrate, trim crap pixels from
; each end and convert to standard FITS format 
       trim = 10
       infiles = file_search(datapath+'spec1d/*.fits',count=nspec)
       outfiles = infiles
       for iobj = 0, nspec-1 do begin
          long_fluxcal, infiles[iobj], sensfunc=sensfuncfile, $
            outfil=outfiles[iobj]
; trim etc
          flux = mrdfits(outfiles[iobj],0,hdr,/silent)
          ferr = mrdfits(outfiles[iobj],1,/silent)
          wave = mrdfits(outfiles[iobj],2,/silent)
          npix = n_elements(wave)
          flux = 1E-17*flux[trim:npix-trim-1]
          ferr = 1E-17*ferr[trim:npix-trim-1]
          wave = wave[trim:npix-trim-1]
          npix = n_elements(wave)
; interpolate over pixels affected by strong sky lines
          mask = ((wave gt 5547.0) and (wave lt 5580.0)) or $
            ((wave gt 6280.0) and (wave lt 6307.0)) or $
            ((wave gt 6345.0) and (wave lt 6368.0))
          flux = djs_maskinterp(flux,mask,wave,/const)
; rebin linearly in wavelength
          dwave = ceil(100D*(max(wave)-min(wave))/(npix-1))/100D
          newwave = dindgen((max(wave)-min(wave))/dwave+1)*dwave+min(wave)
          newflux = rebin_spectrum(flux,wave,newwave)
          newvar = rebin_spectrum(ferr^2,wave,newwave)
          newferr = sqrt(newvar*(newvar gt 0.0)) ; enforce positivity
; build the final header; also grab the redshift from NED
          sxaddpar, hdr, 'CRVAL1', min(wave), ' wavelength at CRPIX1'
          sxaddpar, hdr, 'CRPIX1', 1D, ' reference pixel number'
          sxaddpar, hdr, 'CD1_1', dwave, ' dispersion [Angstrom/pixel]'
          sxaddpar, hdr, 'CDELT1', dwave, ' dispersion [Angstrom/pixel]'
          sxaddpar, hdr, 'CTYPE1', 'LINEAR', ' projection type'
          
          mwrfits, float(newflux), outfiles[iobj], hdr, /create
          mwrfits, float(newferr), outfiles[iobj], hdr
          spawn, 'gzip -f '+outfiles[iobj], /sh
       endfor
    endif

    if keyword_set(qaplot) then begin
       outfiles = file_search(datapath+'spec1d/*.fits.gz',count=nspec)
; make a QAplot of everything       
       psfile = datapath+'spec1d/qa_all.ps'
       im_plotconfig, 8, pos, psfile=psfile
       for iobj = 0, nspec-1 do begin
          flux = mrdfits(outfiles[iobj],0,hdr,/silent)
          ferr = mrdfits(outfiles[iobj],1,/silent)
          wave = make_wave(hdr)
          obj = sxpar(hdr,'object')

          power = ceil(abs(alog10(median(flux))))
          scale = 10.0^power
          ytitle = 'Flux (10^{-'+string(power,format='(I0)')+'} '+flam_units()+')'

          stats = im_stats(flux,sigrej=20.0)
          if strmatch(obj,'*4593*',/fold) then $
            yrange = [stats.min,1.2*stats.maxrej] else $
              yrange = [stats.min,1.1*stats.max]

          djs_plot, wave, scale*flux, xsty=3, ysty=3, ps=10, $
            position=pos, xtitle='Wavelength (\AA)', $
            ytitle=ytitle, yrange=scale*yrange
          legend, [obj], /right, /top, box=0
       endfor
       im_plotconfig, psfile=psfile, /psclose, /gzip
    endif     

;; ##################################################
;; unpack all the projects into dedicated directories
;    if keyword_set(unpack_projects) then begin
;       for inight = 0, nnight-1 do begin
;          planfile = datapath+night[inight]+'/plan_'+night[inight]+'.par'
;          plan1 = yanny_readone(planfile)
;          sci = where(strtrim(plan1.flavor,2) eq 'science')
;          plan1 = plan1[sci]
;          plan1.filename = datapath+night[inight]+'/Science/sci-'+strtrim(plan1.filename,2)
;          if (inight eq 0) then plan = plan1 else $
;            plan = [plan,plan1]
;       endfor
;
;; planets
;       project = 'planets'
;       outdir = datapath+'projects/'+project+'/'
;       if (file_test(outdir,/dir) eq 0) then $
;         spawn, 'mkdir -p '+outdir, /sh
;
;       these = where($
;         strmatch(plan.target,'*G223*',/fold) or $
;         strmatch(plan.target,'*G223*',/fold) or $
;         strmatch(plan.target,'*abell72*',/fold),nthese)
;       if (nthese ne 0) then begin
;          kp09m_2010_unpack, struct[these], outdir=outdir, $
;            indir=datapath+night[inight]+'/Final/', $
;            night=night[inight]
;       endif
;    endif 
       
;; ##################################################
;; unpack all the projects into dedicated directories
;    if keyword_set(unpack_projects) do begin
;; -------------------------
;; planets
;       project = 'planets'
;       outdir = datapath+'projects/'+project+'/'
;       if (file_test(outdir,/dir) eq 0) then $
;         spawn, 'mkdir -p '+outdir, /sh
;       these = where($
;         strmatch(struct.obj,'*abell43*',/fold) or $
;         strmatch(struct.obj,'*ic4593*',/fold) or $
;         strmatch(struct.obj,'*abell39*',/fold) or $
;         strmatch(struct.obj,'*abell72*',/fold),nthese)
;       if (nthese ne 0) then begin
;          kp09m_2010_unpack, struct[these], outdir=outdir, $
;            indir=datapath+night[inight]+'/Final/', $
;            night=night[inight]
;       endif
;       
;       trim = 10
;       if keyword_set(fluxcal) then begin
;          infiles = file_search('Science/sci-*.fits*',count=nspec)
;          keep = where(strmatch(infiles,'*17Jun10_0022*') eq 0,nspec)
;          infiles = infiles[keep]
;
;          info = iforage(infiles)
;          ra = 15D*im_hms2dec(info.ra)
;          dec = im_hms2dec(info.dec)
;          obj = strcompress(info.object,/remove)
;          allgrp = spheregroup(ra,dec,15/3600.0)
;          grp = allgrp[uniq(allgrp,sort(allgrp))]
;          for ig = 0, n_elements(grp)-1 do begin
;             these = where(grp[ig] eq allgrp,nthese)
;             outfile = 'Fspec/'+obj[these[0]]+'.fits'
;             skyfile = 'Fspec/'+obj[these[0]]+'.sky.fits'
;             niceprint, infiles[these]
;             long_coadd, infiles[these], 1, wave=wave, $
;               flux=flux, ivar=ivar, outfil=outfile, $
;               skyfil=skyfile;, /box
;          endfor
;
;; flux-calibrate, trim crap pixels from each end and convert to
;; standard FITS format
;          infiles = file_search('Fspec/*.fits*',count=nspec)
;          nosky = where(strmatch(infiles,'*sky*',/fold) eq 0,nspec)
;          infiles = infiles[nosky]
;          outfiles = 'spec1d/'+file_basename(infiles)
;          for iobj = 0, nspec-1 do begin
;             long_fluxcal, infiles[iobj], sensfunc=sensfuncfile, $
;               outfil=outfiles[iobj]
;; trim, etc.
;             flux = mrdfits(outfiles[iobj],0,hdr,/silent)
;             ferr = mrdfits(outfiles[iobj],1,/silent)
;             wave = mrdfits(outfiles[iobj],2,/silent)
;             npix = n_elements(wave)
;             flux = 1E-17*flux[trim:npix-trim-1]
;             ferr = 1E-17*ferr[trim:npix-trim-1]
;             wave = wave[trim:npix-trim-1]
;             npix = n_elements(wave)
;; interpolate over pixels affected by strong sky lines
;             mask = ((wave gt 5547.0) and (wave lt 5580.0)) or $
;               ((wave gt 6280.0) and (wave lt 6307.0)) or $
;               ((wave gt 6345.0) and (wave lt 6368.0))
;             flux = djs_maskinterp(flux,mask,wave,/const)
;; rebin linearly in wavelength
;             dwave = ceil(100D*(max(wave)-min(wave))/(npix-1))/100D
;             newwave = dindgen((max(wave)-min(wave))/dwave+1)*dwave+min(wave)
;             newflux = rebin_spectrum(flux,wave,newwave)
;             newvar = rebin_spectrum(ferr^2,wave,newwave)
;             newferr = sqrt(newvar*(newvar gt 0.0)) ; enforce positivity
;; build the final header; also grab the redshift from NED
;             sxaddpar, hdr, 'CRVAL1', min(wave), ' wavelength at CRPIX1'
;             sxaddpar, hdr, 'CRPIX1', 1D, ' reference pixel number'
;             sxaddpar, hdr, 'CD1_1', dwave, ' dispersion [Angstrom/pixel]'
;             sxaddpar, hdr, 'CDELT1', dwave, ' dispersion [Angstrom/pixel]'
;             sxaddpar, hdr, 'CTYPE1', 'LINEAR', ' projection type'
;; write out                   
;;            openw, lun, repstr(outfiles,'.fits','.txt'), /get_lun
;;            niceprintf, lun, wave, flux, ferr
;;            free_lun, lun
;
;             mwrfits, float(newflux), outfiles[iobj], hdr, /create
;             mwrfits, float(newferr), outfiles[iobj], hdr
;             spawn, 'gzip -f '+outfiles[iobj], /sh
;          endfor
;; build a QAplot for all the spectra from this night
;          psfile = 'spec1d/qa_'+night[inight]+'.ps'
;          im_plotconfig, 8, pos, psfile=psfile
;          for iobj = 0, nspec-1 do begin
;             flux = mrdfits(outfiles[iobj]+'.gz',0,hdr,/silent)
;             ferr = mrdfits(outfiles[iobj]+'.gz',1,/silent)
;             wave = make_wave(hdr)
;             obj = sxpar(hdr,'object')
;             djs_plot, wave, 1D16*flux, xsty=3, ysty=3, ps=10, $
;               position=pos, xtitle='Wavelength (\AA)', $
;               ytitle='Flux (10^{-16} '+flam_units()+')'
;             legend, [obj], /right, /top, box=0
;          endfor
;          im_plotconfig, psfile=psfile, /psclose, /gzip
;       endif 
;       popd
;    endif 

return
end

