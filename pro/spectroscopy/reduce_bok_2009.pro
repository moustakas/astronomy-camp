;+
; NAME:
;   REDUCE_BOK_2009
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
;   J. Moustakas, 2010 May 25, UCSD
;-

pro reduce_bok_2009, datapath, night=night, fixbadpix=fixbadpix, plan=plan, $
  reduce=reduce, sensfunc=sensfunc, coadd=coadd, fluxcal=fluxcal, $
  clobber=clobber

    if (n_elements(datapath) eq 0) then datapath = getenv('HOME')+'/'+$
      'home/research/data/aycamp/2009/90inch/'
    splog, 'Data path '+datapath

    badpixfile = getenv('AYCAMP_DIR')+'/data/bok_badpix.dat'
    if (file_test(badpixfile) eq 0) then message, $
      'Bad pixel file '+badpixfile+' not found'
    
    if (n_elements(night) eq 0) then night = $
      ['20jun09','22jun09','23jun09']
;     ['20jun09','21jun09','22jun09','23jun09']
;   night = ['21jun09']
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

       sensfuncfile = 'sensfunc_'+night[inight]+'.fits'
       planfile = 'plan_'+night[inight]+'.par'

; ##################################################
; repair bad pixels, clean up headers, and move the spectra to the Raw
; subdirectory for further processing
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
          long_plan, '*.fits.gz', 'Raw/', $
            planfile=planfile
; make an edit script here
          old = yanny_readone(planfile,hdr=hdr)
          new = aycamp_find_calspec(old,radius=radius)
          maxobj = where(strmatch(hdr,'*maxobj*'))
          hdr[maxobj] = 'maxobj 1'
          struct_print, new
          yanny_write, planfile, ptr_new(new), hdr=hdr, /align
       endif 

; ##################################################
; reduce everything
       if keyword_set(reduce) then begin
          long_reduce, planfile, clobber=clobber, /verbose, $
            /nozap, chk=chk
       endif

; ##################################################
; build the sensitivity function for this night
       if keyword_set(sensfunc) then begin
          thisplan = yanny_readone(planfile)
          thisplan.filename = 'Science/std-'+strtrim(thisplan.filename,2)
          
          std = where((strmatch(thisplan.starname,'*...*') eq 0),nstd)
          splog, 'Found '+string(nstd,format='(I0)')+' standards'
          if (nstd ne 0) then begin
             stdfiles = strtrim(thisplan[std].filename,2)
             std_names = strtrim(thisplan[std].starfile,2)
             sens = aycamp_long_sensfunc(stdfiles,sensfuncfile,$
               std_name=std_names,/msk_balm)
          endif
       endif

; ##################################################
; coadd multiple exposures of the same object and flux-calibrate 
       trim = 10
       if keyword_set(fluxcal) then begin
          infiles = file_search('Science/sci-*.fits*',count=nspec)
          info = iforage(infiles)
          ra = 15D*im_hms2dec(info.ra)
          dec = im_hms2dec(info.dec)
          obj = strcompress(info.object,/remove)
          allgrp = spheregroup(ra,dec,10/60.0)
          grp = allgrp[uniq(allgrp,sort(allgrp))]
          for ig = 0, n_elements(grp)-1 do begin
             these = where(grp[ig] eq allgrp,nthese)
             outfile = 'Fspec/'+obj[these[0]]+'.fits'
             long_coadd, infiles[these], 1, outfil=outfile
          endfor
; flux-calibrate, trim crap pixels from each end and convert to
; standard FITS format
          infiles = file_search('Fspec/*.fits*',count=nspec)
          outfiles = 'spec1d/'+file_basename(infiles)
          for iobj = 0, nspec-1 do begin
             long_fluxcal, infiles[iobj], sensfunc=sensfuncfile, $
               outfil=outfiles[iobj]
; trim, etc.
             flux = mrdfits(outfiles[iobj],0,hdr,/silent)
             ferr = mrdfits(outfiles[iobj],1,/silent)
             wave = mrdfits(outfiles[iobj],2,/silent)
             npix = n_elements(wave)
             flux = 1E-17*flux[trim:npix-trim-1]
             ferr = 1E-17*ferr[trim:npix-trim-1]
             wave = wave[trim:npix-trim-1]
             npix = n_elements(wave)
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
             djs_plot, wave, 1D17*flux, xsty=3, ysty=3, ps=10, $
               position=pos, xtitle='Wavelength (\AA)', $
               ytitle='Flux (10^{-17} '+flam_units()+')'
             legend, [obj], /right, /top, box=0
          endfor
          im_plotconfig, psfile=psfile, /psclose, /gzip
       endif
       popd
    endfor ; close night

return
end
