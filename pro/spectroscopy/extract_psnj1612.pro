pro psn_writetxt, spec1dfile
; write out the final spectrum
    
    txtfile = repstr(spec1dfile,'.fits','.txt')
    trim_wave = 5

    flux = mrdfits(spec1dfile,0,hdr,/silent)
    ferr = mrdfits(spec1dfile,1,/silent)
    wave = mrdfits(spec1dfile,2,/silent)
    npix = n_elements(wave)

    flux = 1E-17*flux[trim_wave:npix-trim_wave-1]
    ferr = 1E-17*ferr[trim_wave:npix-trim_wave-1]
    wave = wave[trim_wave:npix-trim_wave-1]
    npix = n_elements(wave)

    openw, lun, txtfile, /get_lun
    aycamp_niceprintf, lun, wave, flux, ferr
    free_lun, lun

return
end

pro extract_psnj1612
; jm11jun27ucsd - extract Psnj1612

    datapath = getenv('AYCAMP_DATA')+'2011/bok/'
    outpath = getenv('AYCAMP_DATA')+'2011/projects/psnj1612/'
    if (file_test(outpath) eq 0) then spawn, 'mkdir -p '+outpath, /sh
    if (file_test(outpath+'Fspec') eq 0) then spawn, 'mkdir -p '+outpath+'Fspec/', /sh
    
    sensfuncfile = datapath+'sensfunc_2011.fits'
    
    root = 'PSNJ1612+14'
    infiles = file_search(datapath+'29jun11/Science/sci-29jun11.035[7-9].fits.gz')
    objid = [1,1,1]
    color = ['red','blue','dark green']
    
    for ii = 0, n_elements(objid)-1 do begin
       suffix = '_'+strtrim(ii+1,2)
       fspecfile = outpath+'Fspec/'+root+suffix+'.fits'
       long_coadd, infiles[ii], objid[ii], wave=wave, flux=flux, $
         ivar=ivar, outfil=fspecfile, /medscale, box=0, check=0

       spec1dfile = outpath+root+suffix+'.fits'
       long_fluxcal, fspecfile, sensfunc=sensfuncfile, outfil=spec1dfile
       psn_writetxt, spec1dfile
    endfor

; make the coadd; leave the last exposure off    
    fspecfile = outpath+'Fspec/'+root+'_coadd.fits'
    long_coadd, infiles, objid, wave=wave, flux=flux, $
      ivar=ivar, outfil=fspecfile, /medscale, box=0, $
      check=0, iref=2

    coaddfile = outpath+root+'_coadd.fits'
    long_fluxcal, fspecfile, sensfunc=sensfuncfile, outfil=coaddfile
    psn_writetxt, coaddfile

; make a QAplot
    coaddfile = file_search(outpath+root+'_coadd.txt')
    allfile = file_search(outpath+root+'_?.txt',count=nall)

    psfile = outpath+'qa_'+root+'.ps'
    aycamp_plotconfig, 8, pos, psfile=psfile
; page 1
    djs_plot, [0], [0], xrange=[3650,6900], yrange=[0.5,12], $
      xsty=1, ysty=1, position=pos, xtitle='Wavelength (\AA)', $
      ytitle='Flux (10^{-16} '+flam_units()+')'
    legend, root, /right, /top, box=0
    readcol, coaddfile, wave, flux, ferr, /silent, format='F,F,F'
    djs_oplot, wave, 1D16*flux, ps=10
    
; page 2    
    djs_plot, [0], [0], xrange=[3650,6900], yrange=[0.5,12], $
      xsty=1, ysty=1, position=pos, xtitle='Wavelength (\AA)', $
      ytitle='Flux (10^{-16} '+flam_units()+')'
    legend, root, /right, /top, box=0
    for ii = 0, nall-1 do begin
       readcol, allfile[ii], wave, flux, ferr, /silent, format='F,F,F'
       djs_oplot, wave, 1D16*flux, ps=10, color=color[ii], thick=1
    endfor
    readcol, coaddfile, wave, flux, ferr, /silent, format='F,F,F'
    djs_oplot, wave, 1D16*flux, ps=10
    aycamp_plotconfig, psfile=psfile, /psclose, /pdf

    pdffile = file_basename(repstr(psfile,'.ps','.pdf'))
    
    pushd, outpath
    spawn, 'tar czvf psnj1612.tar.gz '+pdffile+' '+$
      strjoin(file_search(root+'*.txt'),' '), /sh
    popd
      
return
end
    
