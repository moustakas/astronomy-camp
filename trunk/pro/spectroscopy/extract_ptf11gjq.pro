pro ptf_writetxt, spec1dfile
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

pro extract_ptf11gjq
; jm11jun27ucsd - extract PTF11gjq

    datapath = getenv('AYCAMP_DATA')+'2011/bok/'
    sensfuncfile = datapath+'sensfunc_2011.fits'
    
    root = 'PTF11gjq'
    infiles = file_search(datapath+'27jun11/Science/sci-27jun11.024[5-9].fits.gz')
    objid = [2,1,2,2,2]
    color = ['red','orange','blue','dark green','grey']
    
    for ii = 0, n_elements(objid)-1 do begin
       suffix = '_'+strtrim(ii+1,2)
       fspecfile = datapath+'27jun11/Fspec/'+root+suffix+'.fits'
       long_coadd, infiles[ii], objid[ii], wave=wave, flux=flux, $
         ivar=ivar, outfil=fspecfile, /medscale, box=0, check=0

       spec1dfile = datapath+'27jun11/spec1d/'+root+suffix+'.fits'
       long_fluxcal, fspecfile, sensfunc=sensfuncfile, outfil=spec1dfile
       ptf_writetxt, spec1dfile
    endfor

; make the coadd; leave the last exposure off    
    fspecfile = datapath+'27jun11/Fspec/'+root+'_coadd.fits'
    long_coadd, infiles[0:3], objid[0:3], wave=wave, flux=flux, $
      ivar=ivar, outfil=fspecfile, /medscale, box=0, $
      check=0, iref=1
    
    spec1dfile = datapath+'27jun11/spec1d/'+root+'_coadd.fits'
    long_fluxcal, fspecfile, sensfunc=sensfuncfile, outfil=spec1dfile
    ptf_writetxt, spec1dfile

; make a QAplot
    coaddfile = file_search(datapath+'27jun11/spec1d/'+root+'_coadd.txt')
    allfile = file_search(datapath+'27jun11/spec1d/'+root+'_?.txt',count=nall)

    psfile = datapath+'27jun11/spec1d/qa_'+root+'.ps'
    aycamp_plotconfig, 8, pos, psfile=psfile
; page 1
    djs_plot, [0], [0], xrange=[3650,6900], yrange=[0.05,2], $
      xsty=1, ysty=1, position=pos, xtitle='Wavelength (\AA)', $
      ytitle='Flux (10^{-16} '+flam_units()+')'
    legend, root, /left, /top, box=0
    readcol, coaddfile, wave, flux, ferr, /silent, format='F,F,F'
    djs_oplot, wave, 1D16*flux, ps=10
    
; page 2    
    djs_plot, [0], [0], xrange=[3650,6900], yrange=[0.05,2], $
      xsty=1, ysty=1, position=pos, xtitle='Wavelength (\AA)', $
      ytitle='Flux (10^{-16} '+flam_units()+')'
    legend, root, /left, /top, box=0
    for ii = 0, nall-1 do begin
       readcol, allfile[ii], wave, flux, ferr, /silent, format='F,F,F'
       djs_oplot, wave, 1D16*flux, ps=10, color=color[ii], thick=1
    endfor
    readcol, coaddfile, wave, flux, ferr, /silent, format='F,F,F'
    djs_oplot, wave, 1D16*flux, ps=10
    aycamp_plotconfig, psfile=psfile, /psclose, /pdf

    pushd, datapath+'27jun11/spec1d'
    spawn, 'tar czvf '+datapath+'27jun11/ptf11gjq.tar.gz '+$
      file_basename(repstr(psfile,'.ps','.pdf'))+' '+$
      strjoin(file_search(root+'*.txt'),' '), /sh
    popd
      
    
    
stop    

return
end
    
