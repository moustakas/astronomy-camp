function aycamp_readspec, file, wave=wave
; jm11jun29ucsd - read an 1D extracted spectrum
    
    if (n_elements(file) eq 0) then begin
       doc_library, 'aycamp_plotspec'
       return, -1
    endif

    if (file_test(file) eq 0) then begin
       splog, 'Spectral file '+file+' not found!'
       return, -1
    endif

    flux = mrdfits(file,0,hdr)
    wave = aycamp_make_wave(hdr)

return, flux
end
