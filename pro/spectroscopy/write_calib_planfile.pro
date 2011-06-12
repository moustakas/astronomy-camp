pro write_calib_planfile, planfile, calibplanfile
; support routine for the reduce_bok_???? routines; write out a plan
; file for just the calibration data 
    
    plan = yanny_readone(planfile,hdr=hdr)
    calib = where(strtrim(plan.flavor,2) ne 'science')
    new = plan[calib]

    plotfile = where(strmatch(hdr,'*plotfile*'))
    hdr[plotfile] = 'plotfile '+repstr(calibplanfile,'.par','.ps')
    logfile = where(strmatch(hdr,'*logfile*'))
    hdr[logfile] = 'logfile '+repstr(calibplanfile,'.par','.log')

    splog, 'Writing '+calibplanfile
    yanny_write, calibplanfile, ptr_new(new), hdr=hdr, /align

return
end

