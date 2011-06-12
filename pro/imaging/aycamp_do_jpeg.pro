pro aycamp_do_jpeg, bluefile, greenfile, redfile, $
  scales=scales, jpegfile=jpegfile
; generate the color mosaics

    if (n_elements(scales) eq 0) then scales = [1.0,1.0,1.0]
    
    blue = mrdfits(bluefile)
    red = mrdfits(redfile)
    green = mrdfits(greenfile)
    
    splog, 'Writing '+jpegfile
    nw_rgb_make, red, green, blue, name=jpegfile, scales=scales, $
      nonlinearity=3.0, rebinfactor=2, quality=90

return
end       

