PRO go

hdr1 = strarr(3)
x1   = fltarr(256)
x2   = fltarr(256)
x3   = fltarr(256)
rho  = fltarr(256)
egas = fltarr(256)
vel  = fltarr(256)

junk = ""
time = 0.0e0

infile = FindFile('usraa*.001', Count=numfiles)

thisDevice = !D.NAME
Set_Plot, 'Z', /COPY

Device, Set_Resolution=[640,480], Z_Buffer=0
Erase

colors = FSC_COLOR(['black', 'yellow', 'green', 'red'], Indgen(4)+1)

plot, x1, rho, xtitle='Radius (cm)', $
      title='Pressureless Gravitational Collapse',charsize=1.0, $
      xrange=[0.0,1.0], xstyle=1, xticks=10, xminor=10, $
      yrange=[0.0,1000.0],ystyle=1, yticks=10, yminor=5, $
      position = [0.15, 0.15, 0.85, 0.95], $
      Color=colors[1],Background=colors[0], /nodata
xyouts, 0.05, 0.50, ' Rho (g/cc) ', alignment=0.50, orientation=90, $
        Color=colors[1], charsize=1.5, /normal

for i=0,numfiles-1 do begin

 openr, 1, infile(i)

 x  = 0.0
 y  = 0.0
 z  = 0.0
 d  = 0.0
 tg = 0.0
 tr = 0.0

 readf, 1, hdr1
 reads, hdr1(0), junk, time, format='(A34,x,E15.8)'

 for j=0,255 do begin

  readf, 1, x, y, z, d, eg, v

  x1(j)   = x
  rho(j)  = d

 endfor

 tms_str = STRTRIM(time,2)

 oplot, x1, rho, Color=colors[3]

 close, 1

endfor

image = tvrd()
TVLCT, red, green, blue, /GET
outfile = 'GCOL' + '.png'
write_png, outfile, image, red, green, blue

Device, Z_Buffer=1
SET_PLOT, thisDevice

end
