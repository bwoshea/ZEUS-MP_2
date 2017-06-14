PRO go

hdr1 = strarr(3)
x1   = fltarr(256)
x2   = fltarr(256)
x3   = fltarr(256)
rho  = fltarr(256)
tgas = fltarr(256)
trad = fltarr(256)

junk = ""
time = 0.0e0

infile = FindFile('usraa*', Count=numfiles)

thisDevice = !D.NAME
Set_Plot, 'Z', /COPY

Device, Set_Resolution=[640,480], Z_Buffer=0
Erase

colors = FSC_COLOR(['black', 'yellow', 'green', 'red'], Indgen(4)+1)

plot, x1, tgas, xtitle='Distance (1e10 cm)', $
      title='Super-Critical Radiating Shock',charsize=1.0, $
      xrange=[0.0,10.0], xstyle=1, xticks=10, xminor=10, $
      yrange=[0.0,5000.0],ystyle=1, yticks=5, yminor=10, $
      position = [0.15, 0.15, 0.85, 0.95], $
      Color=colors[1],Background=colors[0], /nodata
xyouts, 0.05, 0.33, 'Tgas (K)', alignment=0.50, orientation=90, $
        Color=colors[3], charsize=1.5, /normal
xyouts, 0.05, 0.50, ' and ', alignment=0.50, orientation=90, $
        Color=colors[1], charsize=1.5, /normal
xyouts, 0.05, 0.67, 'Trad (K)', alignment=0.50, orientation=90, $
        Color=colors[2], charsize=1.5, /normal

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

  readf, 1, x, y, z, d, tg, tr

  x1(j)   = x
  rho(j)  = d
  tgas(j) = tg
  trad(j) = tr

 endfor

 tms_str = STRTRIM(time,2)

 oplot, x1, tgas, Color=colors[3]
 oplot, x1, trad, Color=colors[2]

 close, 1

endfor

image = tvrd()
TVLCT, red, green, blue, /GET
outfile = 'RSHOCK' + '.png'
write_png, outfile, image, red, green, blue

Device, Z_Buffer=1
SET_PLOT, thisDevice

end
