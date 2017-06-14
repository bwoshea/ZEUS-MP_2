;hdfSim - a general plotting routine for visualising 
;         2-d data on a 2-d plot
; hdfsim.pro
; by Charlotte Dirk
; ********************************
pro go

infile = FindFile('hdfaa.*', Count=numfiles)

thisDevice = !D.NAME
Set_Plot, 'Z', /COPY

Device, Set_Resolution=[640,320], Z_Buffer=0
Erase

LoadCT, 13

jstart = 0
jend   = numfiles-1
jskip  = 1

for j=jstart,jend,jskip do begin

 filename = infile(j)
 FileID   = HDF_SD_START(filename, /READ)

 pi      = dataselect(FileID, 24)
 angle   = dataselect(FileID, 25)
 radius  = dataselect(FileID, 26)

 density = dataselect(FileID, 27)

 hdf_sd_end, FileID

 circle2d,radius,angle,pi,density,file='MHD Blast Wave: Density'

 image = tvrd()
 TVLCT, red, green, blue, /GET
 outfile = 'MHDblast' + String(j, Format='(I3.3)') + '.png'
 write_png, outfile, image, red, green, blue

endfor

end

;******************************************** 
;		circle using angles
;******************************************** 
pro circle2d,r,angle,pi,z,file=file;type=type
		; z is the value at each radius
		
z=transpose(z)
nr= size(r)
nr = nr[1];
nz = size(z)

nlevels    = 256
step       = (max(z)-min(z)) / nlevels
userlevels = min(z) + Indgen(nlevels)*step

polar_contour,z,angle,r, Levels=userlevels,/isotropic,/fill,color=110,title=file

return

end
;****************************************** dataselect
Function dataselect,FileID, nn
        sds= HDF_SD_SELECT(FileID, nn)
        HDF_SD_GETDATA, SDS, newData
        return,  newData
end
