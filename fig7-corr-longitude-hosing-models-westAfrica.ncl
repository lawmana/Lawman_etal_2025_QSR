;---EOF analysis of Selected regions

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
  dir = "/glade/u/home/xianwu/ncl/UTAustin/ncl/paleo/herinch/paper/data/"
  diro = "./"

  slev = 0.90
;  percentage = True
  percentage = False
  
  var = "Rainfall"
;  var = "TAS"

  ;regiontype = "land"

  region = "Northern_Western_Africa"
  region1 = "North_Tropical_Atlantic_SST"
  region2 = "North_Atlantic_SST"

  if (region .eq. "Northern_Western_Africa") then
  	latmin = 5.
  	latmax = 15.
  	lonmin = -20.
  	lonmax = 0.
  end if
  
  if (region1 .eq. "North_Tropical_Atlantic_SST") then
  	latmin1 = 12.
  	latmax1 = 22.
    lonmin1 = -80.
  	lonmax1 = -40.
  end if

  if (region2 .eq. "North_Atlantic_SST") then
  	latmin2 = 40.
  	latmax2 = 65.
    lonmin2 = -60.
  	lonmax2 = 0.
  end if

print("Reading data files...")

;---Rainfall   

  f = addfile( dir + "Regridded_precip_18models.nc", "r")
  rh = f->hose
  rc = f->control

  f = addfile( dir + "Regridded_tas_18models.nc", "r")  
  th = f->hose
  tc = f->control

  if(percentage)then
    do i=0,9
      do j=0,dt(0)-1
        do k=0,dt(1)-1   
        if( .not. (ismissing(rc(i,j,k))))then
          if(rc(i,j,k) .eq. 0)then
            rc(i,j,k) = rc@_FillValue
          end if
        end if
        end do
      end do
    end do
    rh = rh/rc
    
    do i=0,9
      do j=0,dt(0)-1
        do k=0,dt(1)-1   
        if( .not. (ismissing(tc(i,j,k))))then
          if(tc(i,j,k) .eq. 0)then
            tc(i,j,k) = tc@_FillValue
          end if
        end if
        end do
      end do
    end do
    th = th/tc
  end if  
    
    x1 = wgt_areaave_Wrap( th(model|:,{lat|latmin1:latmax1}, {lon|lonmin1:lonmax1}) , 1., 1., 0)   
  
    x2 = wgt_areaave_Wrap( th(model|:,{lat|latmin2:latmax2}, {lon|lonmin2:lonmax2}) , 1., 1., 0)   
  
  if( var .eq. "Rainfall")then
    xh = rh
  else
    xh = th
  end if
  
;---
  print ("Correlation analysis")
  
  if(region .eq. "East_Africa")then
    xs = wgt_areaave_Wrap( xh(model|:,{lat|latmins:latmaxs}, {lon|lonmins:lonmaxs}) , 1., 1., 0)   
    xn = wgt_areaave_Wrap( xh(model|:,{lat|latminn:latmaxn}, {lon|lonminn:lonmaxn}) , 1., 1., 0)   
    x = xs
    x = xs-xn 
  else
    x = wgt_areaave_Wrap( xh(model|:,{lat|latmin:latmax}, {lon|lonmin:lonmax}) , 1., 1., 0)   
  end if
  
  x_lon = dim_avg_n_Wrap( xh(model|:,{lat|latmin:latmax}, lon|:) ,1) 
  
  r_lon = x_lon(0:1,:)
  dim = dimsizes(x_lon)
  do i=0,dim(1)-1
    r_lon(0,i) = escorc(x_lon(:,i),x1)
    r_lon(1,i) = escorc(x_lon(:,i),x2)
  end do
  
 
  print ("Plotting...")
  
  wks  = gsn_open_wks( "pdf", get_script_prefix_name())
;  gsn_merge_colormaps(wks,"CBR_drywet","amwg_blueyellowred")
 
;  cmap = gsn_retrieve_colormap( wks)
;  col = (/2,3,4,5,6, 8,9,10,11,12/)
;  lev = fspan( -0.8, 0.8, 9)
;  lab = sprintf( "%3.1f", lev)
;  lbcol = cmap(col,:)

;  col1 = (/2,3,4,5,6,7,8,9,    10,11,12,13,14,15,16,17/) + 11
;  lev1 = fspan( -1.05, 1.05, 15)
;  lab1 = sprintf( "%3.1f", lev1)
;  lbcol1 = cmap(col1,:)
 
  gsn_merge_colormaps(wks,"MPL_BrBG","amwg_blueyellowred")

  cmap = gsn_retrieve_colormap( wks)
  col = (/2,10,17,26,33,42,49,58,   74,82,90,98,106,113,122,129/)
  lev = fspan( -1.4, 1.4, 15)*0.5
  lab = sprintf( "%3.1f", lev)
  lbcol = cmap(col,:)

  col1 = (/2,3,4,5,6,7,8,9,    10,11,12,13,14,15,16,17/) + 128
;  lev1 = (/-5., -3., -1.5, -.8, -0.3, 0.,0.3, .8, 1.5, 3., 5./)
;  lev1 = (/-2.0,-1.5,-1.0,-0.8,-0.6,-0.4,-0.2,0, 0.2,0.4,0.6,0.8,1.0,1.5,2.0/)
;  lev1 = (/-5., -3., -1.5,-1,0 -.8, -0.4,-0.2, 0.,0.2,0.4,.8, 1.0,1.5, 3., 5./)
  lev1 = fspan(-0.7,0.7,15)
  lab1 = sprintf( "%3.1f", lev1)
  lbcol1 = cmap(col1,:) 
;---TAS
  rest = True
  rest@gsnDraw = False
  rest@gsnFrame = False
  rest@gsnRightString = ""
  rest@gsnLeftString = ""
  rest@gsnCenterString = ""
  rest@mpProjection = "Robinson"
  rest@mpPerimOn = False
  rest@mpGridAndLimbOn = True
;  rest@mpCenterLonF = 205
  rest@mpGridLineColor = "transparent"
;  rest@mpGeophysicalLineColor = "transparent"
  rest@cnFillOn = True
  rest@cnLinesOn = False
  rest@cnLevelSelectionMode = "ExplicitLevels"
  rest@cnLevels = lev
  rest@cnFillColors = col
  rest@cnLineLabelsOn = False
  rest@cnInfoLabelOn = False
  rest@lbLabelBarOn = False
  rest@vpWidthF = 0.3
  rest@vpHeightF = 0.2

;Stipple  
  res2 = True
  res2@gsnDraw = False
  res2@gsnFrame = False
  res2@gsnRightString = ""
  res2@gsnLeftString = ""
  res2@gsnCenterString = ""
  res2@cnFillOn = True
  res2@cnLinesOn = False
  res2@cnLevelSelectionMode = "ExplicitLevels"
  res2@cnLevels = (/-0.9,1.2/)
  res2@cnFillColors = (/"gray","transparent"/)
  res2@cnFillPattern = 17
  res2@cnFillScaleF = 0.5
  res2@cnLineLabelsOn = False
  res2@cnInfoLabelOn = False
  res2@lbLabelBarOn = False
  
;---Timeseries
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnRightString = ""
  res@gsnLeftString = ""
  res@gsnCenterString = ""

  res@tmXTOn = False
  res@tmYROn = False
  res@tmXBLabelFontHeightF = 0.01
  res@tmYLLabelFontHeightF = 0.01
  res@tmXBTickSpacingF = 20.
  res@tmXBTickStartF = -40
  res@tmXBMinorLengthF = 0.
  
  res@tiXAxisString = "Longitude"
  res@tiYAxisString = "Correlation"
  res@trXMaxF = 50.
  res@trXMinF = -50.
  res@trYMaxF = -0.2
  res@trYMinF = 1
  res@tiXAxisFontHeightF = 0.011
  res@tiYAxisFontHeightF = 0.011

  res@vpWidthF = 0.3
  res@vpHeightF =  0.25

;---line 
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnRightString = ""
  res@gsnLeftString = ""
  res@gsnCenterString = ""
  res@tmXTOn = False
  res@tmYROn = False
  res@tmXBLabelFontHeightF = 0.01
  res@tmYLLabelFontHeightF = 0.01
  res@tiXAxisFontHeightF = 0.012
  res@tiYAxisFontHeightF = 0.012
  res@xyLineThicknessF = 1.5
  
  
;---Text
  restext = True
  restext@txJust = "TopCenter"
  restext@txFontHeightF = 0.012
;---Textlabel

  restextl = True
  restextl@txJust = "TopLeft"
  restextl@txFontHeightF = 0.012 
  restextl@txFont   = "helvetica-bold"
  

;---Color Bar
  reslb = True
  reslb@lbPerimOn = False
  reslb@lbOrientation = "Horizontal"
  reslb@lbLabelFontHeightF = 0.009
  reslb@lbLabelAlignment = "InteriorEdges"
  reslb@lbMonoFillPattern = True
  reslb@lbFillColors = lbcol
  reslb@vpWidthF = 0.3
  reslb@vpHeightF = 0.045
  
;---Color Bar
  reslb1 = reslb
  delete(reslb1@lbFillColors)
  reslb1@lbFillColors = lbcol1

;---Polyline
;
  respolyline = True
  respolyline@gsLineColor = "white"
  respolyline@gsLineThicknessF = 1.
  
  rest@vpYF = 0.94	
  plot = new( 3, "graphic")

  label = (/"a","b","d","d"/) 

  
  res@vpXF =  0.025 + 0.3 + 0.015/2. - res@vpWidthF/2. 
  res@vpYF = 0.97 - 1 * (rest@vpHeightF + 0.12)
  res@xyLineColor = "black"
  plot(2) = gsn_csm_xy( wks, r_lon&lon,r_lon(0,:), res) 
  res@xyLineColor = "blue"
  ploto = gsn_csm_xy( wks, r_lon&lon,r_lon(1,:), res) 
  overlay(ploto,plot(2))
  restext@txFontHeightF = 0.01
  gsn_text_ndc( wks, label(2),res@vpXF-0.03, res@vpYF + 0.025, restextl)
  modelnumbers = ispan(1,17,1)
  restext@txFontHeightF = 0.011
  
  restext@txJust = "TopLeft"
  gsn_text_ndc( wks, "Corrlation of rainfall (5~S~o~N~N-15~S~o~N~N averaged)", res@vpXF+ 0.01, res@vpYF+0.02, restext)
  
  gsn_text_ndc( wks, "with tropical North Atlantic SST", res@vpXF+ 0.01, res@vpYF-0.01-0.2, restext)
  restext@txFontColor = "blue"
  gsn_text_ndc( wks, "with high-latitude North Atlantic SST", res@vpXF+ 0.01, res@vpYF-0.03-0.2, restext)
 
  print(r_lon)

  do i=0,16
    text = gsn_add_text( wks, plot(2), modelnumbers(i), x1(i), x2(i)-0.05, restext)
  end do

  draw( wks)
  frame( wks)

end



