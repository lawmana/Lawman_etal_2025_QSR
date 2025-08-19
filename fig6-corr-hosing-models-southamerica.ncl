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

  regiontype = "land"
;  regiontype = "ocean"
 
  region = "Brazil"
  region1 = "South_Tropical_Atlantic_SST"
  region2 = "Brazil_Rainfall"
  textx = "Tropical South Atlantic SST (~S~o~N~C)"
  texty = "Brazil Rainfall (mm/day)"
  label = (/"","c","d"/)

;  region = "Central_Andes"
;  region1 = "Eastern_Pacific_Meridional_SST_Gradient"
;  region2 = "Central_Andes_Rainfall"
  
;  region = "Northern_Andes"
;  region1 = "Pacific_Caribbean_SST"
;  region2 = "Northern_Andes_Rainfall"

  region = "Northernmost_South_America"
  region1 = "North_Tropical_Atlantic_SST"
  region2 = "Northernmost_South_America_Rainfall"
  textx = "Tropical North Atlantic SST (~S~o~N~C)"
  texty = "Northern South America Rainfall (mm/day)"
  label = (/"","a","b"/)


  if (region .eq. "Brazil") then
  	latmin = -18.
  	latmax = 0.
    lonmin = -50.
  	lonmax = -35.
  end if
  if (region .eq. "Northernmost_South_America") then
  	latmin = 0.
  	latmax = 10.
    lonmin = -80.
  	lonmax = -50.
  end if
  if ( region .eq. "Central_Andes") then
  	latmin = -25.
  	latmax = -15.
    lonmin = -70.
  	lonmax = -65.
  end if
  if ( region .eq. "Northern_Andes") then
  	latmin = -15.
  	latmax = 0.
    lonmin = -70.
  	lonmax = -80.
  end if

  if (region2 .eq. "Brazil_Rainfall") then
  	latmin2 = -18.
  	latmax2 = 0.
    lonmin2 = -50.
  	lonmax2 = -35.
  end if
  if (region2 .eq. "Northernmost_South_America_Rainfall") then
  	latmin2 = 0.
  	latmax2 = 10.
    lonmin2 = -80.
  	lonmax2 = -50.
  end if
  if (region2 .eq. "Central_Andes_Rainfall") then
  	latmin2 = latmin
  	latmax2 = latmax
    lonmin2 = lonmin
  	lonmax2 = lonmax
  end if
  
  if (region2 .eq. "Northern_Andes_Rainfall") then
  	latmin2 = latmin
  	latmax2 = latmax
    lonmin2 = lonmin
  	lonmax2 = lonmax
  end if
  
  if (region1 .eq. "South_Tropical_Atlantic_SST") then
  	latmin1 = -22.
  	latmax1 = -12.
    lonmin1 = -35.
  	lonmax1 = 10.
  end if
  
  if (region1 .eq. "North_Tropical_Atlantic_SST") then
  	latmin1 = 12.
  	latmax1 = 22.
    lonmin1 = -80.
  	lonmax1 = -40.
  end if
  
  if (region1 .eq. "Eastern_Pacific_Meridional_SST_Gradient") then
  	latmin1 = -25.
  	latmax1 = -5.
    lonmin1 = -80.
  	lonmax1 = -120.
    latmin12 = 5.
  	latmax12 = 25.
    lonmin12 = -80.
  	lonmax12 = -120.
  end if
  
  if (region1 .eq. "Pacific_Caribbean_SST") then
    latmin1 = 10.
  	latmax1 = 30.
    lonmin1 = -120.
  	lonmax1 = -60.
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
    
  if(region1 .eq. "Eastern_Pacific_Meridional_SST_Gradient")then
    x1 = wgt_areaave_Wrap( th(model|:,{lat|latmin1:latmax1}, {lon|lonmin1:lonmax1}) , 1., 1., 0) - wgt_areaave_Wrap( th(model|:,{lat|latmin12:latmax12}, {lon|lonmin12:lonmax12}) , 1., 1., 0)  
    print("Done")
  else
    x1 = wgt_areaave_Wrap( th(model|:,{lat|latmin1:latmax1}, {lon|lonmin1:lonmax1}) , 1., 1., 0)   
  end if
    
  x2 = wgt_areaave_Wrap( rh(model|:,{lat|latmin2:latmax2}, {lon|lonmin2:lonmax2}) , 1., 1., 0)   
  
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
  
  if(regiontype .eq. "land")then
  	f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
    LSMASK = f->LSMASK 
    lsmask = landsea_mask( LSMASK, rh&lat, rh&lon)
    rh1 = rh
    rh1= mask( rh, lsmask.eq.1, True)    ;only land is not masked
    delete( lsmask)
    printVarSummary(rh1)
  	x = wgt_areaave_Wrap( rh1(model|:,{lat|latmin:latmax}, {lon|lonmin:lonmax}) , 1., 1., 0)   
  end if
  print(x)
;---
  rr = rh(model|0,lat|:,lon|:)
  tr = rr

  rr = (/ escorc( x, rh(lat|:,lon|:,model|:)) /)
  tr = (/ escorc( x, th(lat|:,lon|:,model|:)) /)
  
  if(region .eq. "Northernmost_South_America" .or. region .eq. "Maritime_Continent")then
    rr = -rr
    tr = -tr
  end if

  prob = rtest(rr, 18, 0)
  rrtest = rr
  rrtest = mask(rr, prob .le. (1-slev) .or. prob .ge. slev, True)
  rr_sig = rr
  rr_sig = where( ismissing(rrtest), 1., -1.)


  prob = rtest(tr, 18, 0)
  trtest = tr
  trtest = mask(rr, prob .le. (1-slev) .or. prob .ge. slev, True)
  tr_sig = tr
  tr_sig = where( ismissing(trtest), 1., -1.)

  printVarSummary( rr)


  latmintna = 12.
  latmaxtna = 22.
  lonmintna = -80.
  lonmaxtna = -40.

  f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
  LSMASK = f->LSMASK 
  lsmask = landsea_mask( LSMASK, th&lat, th&lon)
  th1 = th
  th1 = mask( th, lsmask.gt.0.5, True)
  delete( lsmask)
  printVarSummary(th1)
  xt1  = wgt_areaave_Wrap( th1(model|:,{lat|latmintna:latmaxtna}, {lon|lonmintna:lonmaxtna}) , 1., 1., 0)   

  print ("Plotting...")
  delete(x)
  
  wks  = gsn_open_wks( "pdf", get_script_prefix_name() + "-"+ region2)
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
 
  ;gsn_merge_colormaps(wks,"MPL_BrBG","amwg_blueyellowred")
  gsn_merge_colormaps(wks,"temp_19lev","amwg_blueyellowred")

  cmap = gsn_retrieve_colormap( wks)
  col = (/2,3,4,5,6,7,8,9,    10,11,12,13,14,15,16,17/)+20
  lev = fspan( -1.4, 1.4, 15) * 3.
  lab = sprintf( "%3.1f", lev)
  lbcol = cmap(col,:)

  col1 = ispan(2,21,1)
;  lev1 = (/-5., -3., -1.5, -.8, -0.3, 0.,0.3, .8, 1.5, 3., 5./)
;  lev1 = (/-2.0,-1.5,-1.0,-0.8,-0.6,-0.4,-0.2,0, 0.2,0.4,0.6,0.8,1.0,1.5,2.0/)
;  lev1 = (/-5., -3., -1.5,-1,0 -.8, -0.4,-0.2, 0.,0.2,0.4,.8, 1.0,1.5, 3., 5./)
  lev1 = fspan(-0.9,0.9,19)
  lab1 = sprintf( "%3.1f", lev1)
  lbcol1 = cmap(col1,:) 

  ncol = dimsizes( col)
  ncol1 = dimsizes(col)

;---TAS
  rest = True
  rest@mpMinLatF = -35.
  rest@mpMaxLatF = 35.
  rest@mpMinLonF = -90.
  rest@mpMaxLonF = 20.
  rest@gsnAddCyclic = False
  rest@gsnMajorLonSpacing = 30.
  rest@gsnMinorLonSpacing = 10.
  rest@gsnMajorLatSpacing = 15.
  rest@gsnMinorLatSpacing = 5.
  rest@tmYLMajorLengthF = 0.003
  rest@tmYLMajorOutwardLengthF = 0.003
  rest@tmXBMajorLengthF = 0.003
  rest@tmXBMajorOutwardLengthF = 0.003
  rest@tmYLMinorLengthF = 0.0025
  rest@tmYLMinorOutwardLengthF = 0.0025
  rest@tmXBMinorLengthF = 0.0025
  rest@tmXBMinorOutwardLengthF = 0.0025
  rest@tmXTOn = False
  rest@tmYROn = False
  rest@tmXBLabelFontHeightF = 0.011
  rest@tmYLLabelFontHeightF = 0.011
  rest@gsnDraw = False
  rest@gsnFrame = False
  rest@gsnRightString = ""
  rest@gsnLeftString = ""
  rest@gsnCenterString = ""
  rest@mpPerimOn = False
  rest@mpGridAndLimbOn = True
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
  rest@vpWidthF = 0.4
  rest@vpHeightF = rest@vpWidthF/2.

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
  res2@cnFillScaleF = 1.
  res2@cnLineLabelsOn = False
  res2@cnInfoLabelOn = False
  res2@lbLabelBarOn = False
  
;---Scatter
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnRightString = ""
  res@gsnLeftString = ""
  res@gsnCenterString = ""
  res@gsnXRefLine = 0.0
  res@gsnYRefLine = 0.0
  res@gsnXRefLineColor = "gray60"
  res@gsnYRefLineColor = "gray60"

  res@tmXTOn = False
  res@tmYROn = False
  res@tmXBLabelFontHeightF = 0.01
  res@tmYLLabelFontHeightF = 0.01

  res@trXMaxF = max(x1)+0.2
  res@trXMinF = min(x1)-0.2
  res@trYMaxF = max(x2)+0.2
  res@trYMinF = min(x2)-0.2

  res@xyMarkerSizeF = 0.014
  res@tiXAxisFontHeightF = 0.011
  res@tiYAxisFontHeightF = 0.011
  res@xyMarkLineMode = "Markers"
  res@xyMarker = 16
  ;res@xyMonoMarker = False
  res@xyMonoMarkerColor = False 
  res@vpWidthF = 0.26
  res@vpHeightF =  res@vpWidthF 


;---Text
  restext = True
  restext@txJust = "TopCenter"
  restext@txFontHeightF = 0.012

;---Text2
  restext2 = True
  restext2@txJust = "TopCenter"
  restext2@txFontHeightF = 0.011
  restext2@txFont = "helvetica-bold"

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
  respolyline@gsLineThicknessF = 2.
  
  rest@vpYF = 0.94	
  plot = new( 3, "graphic")

  ;gsn_text_ndc( wks, "Correlation Maps based on " + region + " " +var +" Index", 0.025 + 0.3 + 0.015/2., 0.97, restext)

  restext@txAngleF = 0

  ;rest@vpXF = 0.025
  ;rest@vpYF = 0.94 - 0 * (rest@vpHeightF + 0.03)
  ;plot(0) = gsn_csm_contour_map( wks, rr(:,:), rest)
  ;plotrstp = gsn_csm_contour(wks, rr_sig(:,:), res2)   ;stipple
  ;overlay( plot(0), plotrstp)
  ;gsn_text_ndc( wks, label(0),rest@vpXF+.01, rest@vpYF , restextl)
  ;gsn_labelbar_ndc( wks, dimsizes(col), lab, rest@vpXF, rest@vpYF-rest@vpHeightF-0.01, reslb)
  
  delete(rest@cnLevels)
  delete(rest@cnFillColors)
  rest@cnLevels = lev1
  rest@cnFillColors = col1
  
  rest@vpXF = 0.05
  rest@vpYF = 0.94 - 0 * (rest@vpHeightF + 0.04)
  plot(1) = gsn_csm_contour_map( wks, tr(:,:), rest)
  plotrstp = gsn_csm_contour(wks, tr_sig(:,:), res2)   ;stipple
  overlay( plot(1), plotrstp)
  gsn_text_ndc( wks, label(1),rest@vpXF, rest@vpYF+0.01, restextl)
  gsn_labelbar_ndc( wks, dimsizes(col1), lab1, rest@vpXF+0.05, rest@vpYF-rest@vpHeightF-0.025, reslb1)
  gsn_text_ndc( wks, "Pearson Correlation Coefficient r", 0.25, rest@vpYF-rest@vpHeightF-0.024*2-0.017, restext2)

  if(var .eq. "Rainfall")then  
    respolyline@gsLineColor = "blue"
    x0 = (/lonmin, lonmax, lonmax, lonmin, lonmin/)
    y0 = (/ latmax, latmax, latmin, latmin, latmax/)
    polyline0 = gsn_add_polyline( wks, plot(1), x0, y0, respolyline)
    respolyline@gsLineColor = "white"  
    x0 = (/lonmin1, lonmax1, lonmax1, lonmin1, lonmin1/)
    y0 = (/ latmax1, latmax1, latmin1, latmin1, latmax1/)
    polyline1 = gsn_add_polyline( wks, plot(1), x0, y0, respolyline) 
    
    if(region1 .eq. "Eastern_Pacific_Meridional_SST_Gradient")then
      respolyline@gsLineColor = "white"
      x0 = (/lonmin12, lonmax12, lonmax12, lonmin12, lonmin12/)
      y0 = (/ latmax12, latmax12, latmin12, latmin12, latmax12/)
      polyline2 = gsn_add_polyline( wks, plot(1), x0, y0, respolyline) 
    end if
  end if
  
  
  res@vpXF =  0.12
  res@vpYF = 0.97 - 1 * (rest@vpHeightF + 0.15)
  res@tiXAxisString = textx
  res@tiYAxisString = texty

  k0 := ispan(0,ncol-1,1)   ;reorder from light to dark colors

  k2 := k0
  print(k2)
  k2(0::2) = k0(ncol/2-1:0) 
  k2(1::2) = k0(ncol/2:ncol-1)
  print(k2)
  print(lev)


  res@xyMarkerColors := col(k2)
  ;res@xyMarkers = (/16,16, 16, 16, 6, 6, 4,4,4,4/)

  xs = x1(0:9)
  ys = x2(0:9)
  zs = xt1(0:9)

  xw = x1(10:17)
  yw = x2(10:17)
  zw = xt1(10:17)

  print(col(k2))
  xns = new( (/ncol,dimsizes(xs)/), "double")
  yns = new( (/ncol,dimsizes(ys)/), "double")
  xnw = new( (/ncol,dimsizes(xw)/), "double")
  ynw = new( (/ncol,dimsizes(yw)/), "double")
 
  do k=0,ncol-1
    if ( k.eq.0) then
      kk = ind( zs.lt.lev(k))
    end if
    if ( k.ge.1 .and. k.le.ncol-2) then
      kk = ind( zs.ge.lev(k-1) .and. zs.lt.lev(k))
    end if
    if ( k.eq.ncol-1) then
      kk = ind( zs.ge.lev(k-1))
    end if
    if ( .not.ismissing(kk(0))) then
      xns(k,kk) = (/ xs(kk) /)
      yns(k,kk) = (/ ys(kk) /)
    end if
    delete( kk)
  end do

  do k=0,ncol-1
    if ( k.eq.0) then
      kk = ind( zw.lt.lev(k))
    end if
    if ( k.ge.1 .and. k.le.ncol-2) then
      kk = ind( zw.ge.lev(k-1) .and. zw.lt.lev(k))
    end if
    if ( k.eq.ncol-1) then
      kk = ind( zw.ge.lev(k-1))
    end if
    if ( .not.ismissing(kk(0))) then
      xnw(k,kk) = (/ xw(kk) /)
      ynw(k,kk) = (/ yw(kk) /)
    end if
    delete( kk)
  end do

  xns = xns(k2,:)
  yns = yns(k2,:)
  xnw = xnw(k2,:)
  ynw = ynw(k2,:)
  print(xns)
  print(yns)
  
  res@xyMarkerThicknessF =3.
  res@xyMarker=4
  plot(2) = gsn_csm_xy( wks, xns, yns, res)

  res@xyMarker=6
  plot_overlay = gsn_csm_xy( wks, xnw, ynw, res)
  overlay(plot(2), plot_overlay)

  restext@txFontHeightF = 0.01
  gsn_text_ndc( wks, label(2),res@vpXF-0.07, res@vpYF + 0.03, restextl)
  rx1x2 = escorc( x1, x2)
  pr = rtest(rx1x2, 18, 0)

  gsn_text_ndc( wks, "r = " + sprintf( "%4.2f", rx1x2)+" (p="+sprintf( "%4.2f", pr)+")", res@vpXF+ 0.06, res@vpYF-0.01, restext)
  modelnumbers = ispan(1,18,1)
  restext@txFontHeightF = 0.009
  ;restext@txFontColor = "gray50"
  do i=0,17
    text = gsn_add_text( wks, plot(2), modelnumbers(i), x1(i), x2(i)+0.05, restext)
  end do


  gsn_labelbar_ndc( wks, dimsizes(col), lab, rest@vpXF+0.05, res@vpYF-res@vpHeightF-0.05, reslb)
  gsn_text_ndc( wks, "marker color denotes TNA SST", 0.25, res@vpYF-res@vpHeightF-0.05-0.04, restext2)

  draw( wks)
  frame( wks)


end



