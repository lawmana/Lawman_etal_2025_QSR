
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
  dir = "/glade/u/home/xianwu/ncl/UTAustin/ncl/paleo/herinch/paper/data/"
  diro = "./"
  
  std = True
  std = False
  
;  percentage = True
  percentage = False
  
  regionr1 = "India_Rainfall"
  region2 = "Tropical_Northern_Atlantic_TAS"
  region3 = "Northern_Hemisphere_TAS"

  regiont1 = "Arabian_Sea_TAS"
  region2 = "Tropical_Northern_Atlantic_TAS"
  region3 = "Northern_Hemisphere_TAS"

  print("Reading data files...")

 
  f = addfile( dir + "Regridded_precip_18models.nc", "r")
  rh = f->hose
  rc = f->control

  f = addfile( dir + "Regridded_tas_18models.nc", "r")  
  th = f->hose
  tc = f->control

  if (any(isnan_ieee(rh))) then
    value = 1.e20
    replace_ieeenan (rh, value, 0)
    rh@_FillValue = value
  end if
  if (any(isnan_ieee(rc))) then
    value = 1.e20
    replace_ieeenan (rc, value, 0)
    rc@_FillValue = value
  end if

  if (any(isnan_ieee(th))) then
    value = 1.e20
    replace_ieeenan (th, value, 0)
    th@_FillValue = value
  end if
  if (any(isnan_ieee(tc))) then
    value = 1.e20
    replace_ieeenan (tc, value, 0)
    tc@_FillValue = value
  end if


  thave = dim_avg_n_Wrap(th, 0)
  
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
  
;---Region1 (Color)

  latmin = 0.
  latmax = 35.
  lonmin = 70.
  lonmax = 90.
  	
;  	latmin = 10.
;  	latmax = 30.
;  	lonmin = 70.
;  	lonmax = 85.
  f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
  LSMASK = f->LSMASK 
  lsmask = landsea_mask( LSMASK, rh&lat, rh&lon)
  rh= mask( rh, lsmask.le.0.5, True)
  delete( lsmask)
  xr1 = wgt_areaave_Wrap( rh(model|:,{lat|latmin:latmax}, {lon|lonmin:lonmax}) , 1., 1., 0)   
  

  latmin = 5.
  latmax = 25.
  lonmin = 50.
  lonmax = 70.
  f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
  LSMASK = f->LSMASK 
  lsmask = landsea_mask( LSMASK, th&lat, th&lon)
  th1 = th
  th1 = mask( th, lsmask.gt.0.5, True)
  delete( lsmask)
  printVarSummary(th1)
  xt1  = wgt_areaave_Wrap( th1(model|:,{lat|latmin:latmax}, {lon|lonmin:lonmax}) , 1., 1., 0)   

  px1 = (/lonmin, lonmax, lonmax, lonmin, lonmin/)
  py1 = (/ latmax, latmax, latmin, latmin, latmax/)

;---Region2
  if ( region2 .eq. "Tropical_Northern_Atlantic_TAS") then
  	;latmin = 12.
  	;latmax = 22.
  	;lonmin = -80.
  	;lonmax = -20.
  	
    latmin = 12.
  	latmax = 22.
    lonmin = -80.
  	lonmax = -40.

    x2 = wgt_areaave_Wrap( th(model|:,{lat|latmin:latmax}, {lon|lonmin:lonmax}) , 1., 1., 0)   
  end if
  
  
  if ( region2 .eq. "Tropical_Pacific_Meridonal_TAS") then
  	latmin = -20.
  	latmax = -5.
  	lonmin = -110.
  	lonmax = -80.

  	latmin1 = 5.
  	latmax1 = 20.
  	lonmin1 = -110.
  	lonmax1 = -80.

    x2 = wgt_areaave_Wrap( th(model|:,{lat|latmin:latmax}, {lon|lonmin:lonmax}) , 1., 1., 0) 
    x2 = x2- wgt_areaave_Wrap( th(model|:,{lat|latmin1:latmax1}, {lon|lonmin1:lonmax1}) , 1., 1., 0)
      
    px3 = (/lonmin1, lonmax1, lonmax1, lonmin1, lonmin1/)
    py3 = (/ latmax1, latmax1, latmin1, latmin1, latmax1/)
  end if 
  
  px2 = (/lonmin, lonmax, lonmax, lonmin, lonmin/)
  py2 = (/ latmax, latmax, latmin, latmin, latmax/)

;---Region3

;  if(region3 .eq. "India_Rainfall")then
;  	latmin = 0.
;  	latmax = 30.
;  	lonmin = 60.
;  	lonmax = 100.
;  	x3 = wgt_areaave_Wrap( rh(model|:,{lat|latmin:latmax}, {lon|lonmin:lonmax}) , 1., 1., 0)   
;  end if

  if ( region3 .eq. "Northern_Hemisphere_TAS") then
  	latmin = 45.
  	latmax = 90.
  	lonmin = -180.
  	lonmax = 180.
    x3 = wgt_areaave_Wrap( th(model|:,{lat|latmin:latmax}, {lon|lonmin:lonmax}) , 1., 1., 0)   
  end if 
  
  px3 =(/lonmin, lonmax, lonmax, lonmin, lonmin/)
  py3 = (/ latmax, latmax, latmin, latmin, latmax/)
  
  
  if(std)then
   x1 = dim_standardize(x1,1)
   x2 = dim_standardize(x2,1)
   x3 = dim_standardize(x3,1)
  end if

  wks  = gsn_open_wks( "pdf", get_script_prefix_name() + where(percentage, "-percentage", ""))

  cgry = "gray50"
;  gsn_define_colormap( wks, "MPL_BrBG")
;  cmap = gsn_retrieve_colormap( wks)
;  col = (/2,10,17,26,33,42,49,58,   74,82,90,98,106,113,122,129/)
;  lev = fspan(-0.7,0.7,15)

;  gsn_define_colormap( wks, "amwg_blueyellowred")
;  cmap = gsn_retrieve_colormap( wks) 
;  col = ispan( 2, 17, 1)
;  lev = fspan( -1.75, 1.75, 15) * 1.   ;-2.6 < n0 < 2.6
 
  gsn_merge_colormaps(wks,"MPL_BrBG","amwg_blueyellowred")
  cmap = gsn_retrieve_colormap( wks)
  col = (/2,10,17,26,33,42,49,58,   74,82,90,98,106,113,122,129/)
  lev = fspan(-0.7,0.7,15)
  lab = sprintf( "%3.1f", lev)
  lbcol = cmap(col,:)

  col1 = (/2,3,4,5,6,7,8,9,    10,11,12,13,14,15,16,17/) + 128
  lev1 = fspan( -1.75, 1.75, 15) * 1.
  lab1 = sprintf( "%3.1f", lev1)
  lbcol1 = cmap(col1,:)
  ncol = dimsizes( col)
  ncol1 = dimsizes(col)
;---Scatter
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnRightString = ""
  res@gsnLeftString = ""
  res@gsnCenterString = ""
  res@gsnXRefLineColor = "Grey"
  res@gsnYRefLineColor = "Grey"
  res@gsnXRefLine = 0.0
  res@gsnYRefLine = 0.0
  res@tmXTOn = False
  res@tmYROn = False
  res@tmXBLabelFontHeightF = 0.012
  res@tmYLLabelFontHeightF = 0.012
  res@trXMaxF = max(x2)+0.5
  res@trXMinF = min(x2)-0.2
  res@trYMaxF = max(x3)+0.5
  res@trYMinF = min(x3)-0.5
  res@xyMarkerSizeF = 0.016
  res@tiXAxisFontHeightF = 0.012
  res@tiYAxisFontHeightF = 0.012
  res@xyMarkLineMode = "Markers"
  res@xyMarker = 16
  res@xyMonoMarkerColor = False

;---Color Bar
  reslb = True
  reslb@lbPerimOn = False
  reslb@lbOrientation = "Horizontal"
  reslb@lbLabelFontHeightF = 0.01
  reslb@lbLabelAlignment = "InteriorEdges"
  reslb@lbMonoFillPattern = True
  print(lbcol)
;  reslb@lbFillColors = lbcol(0:8,:)
  reslb@lbFillColors = lbcol
  reslb@vpWidthF = 0.3
  reslb@vpHeightF = 0.04
  
;---Color Bar
  reslb1 = reslb
  delete(reslb1@lbFillColors)
  reslb1@lbFillColors = lbcol1
  
;---Text
  restext = True
  restext@txJust = "TopCenter"
  restext@txFontHeightF = 0.012

;---Textlabel
  restextl = True
  restextl@txJust = "TopLeft"
  restextl@txFontHeightF = 0.012 
  restextl@txFont   = "helvetica-bold"

   
  modelnumbers = ispan(1,18,1)
  res@vpYF = 0.94	
  res@vpHeightF = 0.3
  res@vpWidthF = 0.3
  plot = new( 1, "graphic")

;  gsn_text_ndc( wks, "Scatter Plot", 0.025 + 0.3 + 0.015/2., res@vpYF + 0.04, restext)

  restext@txAngleF = 0
  k0 = ispan(0,ncol-1,1)   ;reorder from light to dark colors

  k2 = k0
  print(k2)
  k2(0::2) = k0(ncol/2-1:0) 
  k2(1::2) = k0(ncol/2:ncol-1)
  print(k2)
  print(lev)

  res@xyMarkerColors = col(k2)
  print(col(k2))
  
;---Rainfall
  res@vpXF = 0.2

  xs = x2(0:9)
  ys = x3(0:9)
  zs = xr1(0:9)

  xw = x2(10:17)
  yw = x3(10:17)
  zw = xr1(10:17)


  res@vpWidthF = 0.32
  res@vpHeightF = 0.32
  res@vpXF = 0.45 - res@vpWidthF

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
  print(xnw)
  print(ynw)

  res@tiXAxisString = ""
  res@tiYAxisString = ""
  res@tiXAxisString = "Tropical North Atlantic SST (~S~o~N~C)"
  res@tiYAxisString = "Northern Hemisphere Surface Temperature (~S~o~N~C)"

  res@xyMarkerThicknessF =5.
  res@xyMarker=4
  plot = gsn_csm_xy( wks, xns, yns, res)
  res@xyMarker=6
  plot_overlay = gsn_csm_xy( wks, xnw, ynw, res)
  overlay(plot, plot_overlay)
 
  restext@txFontHeightF = 0.01
  restext@txJust = "TopRight"
;  gsn_text_ndc( wks,     "r(x,y) = " + sprintf( "%4.2f", escorc( x2, x3)), res@vpXF+ 0.15, res@vpYF-0.01, restext)
;  gsn_text_ndc( wks,     "p(x,y) = " + sprintf( "%4.2f", rtest( escorc( x2, x3), 14, 0)), res@vpXF+ 0.15, res@vpYF-0.025, restext)
;  gsn_text_ndc( wks, "r(x,color) = " + sprintf( "%4.2f", escorc( x2, x1)), res@vpXF+ 0.15, res@vpYF-0.04, restext)
;  gsn_text_ndc( wks, "p(x,color) = " + sprintf( "%4.2f", rtest( escorc( x2, x1), 14, 0)), res@vpXF+ 0.15, res@vpYF-0.055, restext)
;  gsn_text_ndc( wks, "r(y,color) = " + sprintf( "%4.2f", escorc( x3, x1)), res@vpXF+ 0.15, res@vpYF-0.07, restext)
;  gsn_text_ndc( wks, "p(y,color) = " + sprintf( "%4.2f", rtest( escorc( x3, x1), 14, 0)), res@vpXF+ 0.15, res@vpYF-0.085, restext)
  r1 = escorc( x2, x3)
  pr1 = rtest(r1, 18, 0)
  r2 = escorc( x2, xr1)
  pr2 = rtest(r2, 18, 0)
  r3 = escorc( x3, xr1)
  pr3 = rtest(r3, 18, 0)
  gsn_text_ndc( wks, "e", res@vpXF, res@vpYF+0.015, restextl)
  gsn_text_ndc( wks,     "r(x,y) = " + sprintf( "%4.2f", r1)+" (p="+sprintf( "%4.2f", pr1)+")", res@vpXF+ 0.16, res@vpYF-0.01, restext)
  gsn_text_ndc( wks, "r(x,color) = " + sprintf( "%4.2f", r2)+" (p="+sprintf( "%4.2f", pr2)+")", res@vpXF+ 0.16, res@vpYF-0.03, restext)
  gsn_text_ndc( wks, "r(y,color) = " + sprintf( "%4.2f", r3)+" (p="+sprintf( "%4.2f", pr3)+")", res@vpXF+ 0.16, res@vpYF-0.05, restext)

  printVarSummary(x2)
  printVarSummary(x3)
 
  restext@txJust = "CenterCenter"
  restext@txFontHeightF = 0.008
  do i=0,17
    if(i .eq. 3)then 
      text = gsn_add_text( wks, plot, modelnumbers(i), x2(i), x3(i), restext)
    else if (i .eq. 10)then
      text = gsn_add_text( wks, plot, modelnumbers(i), x2(i), x3(i), restext)
    else
      text = gsn_add_text( wks, plot, modelnumbers(i), x2(i), x3(i), restext)
    end if
    end if
  end do

  restext@txJust = "TopCenter"
  restext@txFontHeightF = 0.011
  restext@txFont = "helvetica-bold"
  gsn_text_ndc( wks, "marker color denotes India rainfall (mm/day)", 0.29, res@vpYF-res@vpHeightF-0.1, restext)
  gsn_labelbar_ndc( wks, dimsizes(col), lab, res@vpXF+0.01, res@vpYF-res@vpHeightF-0.06, reslb)
  
  ;---SST
  delete( xs)
  delete( ys)
  delete( zs)
  delete( xns)
  delete( yns) 
  delete( xw)
  delete( yw)
  delete( zw)
  delete( xnw)
  delete( ynw)



  draw( wks)
  frame( wks)
  delete( plot)
end



