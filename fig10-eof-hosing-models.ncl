;---EOF analysis of Selected regions

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
  dir = "/glade/u/home/xianwu/ncl/UTAustin/ncl/paleo/herinch/paper/data/"
  diro = "./"

  model_names = (/"HadCM3-0.4Sv", "CCSM3-1Sv-NA", "CESM1-0.2Sv",\
  "CCSM3-0.2Sv-GIN-LR", "IPSL-CM4", "CCSM3-TraCE-MWF", "CESM1-iTraCE-MWF","IPSL-CM5", "CCSM3-0.2Sv-GIN-HR", "CESM1-0.15Sv",\
  "CCSM4-0.1Sv-RB", "COSMOS-W", "COSMOS-S", "CCSM3-0.1Sv-NA", "MIROC-W", \
     "HadCM3-0.1Sv", "CESM1-0.1Sv", "MIROC-S","Mean of 18 Model Experiments"/)
  slev = 0.90
;  percentage = True
  percentage = False
  
  var = "Rainfall"
;  var = "TAS"

  regiontype = "land"
;  regiontype = "ocean"
 
  region = "Maritime_Continent"
  ;region = "East_Asia"
  
  region1 = "Southeast_Equatorial_Pacific_Rainfall"
  region2 = "Maritime_Continent_Rainfall"
  region3 = "Indian_Ocean_Rainfall"
  if (region .eq. "Maritime_Continent") then
  	latmin = -10.
  	latmax = 10.
    lonmin = 95.
  	lonmax = 125.
  end if
  
  if (region .eq. "East_Asia") then
  	latmin = 30.
  	latmax = 80.
    lonmin = 60.
  	lonmax = 140.
  end if

print("Reading data files...")

;---Rainfall   

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
    
;---Intermodel EOF analysis
  if(var .eq. "Rainfall")then
    x = rh
  end if

;  if(regiontype .eq. "land")then
;  	f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
;    LSMASK = f->LSMASK 
;    lsmask = landsea_mask( LSMASK, rh&lat, rh&lon)
;    rh1 = rh
;    rh1= mask( rh, lsmask.eq.1, True)    ;only land is not masked
;    delete( lsmask)
;    printVarSummary(rh1)
;  	x = wgt_areaave_Wrap( rh1(model|:,{lat|latmin:latmax}, {lon|lonmin:lonmax}) , 1., 1., 0)   
;  end if

  if(region .ne. "Global")then
    x := x(model|:,{lat|latmin:latmax},{lon|lonmin:lonmax})
  end if
  
  eof = x(0:2,:,:)
  eof(:,:,:) = eofunc_n_Wrap( x(:,:,:), 3, False, 0)


  ;eof(1,:,:) = -eof(1,:,:)

  
  printVarSummary(eof)
  print(eof@pcvar)

  pcvar = eof@pcvar
  print(pcvar)
  
  ts = eofunc_ts_n_Wrap( x, eof, False, 0)
  printVarSummary(ts)
  
;---Regression
  
  ts_std = dim_standardize( ts, 1)  ;evn, S
  rr = rh(0:2,:,:) 
  rt = rr
  rr(:,:,:) = (/ regCoef_n( ts_std (:,:),rh(:,:,:),1,0) /)
  rt(:,:,:) = (/ regCoef_n( ts_std (:,:),th(:,:,:),1,0) /)  
  
  rt!0="env"
  rr!0="env"
  
  ts_std!1 = "model"
  ts_std&model = ispan(1,18,1)
  
  rt&env = ispan(1,3,1)
  rr&env = ispan(1,3,1)

  printVarSummary(rr)
  
  
  print ("Plotting...")
  
  wks  = gsn_open_wks( "pdf", get_script_prefix_name() + "-"+ region)
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
  rest@mpCenterLonF = 205
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
  res2@cnFillScaleF = 0.5
  res2@cnLineLabelsOn = False
  res2@cnInfoLabelOn = False
  res2@lbLabelBarOn = False
  
;---PC  
  res = True
  res@gsnFrame = False
  res@gsnDraw = False
  res@gsnRightString = ""
  res@gsnLeftString = ""
  res@gsnCenterString  = "" 
  res@gsnXYBarChart = True

  res@gsnXRefLineColor = "Gray"
  res@gsnXRefLineDashPattern = 2
  res@gsnYRefLine = 0.
  res@gsnAboveYRefLineColor = "red"
  res@gsnBelowYRefLineColor = "blue"
;  res@gsnXYBarChartBarWidth = 0.6
  res@tmXTOn = False
  res@tmYROn = False
;  res@tmXBMode = "Explicit"
;  res@tmXBLabels = modelname
;  res@tmXBValues = ispan( 1, modelnum, 1) 
  res@tmXBLabelAngleF = 90.
  res@tmXBMinorPerMajor = 9
  res@tmXBLabelDeltaF = -0.4
  res@tmXBLabelFontHeightF = 0.008
  res@tmXBMajorLengthF = 0.003
  res@tmXBMajorOutwardLengthF = 0.003
  res@tmXBMinorLengthF = 0.003
  res@tmXBMinorOutwardLengthF = 0.003
  res@tmXBLabelJust = "CenterCenter"
  res@tmYLMajorLengthF = 0.005
  res@tmYLMajorOutwardLengthF = 0.005
  res@tmYLMinorLengthF = 0.005
  res@tmYLMinorOutwardLengthF = 0.005
  res@tmYLMinorOn = False
  res@tmYLLabelDeltaF = -0.4
  res@tmYLLabelFontHeightF = 0.01
  res@tmYLPrecision = 1.

  res@trYMinF = -4.
  res@trYMaxF =4.
  res@tiXAxisString = ""
  res@tiYAxisString = ""
  res@xyLineColor = "transparent"
  res@vpWidthF = 0.26
  res@vpHeightF = 0.11
  res@vpXF = 0.55 
  

;---Text
  restext = True
  restext@txJust = "TopCenter"
  restext@txFont   = "helvetica-bold"
  restext@txFontHeightF = 0.012
;---Textlabel
  restextl = True
  restextl@txJust = "TopCenter"
  restextl@txFontHeightF = 0.012 
  restextl@txFont   = "helvetica-bold"
  
;---Textlabel
  restextl2 = restextl
  restextl2@txFontHeightF = 0.008
  
;---Color Bar
  reslb = True
  reslb@lbPerimOn = False
  reslb@lbOrientation = "Horizontal"
  reslb@lbLabelFontHeightF = 0.009
  reslb@lbLabelAlignment = "InteriorEdges"
  reslb@lbMonoFillPattern = True
  reslb@lbFillColors = lbcol
  reslb@vpWidthF = 0.3
  reslb@vpHeightF = 0.04
  
;---Color Bar
  reslb1 = reslb
  delete(reslb1@lbFillColors)
  reslb1@lbFillColors = lbcol1

;---Polyline
;
  respolyline = True

  
  rest@vpYF = 0.94	
  plot1 = new( 4, "graphic")
  plot2 = new( 4, "graphic")
  plot3 = new( 4, "graphic")

  label1 = (/"a","d"/)
  label2 = (/"b","e"/)
  label3 = (/"c","f"/)

  gsn_text_ndc( wks, "Inter-model EOF Analysis of Maritime Continent " +var, 0.025 + 0.3 + 0.015/2., 0.97, restext)

  restext@txAngleF = 0
  
  do i=0,1
    rest@vpXF = 0.025
    rest@vpYF = 0.94 - i * (rest@vpHeightF + 0.03)
    plot1(i) = gsn_csm_contour_map( wks, rr(i,:,:), rest)
    gsn_text_ndc( wks, label1(i),rest@vpXF+.01, rest@vpYF , restextl)
  end do
  
  gsn_labelbar_ndc( wks, dimsizes(col), lab, rest@vpXF, rest@vpYF-rest@vpHeightF-0.01, reslb)
  gsn_text_ndc( wks, "Regression of PREC changes on EOF PCs (mm/day/std)", rest@vpXF+rest@vpWidthF/2., rest@vpYF-rest@vpHeightF-0.045, restextl2)


  if(var .eq. "Rainfall")then  
    x0 = (/lonmin, lonmax, lonmax, lonmin, lonmin/)
    y0 = (/latmax, latmax, latmin, latmin, latmax/)
    respolyline@gsLineColor = "White"
    respolyline@gsLineThicknessF = 3.
    polyline0 = gsn_add_polyline( wks, plot1(0), x0, y0, respolyline)
    respolyline@gsLineColor = "Blue"
    respolyline@gsLineThicknessF = 1.
    polyline1 = gsn_add_polyline( wks, plot1(0), x0, y0, respolyline)
  end if
  
  delete(rest@cnLevels)
  delete(rest@cnFillColors)
  rest@cnLevels = lev1
  rest@cnFillColors = col1
  
  do i=0,1
    rest@vpXF = 0.025 + 0.3 + 0.015
    rest@vpYF = 0.94 - i * (rest@vpHeightF + 0.03)
    plot2(i) = gsn_csm_contour_map( wks, rt(i,:,:), rest)
    gsn_text_ndc( wks, label2(i),rest@vpXF+.01, rest@vpYF, restextl)
  end do
  gsn_labelbar_ndc( wks, dimsizes(col1), lab1, rest@vpXF, rest@vpYF-rest@vpHeightF-0.01, reslb1)
  gsn_text_ndc( wks, "Regression of TS changes on EOF PCs (~S~o~N~C/std)", rest@vpXF+rest@vpWidthF/2., rest@vpYF-rest@vpHeightF-0.045, restextl2)

  do i=0,1
    res@vpXF = 0.025 + 0.3*2+0.04
    res@vpYF = 0.94 - i * (rest@vpHeightF + 0.03)-0.02
    print(ts_std)
    gsn_text_ndc( wks, label3(i),res@vpXF, res@vpYF+0.02 , restextl)
    PCname = "  standardized PC"+(i+1)+" ("+sprintf("%2.0f",pcvar(i)) + "%)"
    if(i .eq. 1)then
      res@tmXBOn = True
      res@tmXBMode = "Explicit"
      res@tmXBLabels = model_names
      res@tmXBValues = ispan( 1, 18, 1) 
    else
      res@tmXBOn = False

    end if
    gsn_text_ndc( wks,PCname, res@vpXF+0.1, res@vpYF+0.02, restext)
    plot3(i) = gsn_csm_xy( wks, ts_std&model, ts_std(i,:), res)
  end do

  draw( wks)
  frame( wks)

end



