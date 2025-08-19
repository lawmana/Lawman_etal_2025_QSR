load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  dir = "/glade/u/home/xianwu/ncl/UTAustin/ncl/paleo/herinch/paper/data/"


;------------------------------------------------------------------
  print ("Reading data files...")
;    
  model_names = (/"HadCM3-0.4Sv", "CCSM3-1Sv-NA", "CESM1-0.2Sv",\
  "CCSM3-0.2Sv-GIN-LR", "IPSL-CM4", "CCSM3-TraCE-MWF",\
   "CESM1-iTraCE-MWF","IPSL-CM5", \
  "CCSM3-0.2Sv-GIN-HR", "CESM1-0.15Sv",\
  "CCSM4-0.1Sv-RB", "COSMOS-W", "COSMOS-S", "CCSM3-0.1Sv-NA", "MIROC-W", \
   "HadCM3-0.1Sv", "CESM1-0.1Sv", "MIROC-S"/)
	 
  g = addfile(dir + "Regridded_tas_18models.nc", "r")
  hose_tas = g->hose
  control_tas = g->control

  printVarSummary(hose_tas)
  g = addfile(dir + "Regridded_precip_18models.nc", "r")
  hose_pr = g->hose
  control_pr = g->control

  if (any(isnan_ieee(hose_tas))) then
    value = 1.e20
    replace_ieeenan (hose_tas, value, 0)
    hose_tas@_FillValue = value
  end if
  if (any(isnan_ieee(control_tas))) then
    value = 1.e20
    replace_ieeenan (control_tas, value, 0)
    control_tas@_FillValue = value
  end if
  
  if (any(isnan_ieee(hose_pr))) then
    value = 1.e20
    replace_ieeenan (hose_pr, value, 0)
    hose_pr@_FillValue = value
  end if
  if (any(isnan_ieee(control_pr))) then
    value = 1.e20
    replace_ieeenan (control_pr, value, 0)
    control_pr@_FillValue = value
  end if

  pr = control_pr(model|0:1, lat|:, lon|:)  
  tas = pr
  
  pr(0,:,:) = dim_avg_n(hose_pr(0:9,:,:), 0)
  tas(0,:,:) = dim_avg_n(hose_tas(0:9,:,:), 0)

  pr(1,:,:) = dim_avg_n(hose_pr(10:17,:,:), 0) 
  tas(1,:,:) = dim_avg_n(hose_tas(10:17,:,:), 0)
  
  pr_sig = pr
  tas_sig = tas 

  printVarSummary(tas_sig)
  
  do i = 0, 95 
  do j = 0, 143
  	pr_sig(0,i,j) = where(num(hose_pr(0:9, i, j).ge.0).ge.9 .or. num(hose_pr(0:9, i, j).le.0).ge.9 , pr(0,i,j), pr_sig@_FillValue)
  	pr_sig(1,i,j) = where(num(hose_pr(10:17, i, j).ge.0).ge.7 .or. num(hose_pr(10:17, i, j).le.0).ge.7 , pr(0,i,j), pr_sig@_FillValue)

  	tas_sig(0,i,j) = where(num(hose_tas(0:9, i, j).ge.0).ge.9 .or. num(hose_tas(0:9, i, j).le.0).ge.9 , tas(0,i,j), tas_sig@_FillValue)
  	tas_sig(1,i,j) = where(num(hose_tas(10:17, i, j).ge.0).ge.7 .or. num(hose_tas(10:17, i, j).le.0).ge.7 , tas(0,i,j), tas_sig@_FillValue)
  end do
  end do
   
  data = asciiread( dir+"DiNezio_etal_HS1_database_2025_simplified.csv", (/157, 5/), "float")
  x = data(:,0)
  printVarSummary(data)
  print(data)
  
  x_flag = x 
  x_flag(:) = where(x.ge.(6/100.), 1., 0.)
  x_flag(:) = where(x.le.(-6/100.), -1., x_flag(:))
  x_flag(:) = where(abs(x).lt.(6/100.), 0., x_flag(:))
;  
;---Figures  
;
  wks  = gsn_open_wks( "pdf", get_script_prefix_name())
  gsn_merge_colormaps(wks,"MPL_BrBG","amwg_blueyellowred")

  cmap = gsn_retrieve_colormap( wks)
  col = (/2,10,17,26,33,42,49,58,   74,82,90,98,106,113,122,129/)
  lev = fspan( -1.4, 1.4, 15)
  lab = sprintf( "%3.1f", lev)
  lbcol = cmap(col,:)

  col1 = (/2,3,4,5,6,7,8,9,    10,11,12,13,14,15,16,17/) + 128
;  lev1 = (/-5., -3., -1.5, -.8, -0.3, 0.,0.3, .8, 1.5, 3., 5./)
;  lev1 = (/-2.0,-1.5,-1.0,-0.8,-0.6,-0.4,-0.2,0, 0.2,0.4,0.6,0.8,1.0,1.5,2.0/)
  lev1 = (/-5., -3., -1.5,-1,0 -.8, -0.4,-0.2, 0.,0.2,0.4,.8, 1.0,1.5, 3., 5./)
;  lev1 = fspan(-0.7,0.7,15)*5
  lab1 = sprintf( "%3.1f", lev1)
  lbcol1 = cmap(col1,:)

;---SST Anomaly 
  rest = True
  rest@gsnDraw = True
  rest@gsnFrame = False
  rest@gsnRightString = "" 
  rest@gsnLeftString = ""
  rest@gsnCenterString = "" 
  rest@gsnAddCyclic = True
  rest@mpProjection = "Robinson"
  rest@mpCenterLonF = 275.
  rest@mpPerimOn = False
  rest@mpGridAndLimbOn = True 
  rest@mpGridLineColor = "transparent"
  rest@mpGeophysicalLineColor = "black"
  rest@mpLandFillColor = "transparent"
  rest@cnFillOn = True 
;  rest@cnFillMode = "RasterFill" 
  rest@cnFillOpacityF = 1.
  rest@cnLinesOn = False
  rest@cnLevelSelectionMode = "ExplicitLevels"
  rest@cnLevels = lev
  rest@cnFillColors = col 
  rest@cnFillDrawOrder = "PreDraw"
  rest@cnLineLabelsOn = False
  rest@cnInfoLabelOn = False
  rest@lbLabelBarOn = False
  rest@tmYLLabelFontHeightF = 0.01
  rest@tmXBLabelFontHeightF = 0.01
  rest@tmXTOn = False
  rest@tmYROn = False 
  rest@tmXBOn = False
  rest@tmYLOn = False

;---SST_ci
  rest1 = True
  rest1@gsnDraw = False
  rest1@gsnFrame = False
  rest1@gsnRightString = "" 
  rest1@gsnLeftString = ""
  rest1@gsnCenterString = ""
  rest1@tmXTOn = False
  rest1@tmYROn = False 
  rest1@tmYLLabelFontHeightF = 0.008
  rest1@tmXBLabelFontHeightF = 0.008

  rest1@cnFillOn = True
;  rest1@cnConstFEnableFill = True
  rest1@cnMonoFillPattern = True
  rest1@cnFillPattern = 17
  rest1@cnFillScaleF = 0.5
  rest1@cnFillDotSizeF = 0.0002
  rest1@cnMonoFillColor = True 
  rest1@cnFillDrawOrder = "PostDraw"
  rest1@cnLinesOn = False
  rest1@cnFillColor = "gray60"
  rest1@cnLineLabelsOn = False
  rest1@cnInfoLabelOn = False
  rest1@lbLabelBarOn = False
  rest1@tmXBMajorLengthF = 0.009
  rest1@tmXBMajorOutwardLengthF = 0.009
  rest1@tmXBMinorPerMajor = 1  


;---Text
  restext = True
  restext@txJust = "TopCenter"
  restext@txFontHeightF = 0.012 
  restext@txFont   = "helvetica-bold"
;---Textlabel

  restextl = True
  restextl@txJust = "TopLeft"
  restextl@txFontHeightF = 0.012 
  restextl@txFont   = "helvetica-bold"
;---Polyline 
;
  respolyline = True
  respolyline@gsLineColor = "white"
  respolyline@gsLineThicknessF = 1.
  
;---polygon---
  
  polyres                   = True
  polyres@gsMarkerIndex     = 16          ; polymarker style
  polyres@gsMarkerSizeF     = 7          ; polymarker size

  ;---ploycover---
  polyres1                   = True
  polyres1@gsMarkerIndex     = 16          ; polymarker style
  polyres1@gsMarkerSizeF     = 12          ; polymarker size
  polyres1@gsMarkerColor = "black"  	

  ;---ploycover---
  polyres2                   = True
  polyres2@gsMarkerIndex     = 4          ; polymarker style
  polyres2@gsMarkerSizeF     = 8          ; polymarker size
  polyres2@gsMarkerColor = "black"  	

  ;---poly 
  polyresl = True
  polyresl@gsLineColor = "black" 
  polyresl@gsLineThicknessF = 1.  

;---Color Bar  
  reslb = True 
  reslb@lbPerimOn = False
  reslb@lbOrientation = "Horizontal"
  reslb@lbLabelFontHeightF = 0.012
  reslb@lbLabelAlignment = "InteriorEdges"
  reslb@lbMonoFillPattern = True
  reslb@lbFillColors = lbcol
  reslb@vpWidthF = 0.6
  reslb@vpHeightF = 0.04
  
;---Color Bar
  reslb1 = reslb
  delete(reslb1@lbFillColors)
  reslb1@lbFillColors = lbcol1
   

;--------------------------------------------------------
  rest@vpWidthF = 0.45
  rest@vpHeightF = 0.225 
  rest@vpYF = 0.98     
  rest@vpXF = 0.1  
   
  label = (/"b","d","a","c"/) 
  labelname = (/"Strong/moderate TNA cooling","Weak TNA Cooling"/)

  plot = new(20, "graphic")   
  polyline1 = new(2, "graphic")   

  latmin1 = 12.
  latmax1 = 22.
  lonmin1 = -80.
  lonmax1 = -40.
     
  do j = 0, 1
  	rest@vpYF = 0.98 - 1 * (rest@vpHeightF + 0.03)-0.055
  	rest@vpXF = 0.04 + j * (rest@vpWidthF + 0.02) 
    gsn_text_ndc( wks, label(j),rest@vpXF+0.01, rest@vpYF , restextl)
    ;gsn_text_ndc( wks, labelname(j),rest@vpXF+rest@vpWidthF/2., rest@vpYF + 0.02, restext)
  	plot(j) = gsn_csm_contour_map( wks,pr(j,:,:), rest)
 	  plot1 = gsn_csm_contour( wks, pr_sig(j,:,:), rest1)
  	overlay(plot(j), plot1) 
  end do  
  restext@txFontHeightF = 0.011
  gsn_labelbar_ndc( wks, dimsizes(col), lab, 0.2, rest@vpYF-rest@vpHeightF-0.01, reslb)
  gsn_text_ndc( wks, "Precipitation change due to freshwater hosing (mm/day)", .5, rest@vpYF-rest@vpHeightF-0.05, restext)



  plotd = new(1500, "graphic")        
  plotd1 = new(1500, "graphic")    
  plotd2 = new(1500, "graphic")      
  
  plotd = new(1500, "graphic")        
  plotd1 = new(1500, "graphic")    
  plotd2 = new(1500, "graphic")      
  
  do k = 0, 1 
  do i = 0, 153
   
  	if( data(i, 4).eq.1) then
      ;triangle
  		polyres@gsMarkerIndex     = NhlNewMarker(wks, "u", 34, 0., 0., 1., 1.5, 0.)      
  		polyres@gsMarkerSizeF     = 6          ; polymarker size
  		polyres1@gsMarkerIndex     = NhlNewMarker(wks, "u", 34, 0., 0.00008, 1., 1.5, 0.)        
  		polyres1@gsMarkerSizeF     = 8          ; polymarker size
  	else if (data(i, 4).eq.4) then
      ;dots
  		polyres@gsMarkerIndex     = 16          ; polymarker style
  		polyres@gsMarkerSizeF     = 5  ;5          ; polymarker size
  		polyres1@gsMarkerIndex     = 16          ; polymarker style   
  		polyres1@gsMarkerSizeF     = 7  ;10          ; polymarker size
  	else if (data(i, 4).eq.2) then
      ;square
  		polyres@gsMarkerIndex     = NhlNewMarker(wks, "y", 35, 0., 0., 1., 1., 0.)          ; polymarker style
  		polyres@gsMarkerSizeF     = 5  ;5          ; polymarker size
  		polyres1@gsMarkerIndex     = NhlNewMarker(wks, "y", 35, 0., 0., 1., 1., 0.)          ; polymarker style
  		polyres1@gsMarkerSizeF     = 7  ;9          ; polymarker size  
  	elseif (data(i, 4).eq.3)
      ;dimond
  		polyres@gsMarkerIndex     = NhlNewMarker(wks, "q", 35, 0., 0., 1., 1.5, 0.)          ; polymarker style
  		polyres@gsMarkerSizeF     = 4 ;4         ; polymarker size
  		polyres1@gsMarkerIndex     = NhlNewMarker(wks, "q", 35, 0., 0., 1., 1.5, 0.)          ; polymarker style
  		polyres1@gsMarkerSizeF     = 6  ;8          ; polymarker size 
    end if 
    end if    
    end if    
       
  	if( data(i, 3).eq.1.) then
  		polyres1@gsMarkerColor =  "white"   	
  		polyres@gsMarkerColor = "blue" 
  		plotd1(k*157+i) = gsn_add_polymarker( wks, plot(k), data(i, 2), data(i, 1), polyres1) 
  		plotd(k*157+i) = gsn_add_polymarker( wks, plot(k), data(i, 2), data(i, 1), polyres) 
  	else if( data(i, 3).eq.-1.) then
  		polyres1@gsMarkerColor = "white"  
  		polyres@gsMarkerColor = "saddle brown" 
  		plotd1(k*157+i) = gsn_add_polymarker( wks, plot(k), data(i, 2), data(i, 1), polyres1)  	
  		plotd(k*157+i) = gsn_add_polymarker( wks, plot(k), data(i, 2), data(i, 1), polyres)  	
  	else if( data(i, 3).eq.0.) then
  		polyres1@gsMarkerColor = "black"  	
  		polyres@gsMarkerColor = "white"
  		plotd1(k*157+i) = gsn_add_polymarker( wks, plot(k), data(i, 2), data(i, 1), polyres1)  	
  		plotd(k*157+i) = gsn_add_polymarker( wks, plot(k), data(i, 2), data(i, 1), polyres)  	
    else
  		polyres@gsMarkerColor = "black"
  		plotd1(k*157+i) = gsn_add_polymarker( wks, plot(k), data(i, 2), data(i, 1), polyres1)  	
  		plotd(k*157+i) = gsn_add_polymarker( wks, plot(k), data(i, 2), data(i, 1), polyres)      
  	end if 
  	end if   
  	end if      
  	      
  end do
  end do   
  


  delete(rest@cnLevels)
  delete(rest@cnFillColors)
  rest@cnLevels = lev1
  rest@cnFillColors = col1
  do j = 0, 1
  	rest@vpYF = 0.98 - 0 * (rest@vpHeightF + 0.03)
  	rest@vpXF = 0.04 + j * (rest@vpWidthF + 0.02)
    gsn_text_ndc( wks, label(j+2),rest@vpXF+0.01, rest@vpYF , restextl)
    gsn_text_ndc( wks, labelname(j),rest@vpXF+rest@vpWidthF/2., rest@vpYF + 0.02, restext)
  	plot(j+2) = gsn_csm_contour_map( wks,tas(j,:,:), rest)
   	plot1 = gsn_csm_contour( wks, tas_sig(j,:,:), rest1)
  	overlay(plot(j+2), plot1)
    respolyline@gsLineColor = "white"  
    x0 = (/lonmin1, lonmax1, lonmax1, lonmin1, lonmin1/)
    y0 = (/ latmax1, latmax1, latmin1, latmin1, latmax1/)
    polyline1(j) = gsn_add_polyline( wks, plot(j+2), x0, y0, respolyline)  
  end do  
  gsn_labelbar_ndc( wks, dimsizes(col), lab1, 0.2, rest@vpYF-rest@vpHeightF-0.01, reslb1)
  gsn_text_ndc( wks, "Surface temperature changes due to freshwater hosing (~S~o~N~C)", .5, rest@vpYF-rest@vpHeightF-0.05, restext)

  draw(wks)   
end
