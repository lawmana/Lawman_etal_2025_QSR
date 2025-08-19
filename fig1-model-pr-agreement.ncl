load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  dir = "/glade/u/home/xianwu/ncl/UTAustin/ncl/paleo/herinch/paper/data/"

;------------------------------------------------------------------
  print ("Reading data files...")
;    
  model_names = (/"Mean of 18 Model Experiments","HadCM3-0.4Sv", "CCSM3-1Sv-NA", "CESM1-0.2Sv",\
                  "CCSM3-0.2Sv-GIN-LR", "IPSL-CM4", "CCSM3-TraCE-MWF", "CESM1-iTraCE-MWF","IPSL-CM5", "CCSM3-0.2Sv-GIN-HR", "CESM1-0.15Sv",\
                  "CCSM4-0.1Sv-RB", "COSMOS-W", "COSMOS-S", "CCSM3-0.1Sv-NA", "MIROC-W", "HadCM3-0.1Sv", "CESM1-0.1Sv", "MIROC-S"/)
	 
  g = addfile(dir + "Regridded_precip_18models.nc", "r")
  hose = g->hose
  control = g->control
  printVarSummary(hose)

  control = where(control.eq.0, 1e-10, control)  ; Replace zeros with a small value to avoid division by zero

  hose_avg = dim_avg_n_Wrap(hose,0)
  control_avg = dim_avg_n_Wrap(control, 0)
  control_avg = where(control_avg.eq.0, 1e-10, control_avg)  ; Same for control_avg

  ratio = hose
  ratio = hose/control

  ratio_avg = hose_avg
  ratio_avg = hose_avg/control_avg

  hose_num = hose_avg
  
  printVarSummary(hose_num)
  
  threshold = 0.05
  do i = 0, 95 
  do j = 0, 143
    if(ratio_avg(i,j) .gt. threshold)then
      hose_num(i,j) = num(ratio(:,i,j).gt.threshold)
    else if(ratio_avg(i,j) .lt. threshold)then
      hose_num(i,j) = num(ratio(:,i,j).lt.-threshold)
    else
      hose_num(i,j) = num(ratio(:,i,j).ge.-threshold .and. ratio(:,i,j).le.threshold )
  	;hose_num(i,j) = where(hose_avg(i,j)/.gt.0, num(hose(:,i,j).gt.0), num(hose(:,i,j).lt.0))
    end if
    end if
  end do
  end do
  printVarSummary(hose_num)
  ;hose_num := lonPivot(hose_num, 1)
  ;hose_num := linint2_Wrap(hose_num&lon, hose_num&lat, hose_num, True, lonGlobeF(288, "lon", "longitude", "degrees_east"), hose_num&lat, 0)
;
;---Figures
; 
  wks  = gsn_open_wks( "pdf", get_script_prefix_name() )
  ;gsn_merge_colormaps(wks,"WhiteBlue","amwg_blueyellowred")
  ;col = (/15,45,75,105,135,165,195,225,255/)
  gsn_define_colormap(wks, "CBR_wet")
  ;gsn_define_colormap(wks, "CBR_drywet")
  cmap = gsn_retrieve_colormap( wks)
  ;col = (/2,4,6,7,8,9,10,11,12/)
  col = (/12,11,10,9,8,7,6,4,2/)
  ;lev = fspan( 9,16,8)
  lev = fspan(2,16, 8)
  lab = sprintf( "%3.0f", lev)
  
  ;col = ispan(2,12,1)
  ;lev = (/-0.45,-0.35,-0.25,-0.15,-0.05,0.05,0.15,0.25,0.35,0.45/)
  ;lab = sprintf( "%3.2f", lev)
  
  
  lbcol = cmap(col,:)

;---SST Anomaly
  rest = True
  rest@gsnDraw = False
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
  rest@mpGeophysicalLineColor = "Black"
  rest@mpGeophysicalLineThicknessF = 2.
  rest@mpLandFillColor = "transparent"
  rest@cnFillOn = True
  rest@cnFillOpacityF = 1.
  rest@cnLinesOn = False
  rest@cnLevelSelectionMode = "ExplicitLevels"
  rest@cnLevels = lev
  rest@cnFillColors = col
  ;rest@cnFillDrawOrder = "PreDraw"
  rest@cnLineLabelsOn = False
  rest@cnInfoLabelOn = False
  rest@lbLabelBarOn = False
  rest@tmYLLabelFontHeightF = 0.01
  rest@tmXBLabelFontHeightF = 0.01
  rest@tmXTOn = False
  rest@tmYROn = False 
  rest@tmXBOn = False
  rest@tmYLOn = False
      
;  rest@pmTickMarkDisplayMode = "Never"    
;  rest@mpLimitMode = "LatLon"
;  rest@mpMinLatF = -40. 
;  rest@mpMaxLatF = 40. 
;  rest@mpMinLonF = -90.
;  rest@mpMaxLonF = 50.
  


;---Text
  restext = True
  restext@txJust = "TopCenter"
  restext@txFontHeightF = 0.011

;---Color Bar
  reslb = True
  reslb@lbPerimOn = False
  reslb@lbOrientation = "Horizontal"
  reslb@lbLabelFontHeightF = 0.01
  reslb@lbLabelAlignment = "InteriorEdges"
  reslb@lbMonoFillPattern = True
  reslb@lbFillColors = lbcol
  reslb@vpWidthF = 0.4
  reslb@vpHeightF = 0.04
   
;---polygon--- 
  
  polyres                   = True
  polyres@gsMarkerIndex     = 16          ; polymarker style
  polyres@gsMarkerSizeF     = 6          ; polymarker size

;---ploycover---
  polyres1                   = True
  polyres1@gsMarkerIndex     = 4          ; polymarker style
  polyres1@gsMarkerSizeF     = 6          ; polymarker size
  polyres1@gsMarkerColor = "black"  	

  rest@vpWidthF = 0.6
  rest@vpHeightF = 0.3
  rest@vpYF = 0.98 
  rest@vpXF = 0.1
  
  labeli = model_names  
  plot = new(20, "graphic") 
  

    	rest@vpYF = 0.98 
    	rest@vpXF = 0.2
  
      restext@txFont = "helvetica-bold"
  	  ;gsn_text_ndc( wks, labeli(i*3+j), rest@vpXF+rest@vpWidthF/2., rest@vpYF+0.014, restext)
    	plot(0) = gsn_csm_contour_map( wks, hose_num, rest)
    
  gsn_labelbar_ndc( wks, dimsizes(col), lab, 0.3, rest@vpYF-rest@vpHeightF-0.025, reslb)
  restext@txFont = "helvetica-bold"
  gsn_text_ndc( wks, "Number of simulations with the same sign as the multi-model mean", .5, rest@vpYF-rest@vpHeightF-0.065, restext)
  gsn_text_ndc( wks, "Weaker agreement", .3, rest@vpYF-rest@vpHeightF-0.01, restext)
  gsn_text_ndc( wks, "Stronger agreement", .7, rest@vpYF-rest@vpHeightF-0.01, restext)
  draw( wks)
  frame( wks)

end
