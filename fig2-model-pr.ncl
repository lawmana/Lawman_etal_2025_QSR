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
  "CCSM3-0.2Sv-GIN-LR", "IPSL-CM4", "CCSM3-TraCE-MWF", "CESM1-iTraCE-MWF",\
  "IPSL-CM5", "CCSM3-0.2Sv-GIN-HR", "CESM1-0.15Sv",\
  "CCSM4-0.1Sv-RB", "COSMOS-W", "COSMOS-S", "CCSM3-0.1Sv-NA", "MIROC-W", \
  "HadCM3-0.1Sv", "CESM1-0.1Sv", "MIROC-S","Mean of 18 Model Experiments"/)

  g = addfile(dir + "Regridded_precip_18models.nc", "r")
  hose = g->hose
  control = g->control
  if (any(isnan_ieee(hose))) then
      value = 1.e20
      replace_ieeenan (hose, value, 0)
      hose@_FillValue = value
  end if
  if (any(isnan_ieee(control))) then
      value = 1.e20
      replace_ieeenan (control, value, 0)
      control@_FillValue = value
  end if
  printVarSummary(hose)
  hose_avg = dim_avg_n_Wrap(hose,0)
  printVarSummary(hose_avg)
  hose := array_append_record(hose,conform(hose(0:0,:,:),hose_avg,(/1,2/)),0)
  printVarSummary(hose)
  hosenew = hose
  ;hosenew(0,:,:) = hose(17,:,:)
  ;hosenew(1:17,:,:) = hose(0:16,:,:)

  print(hose_avg)

;---Figures
; 
  wks  = gsn_open_wks( "pdf", get_script_prefix_name() )
  gsn_merge_colormaps(wks,"MPL_BrBG","amwg_blueyellowred")

  cmap = gsn_retrieve_colormap( wks)
  col = (/2,10,17,26,33,42,49,58,   74,82,90,98,106,113,122,129/)
  lev = fspan( -1.4, 1.4, 15)
  lab = sprintf( "%3.1f", lev)
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
  rest@mpGeophysicalLineColor = "black"
  rest@mpLandFillColor = "transparent"
  rest@cnFillOn = True
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
      
;  rest@pmTickMarkDisplayMode = "Never"    
;  rest@mpLimitMode = "LatLon"
;  rest@mpMinLatF = -40. 
;  rest@mpMaxLatF = 40. 
;  rest@mpMinLonF = -90.
;  rest@mpMaxLonF = 50.
  


;---Text
  restext = True
  restext@txJust = "TopCenter"
  restext@txFontHeightF = 0.01

;---Color Bar
  reslb = True
  reslb@lbPerimOn = False
  reslb@lbOrientation = "Horizontal"
  reslb@lbLabelFontHeightF = 0.01
  reslb@lbLabelAlignment = "InteriorEdges"
  reslb@lbMonoFillPattern = True
  reslb@lbFillColors = lbcol
  reslb@vpWidthF = 0.5
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

  rest@vpWidthF = 0.2
  rest@vpHeightF = 0.1
  rest@vpYF = 0.98 
  rest@vpXF = 0.18
  
  labeli = model_names  
  plot = new(20, "graphic") 
  
  do i = 0, 6
    if(i .le. 5) then
      do j = 0, 2
    	rest@vpYF = 0.98 - i * (rest@vpHeightF + 0.02)
    	rest@vpXF = 0.18 + j * (rest@vpWidthF + 0.01)
  	    gsn_text_ndc( wks, labeli(i*3+j), rest@vpXF+rest@vpWidthF/2., rest@vpYF+0.014, restext)
    	plot(i*3+j) = gsn_csm_contour_map( wks, hosenew(i*3+j,:,:), rest)
      end do
  	else
      do j = 0, 0
    	rest@vpYF = 0.98 - i * (rest@vpHeightF + 0.02)
    	rest@vpXF = 0.18 + (j+1) * (rest@vpWidthF + 0.01)
 
      restext@txFont = "helvetica-bold"

      gsn_text_ndc( wks, labeli(i*3+j), rest@vpXF+rest@vpWidthF/2., rest@vpYF+0.014, restext)
    	plot(i*3+j) = gsn_csm_contour_map( wks, hosenew(i*3+j,:,:), rest)
      end do
  	end if
  end do
  gsn_labelbar_ndc( wks, dimsizes(col), lab, 0.25, rest@vpYF-rest@vpHeightF-0.005, reslb)
  restext@txFont = "helvetica-bold"

  gsn_text_ndc( wks, "Precipitation change due to freshwater hosing (mm/day)", .5, rest@vpYF-rest@vpHeightF-0.045, restext)

  draw( wks)
  frame( wks)

end
