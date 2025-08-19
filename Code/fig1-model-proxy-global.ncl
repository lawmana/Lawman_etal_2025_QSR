
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  dir = "/glade/u/home/xianwu/ncl/UTAustin/ncl/paleo/herinch/paper/data/"

  modelindex = (/3,4,5,6,7,10/)
  ;modelindex = (/3,4,5,6,7/)

  modelnum = dimsizes(modelindex)
  makedata = True
  
;------------------------------------------------------------------
  print ("Reading data files...")
;    

if(makedata) then
 
  region = (/  "Australia_Rainfall",\                 ;4
              "North_Eq_Indian_Ocean_Rainfall", \  ;5
              "East_Asia_Rainfall",\  ;6
              "Maritime_Continent_Rainfall", \  ;7              
              "Caribbean_Rainfall",\  ;10
              "Brazil_Rainfall",\  ;11
              "Western_US_Rainfall",\  ;12
              "Northern_Andes_Rainfall",\  ;13
              "Southern_Andes_Rainfall"/)       
              
  location = (/ (/-25., -15., 110., 150./),\
                (/0.,30.,60.,100./),\  
                (/22.,35.,95.,120./),\ 
                (/-5.,10.,95.,125./),\    
                (/0.,28.,-100.,-55./),\
                (/-18.,3.,-55.,-30./),\
                (/35.,45.,-125.,-110./),\
                (/-15.,0.,-80.,-75./),\
                (/-25.,-15.,-70.,-65./)/)

   
  go = addfile(dir + "Regridded_precip_18models.nc", "r")
  ave_h = go->hose(model|:, lat|:, lon|:)
  ave_ci = ave_h(model|0, lat|:,lon|:)
  
  do i = 0, 95 
  do j = 0, 143
  	ave_ci(i,j) = where(( all(ave_h(:,i,j).ge.0.) .or. all(ave_h(:,i,j).le.0.)), ave_h(0,i,j), ave_h@_FillValue)
  end do
  end do
  
  ave_h := dim_avg_n_Wrap(ave_h(model|modelindex, lat|:, lon|:), 0)
  
  ave_c = go->control(model|:, lat|:, lon|:)
  ave_c := dim_avg_n_Wrap(ave_c, 0)
  
  p_ave = ave_h
  p_ave = ave_h/ave_c
   
  ;data = asciiread( dir+"HS1database_2023_simplified.csv", (/142, 6/), "float")
  data = asciiread( dir+"DiNezio_etal_HS1_database_2025_simplified.csv", (/157, 5/), "float")
  x = data(:,0)
  printVarSummary(data)

  
  printVarSummary(p_ave)
  do i = 0, 153
    ;print(i)
    ;print(data(i,1))
    ;print(data(i,2))
  	x(i) = p_ave({lat|data(i,1)}, {lon|data(i,2)}) 
  end do
  
  x_flag = x 
  x_flag(:) = where(x.ge.(6/100.), 1., 0.)
  x_flag(:) = where(x.le.(-6/100.), -1., x_flag(:))
  x_flag(:) = where(abs(x).lt.(6/100.), 0., x_flag(:))
 

;---Figures 
;
  wks  = gsn_open_wks( "pdf", get_script_prefix_name()+"-"+modelnum)
  cgrn = "Forest Green"
  cbrw = "saddlebrown"  
  cgry = "Gray"
;  gsn_merge_colormaps(wks,"amwg_blueyellowred",(/cgrn,cbrw,cgry/));
;
;  lev = fspan( -0.35, 0.35, 15) 
;  col = ispan( 2, 17, 1) 

;  gsn_merge_colormaps(wks,"MPL_BrBG",(/cgrn,cbrw,cgry/))
;  lev = (/-1.5, -.8, -0.5, -0.3,-0.1,  0.,0.1,0.3,0.5,.8,1.5/)
;  col = (/4,  18, 30,  35, 42,  57,    73, 83, 98, 109, 117, 127/)

  gsn_merge_colormaps(wks,"MPL_BrBG","amwg_blueyellowred")

  cmap = gsn_retrieve_colormap( wks)
  
;  lev = (/-1.3,-1.1,-0.9, -.7, -0.5, -0.3,-0.1,  0.,0.1,0.3,0.5,.7,0.9,1.1,1.3/)

  lev = fspan(-1.4,1.4,15)
  col = (/2,10,17,26,33,42,49,58,   74,82,90,98,106,113,122,129/)


;  lev = (/0.,1./)
;  col = (/0,0,0/)
   
  lab = sprintf( "%3.1f", lev)  
;  lab(0::2) = "" 

  cmap = gsn_retrieve_colormap( wks)
  lbcol = cmap(col,:)

;---SST Anomaly
  rest = True 
  rest@gsnDraw = False
  rest@gsnFrame = False 
  rest@gsnRightString = ""
  rest@gsnLeftString = ""
  rest@gsnCenterString = "" 
  rest@gsnAddCyclic = True 
;  rest@mpProjection = "CylindricalEqualArea"
  rest@mpCenterLonF = 0.   
  rest@mpPerimOn = False      
;  rest@mpOutlineOn = False 
  rest@mpGridAndLimbOn = False 
  rest@mpGridLineColor = "transparent"
  rest@mpGeophysicalLineColor = "black"   
  rest@mpLandFillColor = "transparent"
  rest@mpFillDrawOrder = "PreDraw"   
  rest@cnFillOn = True
  rest@cnFillOpacityF = 0.8  
  rest@cnLinesOn = False
  rest@cnLevelSelectionMode = "ExplicitLevels"
  rest@cnLevels = lev
  rest@cnFillColors = col
  rest@cnFillDrawOrder = "PreDraw"
  rest@cnLineLabelsOn = False 
  rest@cnInfoLabelOn = False
  rest@lbLabelBarOn = False  
  rest@tmYLLabelFontHeightF = 0.012
  rest@tmXBLabelFontHeightF = 0.012
  rest@tmXTOn = False 
  rest@tmYROn = False
  rest@tmXBOn = True
  rest@tmYLOn = True  
      
;  rest@pmTickMarkDisplayMode = "Never"    
  rest@mpLimitMode = "LatLon"
  rest@mpMinLatF = -40.  
  rest@mpMaxLatF = 40.   
  rest@mpMinLonF = -180.      
  rest@mpMaxLonF = 180.   
     
 
;---Text
  restext = True
  restext@txJust = "TopCenter"  
  restext@txFontHeightF = 0.011  

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

;---attach 
  attachres1                     = True 
;  attachres1@gsnAttachPlotsYAxis = True  ;; attaches along x-axis
  attachres1@gsnAttachBorderOn   = False ;; No border please
  attachres2                     = True 
;  attachres2@gsnAttachPlotsYAxis = True  ;; attaches along x-axis
  attachres2@gsnAttachBorderOn   = False ;; No border please
   
  
  rest@vpWidthF = 0.9 
  rest@vpHeightF = 0.3   
  rest@vpYF = 0.98 
  rest@vpXF = 0.05
  
  labeli = (/"CCSM_MARUM", "CCSM_NCAR", "HadCM3_S", "HadCM3_W","MIROC_S","MIROC_W","COSMOS_S","COSMOS_W","IPSL","CESM(PDN)","Ave_wTAc"/)
  plot = new(10, "graphic")
  
  rest@vpYF = 0.98 - 0 * (rest@vpHeightF + 0.03)
  rest@vpXF = 0.05 + 0 * (rest@vpWidthF + 0.01)
  plot(0) = gsn_csm_contour_map( wks, ave_h, rest)
 
  plotl1 = gsn_add_polyline( wks, plot(0), (/0.,-179.9/), (/0.,0./), polyresl)
  plotl2 = gsn_add_polyline( wks, plot(0), (/0.,179.9/), (/0.,0./), polyresl)
  
 
  plotd = new(1500, "graphic")        
  plotd1 = new(1500, "graphic")    
  plotd2 = new(1500, "graphic")      
  
  do k = 0, 0  
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
  
;  rest@vpYF = 0.98 - 1 * (rest@vpHeightF + 0.03)
;  rest@vpXF = 0.05 + 0 * (rest@vpWidthF + 0.01)
;    
;  plot(1) = gsn_csm_contour_map( wks, ave_h, rest)
   
;  do i = 0, 73
;  	dum = gsn_add_text( wks, plot(1), (/data(i*2, 0)/), data(i*2, 2), data(i*2, 1), restext)
;  end do

;  rest@vpYF = 0.98 - 2 * (rest@vpHeightF + 0.03)
;  rest@vpXF = 0.05 + 0 * (rest@vpWidthF + 0.01)
    
;  plot(2) = gsn_csm_contour_map( wks, ave_h, rest)
   
;  do i = 0, 72
;  	dum = gsn_add_text( wks, plot(2), (/data(i*2+1, 0)/), data(i*2+1, 2), data(i*2+1, 1), restext)
;  end do
  
  
;  do kk = 0, 8 
;  
;  	px1 = (/location(kk,2), location(kk,3), location(kk,3), location(kk,2), location(kk,2)/)
;  	py1 = (/ location(kk,1), location(kk,1), location(kk,0), location(kk,0), location(kk,1)/)
;
;    plotl(kk) = gsn_add_polyline( wks, plot(1), px1, py1, polyresl)
;    plotl1(kk) = gsn_add_polyline( wks, plot(0), px1, py1, polyresl)
;    
;  end do 
  
  gsn_labelbar_ndc( wks, dimsizes(col), lab, 0.2, 0.69, reslb) 
  restext@txFont = "helvetica-bold"
  gsn_text_ndc( wks, "Precipitation change due to freshwater hosing (mm/day)", .5, 0.65, restext)


  draw( wks) 
  frame( wks)    
  end if

  printVarSummary(ave_h)
end
