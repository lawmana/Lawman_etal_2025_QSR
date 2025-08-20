;---Spatial Maps of Averaged variables based on Model IPSL, CESM02, CCSM

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
  dir = "/glade/u/home/xianwu/ncl/UTAustin/ncl/paleo/herinch/paper/data/seasonal/"
  dir2 = "/glade/u/home/xianwu/ncl/UTAustin/ncl/paleo/herinch/paper/data/"

  diro = "./"
  region = "Indian"
;  region = "global"
;  region = "SouthAmeria"

;  region = "Indian"
;  region = "global"
;  region = "SouthAmeria"
  ;region = "WesternAfrica"

if(region .eq. "Indian")then
  latmin = -40.
  latmax = 40.
  lonmin = 10.
  lonmax = 160.
end if

if(region .eq. "WesternAfrica")then
latmin = -20.
latmax = 40.
lonmin = -50.
lonmax = 50.
end if


;---TS & Rainfall

  f = addfile( dir + "Regridded_monthly_pr_tas.nc", "r")
  t = f->tas((/0,1,3/),:,:,:)
  r = f->pr((/0,1,3/),:,:,:)
  lat = f->lat
  lon = f->lon
  
;---UV
  f = addfile( dir + "Regridded_monthly_u_v.nc", "r")
  u = f->u((/0,1,3/),:,:,:)
  v = f->v((/0,1,3/),:,:,:)
  
;---SLP
  f = addfile( dir + "Regridded_monthly_ipsl.nc", "r")
  p1 = f->slp

  f = addfile( dir + "Regridded_monthly_cesm02.nc", "r")
  p2 = f->slp

  f = addfile( dir + "Regridded_monthly_ccsm.nc", "r")
  p3 = f->slp

;---iTraCE
  f = addfile( dir +"b.e13.Bi1850C5.f19_g16.17-16ka.itrace.ice_ghg_orb_wtr.01.cam.h0.clim.300001-499912.nc", "r")
  th = f->TS
  rh = f->PRECT
  uh = f->UBOT
  vh = f->VBOT
  ph = f->PSL
  
  f = addfile( dir +"b.e13.Bi1850C5.f19_g16.17-16ka.itrace.ice_ghg_orb.01.cam.h0.clim.300001-499912.nc", "r")
  tc = f->TS
  rc = f->PRECT
  uc = f->UBOT
  vc = f->VBOT
  pc = f->PSL

  t4 = th-tc
  r4 = rh-rc
  u4 = uh-uc
  v4 = vh-vc  
  p4 = ph-pc
  
;---4 model average

  tm = dim_sum_n_Wrap(t,0)
  rm = dim_sum_n_Wrap(r,0)
  um = dim_sum_n_Wrap(u,0)
  vm = dim_sum_n_Wrap(v,0)
  
  tm = (tm+t4)/4.
  rm = (rm+r4)/4.
  um = (um+u4)/4.
  vm = (vm+v4)/4.

  pm = p1
  pm = (p1+p2+p3+p4)/4.
  pm = pm/100.
  
;---3 month running mean
  tmm = runave_n_Wrap(tm, 3, -1, 0)
  rmm = runave_n_Wrap(rm, 3, -1, 0)
  umm = runave_n_Wrap(um, 3, -1, 0)
  vmm = runave_n_Wrap(vm, 3, -1, 0)
  pmm = runave_n_Wrap(pm, 3, -1, 0)
    
  printVarSummary(tmm)

  tmm_sel = tm(0:4,:,:)
  tmm_sel(1:4,:,:) = tmm(0::3,:,:)
  tmm_sel(0,:,:) = dim_avg_n_Wrap(tmm,0)

  rmm_sel = rm(0:4,:,:)
  rmm_sel(1:4,:,:) = rmm(0::3,:,:)
  rmm_sel(0,:,:) = dim_avg_n_Wrap(rmm,0)

  umm_sel = um(0:4,:,:)
  umm_sel(1:4,:,:) = umm(0::3,:,:)
  umm_sel(0,:,:) = dim_avg_n_Wrap(umm,0)

  vmm_sel = vm(0:4,:,:)
  vmm_sel(1:4,:,:) = vmm(0::3,:,:)
  vmm_sel(0,:,:) = dim_avg_n_Wrap(vmm,0)

  pmm_sel = pm(0:4,:,:)
  pmm_sel(1:4,:,:) = pmm(0::3,:,:)
  pmm_sel(0,:,:) = dim_avg_n_Wrap(pmm,0)

;---Proxy Data

  go = addfile(dir2 + "Regridded_precip_18models.nc", "r")
  ave_h = go->hose(model|:, lat|:, lon|:)
  ave_ci = ave_h(model|0, lat|:,lon|:)

  do i = 0, 95 
  do j = 0, 143
    ave_ci(i,j) = where(( all(ave_h(:,i,j).ge.0.) .or. all(ave_h(:,i,j).le.0.)), ave_h(0,i,j), ave_h@_FillValue)
  end do
  end do
  modelindex = (/5, 6, 9, 3, 4, 12/)
  ave_h := dim_avg_n_Wrap(ave_h(model|modelindex, lat|:, lon|:), 0)

  ave_c = go->control(model|:, lat|:, lon|:)
  ave_c := dim_avg_n_Wrap(ave_c, 0)

  p_ave = ave_h
  p_ave = ave_h/ave_c

  data = asciiread( dir2+"DiNezio_etal_HS1_database_2025_simplified.csv", (/157, 5/), "float")
  x = data(:,0)
  printVarSummary(data)

  printVarSummary(p_ave)
  do i = 0, 153
    print(i)
    print(data(i,1))
    print(data(i,2))
    x(i) = p_ave({lat|data(i,1)}, {lon|data(i,2)}) 
  end do

  x_flag = x 
  x_flag(:) = where(x.ge.(6/100.), 1., 0.)
  x_flag(:) = where(x.le.(-6/100.), -1., x_flag(:))
  x_flag(:) = where(abs(x).lt.(6/100.), 0., x_flag(:))


;---Plot
  
  wks  = gsn_open_wks( "pdf", get_script_prefix_name() +"-" + region)
  gsn_merge_colormaps(wks,"MPL_BrBG","amwg_blueyellowred")

;  cmap = gsn_retrieve_colormap( wks)
;  col = (/2,10,17,26,33,42,49,58,   74,82,90,98,106,113,122,129/)
;  lev = fspan( -2.1, 2.1, 15)
;  lab = sprintf( "%3.1f", lev)
;  lbcol = cmap(col,:)

;  col1 = (/2,3,4,5,6,7,8,9,    10,11,12,13,14,15,16,17/) + 128
;  lev1 = (/-2.0,-1.5,-1.0,-0.8,-0.6,-0.4,-0.2,0, 0.2,0.4,0.6,0.8,1.0,1.5,2.0/)
;  lab1 = sprintf( "%3.1f", lev1)
;  lbcol1 = cmap(col1,:)
 
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
;---TAS
  rest = True
  rest@mpMinLatF = latmin
  rest@mpMaxLatF = latmax
  rest@mpMinLonF = lonmin
  rest@mpMaxLonF = lonmax
  rest@gsnAddCyclic = False
  rest@gsnMajorLonSpacing = 30.
  rest@gsnMinorLonSpacing = 10.
  rest@gsnMajorLatSpacing = 30.
  rest@gsnMinorLatSpacing = 10.
  rest@tmYLMajorLengthF = 0.0025
  rest@tmYLMajorOutwardLengthF = 0.0025
  rest@tmXBMajorLengthF = 0.0025
  rest@tmXBMajorOutwardLengthF = 0.0025
  rest@tmYLMinorLengthF = 0.0025
  rest@tmYLMinorOutwardLengthF = 0.0025
  rest@tmXBMinorLengthF = 0.0025
  rest@tmXBMinorOutwardLengthF = 0.0025
  rest@tmXTOn = False
  rest@tmYROn = False
  rest@tmXBLabelFontHeightF = 0.01
  rest@tmYLLabelFontHeightF = 0.01
  rest@tmXBLabelJust = "TopCenter"
  rest@tmYLLabelJust = "CenterRight"
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
  rest@vpWidthF = 0.3
  rest@vpHeightF = rest@vpWidthF/2.

;---Wind
  resw = True
  resw@gsnDraw = False
  resw@gsnFrame = False
  resw@gsnLeftString = ""
  resw@gsnRightString = ""
  resw@gsnCenterString = ""
  resw@vcGlyphStyle = "FillArrow" 
  resw@vcFillArrowFillColor = "Black" 
  resw@vcFillArrowEdgeColor = "White" 
;  resw@vcMapDirection = False
  resw@vcRefLengthF = 0.015
  resw@vcRefAnnoOn = True 
  resw@vcRefAnnoSide = "Bottom" 
  resw@vcRefAnnoFontHeightF = 0.012
  resw@vcRefAnnoOrthogonalPosF = -0.14
  resw@vcRefAnnoString2On = False
  resw@vcRefMagnitudeF = 1.
  
;---SLP
  resp = True
  resp@gsnDraw = False
  resp@gsnFrame = False
  resp@gsnLeftString = ""
  resp@gsnCenterString = ""
  resp@gsnRightString = ""
  resp@cnLineLabelsOn = False
  resp@cnInfoLabelOn = False
  resp@cnLevelSpacingF = 0.5
  
;---Text
  restext = True
  restext@txJust = "TopCenter"
  restext@txFontHeightF = 0.011


;---Text1
  restx = True
  restx@txPerimOn = True
  restx@txBackgroundFillColor = "White"
  restx@txFontHeightF = 0.011
;  restx@txConstantSpacingF = 1.05
  restx@txJust = "CenterLeft"
  
;---Text2
  restext2 = True
  restext2@txJust = "TopLeft"
  restext2@txFontHeightF = 0.011
  restext2@txFont = "helvetica-bold"

;---Textlabel
restextl = True
restextl@txJust = "TopLeft"
restextl@txFontHeightF = 0.011
restextl@txFont   = "helvetica-bold"

;---Color Bar
  reslb = True
  reslb@lbPerimOn = False
  reslb@lbOrientation = "Horizontal"
  reslb@lbLabelFontHeightF = 0.009
  reslb@lbLabelAlignment = "InteriorEdges"
  reslb@lbMonoFillPattern = True
  reslb@lbFillColors = lbcol
  reslb@vpWidthF = 0.45
  reslb@vpHeightF = 0.04
  
;---Color Bar
  reslb1 = reslb
  delete(reslb1@lbFillColors)
  reslb1@lbFillColors = lbcol1

  
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

;---Polyline
  polyresl = True
  polyresl@gsLineColor = "black" 
  polyresl@gsLineThicknessF = 1.  

  plot = new( 10, "graphic")
  plotl = new(10,"graphic")
  plotl2 = new(10,"graphic")

  gsn_text_ndc( wks, "Average for IPSL-CM5, CCSM3-1Sv-NA, CESM1-0.2Sv, and CESM1-iTraCE-MWF", 0.5, 1.0 , restext)
  gsn_text_ndc( wks, "Rainfall/Surface Wind", 0.17 + rest@vpWidthF/2., 0.985 , restext)
  gsn_text_ndc( wks, "Surface Temperature/SLP", 0.53 + rest@vpWidthF/2., 0.985, restext)

  restext@txAngleF = 0
  restext@txFontHeightF = 0.014
  labeli = (/"Annual","DJF","MAM","JJA","SON"/)
  label = (/"a","b","c","d"/) 

  
  hnum = (/3,1/)
;---Rainfall, Wind
  do h=0,1
    hs = hnum(h)
    rest@vpXF = 0.17
    rest@vpYF = 0.97 - h * (rest@vpHeightF + 0.03)
    plot(h) = gsn_csm_contour_map( wks, rmm_sel(hs,:,:), rest)
    plotw = gsn_csm_vector( wks, umm_sel(hs,::2,::3), vmm_sel(hs,::2,::3), resw)
    overlay(plot(h), plotw)
    gsn_text_ndc( wks, label(h),rest@vpXF-.015, rest@vpYF + 0.01, restextl)
    text = gsn_add_text(wks,plot(h),labeli(hs),lon({lonmin+1}),lat({latmax-5.9}),restx)
    plotl(h) = gsn_add_polyline( wks, plot(h), (/0.,-179.9/), (/0.,0./), polyresl)
    plotl2(h) = gsn_add_polyline( wks, plot(h), (/0.,179.9/), (/0.,0./), polyresl)
  end do
  gsn_labelbar_ndc( wks, dimsizes(col), lab, 0.15 , rest@vpYF-rest@vpHeightF-0.024, reslb)
  gsn_text_ndc( wks, "Precipitation changes (mm/day)", .605, rest@vpYF-rest@vpHeightF-0.024-0.013, restext2)
  delete(rest@cnLevels)
  delete(rest@cnFillColors)
  rest@cnLevels = lev1
  rest@cnFillColors = col1

;---TS, SLP
  do h =0,1
    rest@vpXF = 0.53
    rest@vpYF = 0.97 - h * (rest@vpHeightF + 0.03)
    hs = hnum(h)
    plot(h+5) = gsn_csm_contour_map( wks, tmm_sel(hs,:,:), rest)
    
    resp@cnLineThicknessF = 3.5
    resp@cnLineColor = "White"  
    plotp = gsn_csm_contour( wks, pmm_sel(hs,:,:), resp)
    overlay( plot(h+5), plotp)

    resp@cnLineThicknessF = 1.5
    resp@gsnContourNegLineDashPattern = 1
    resp@cnLineColor = "Black"  
    plotp = gsn_csm_contour( wks, pmm_sel(hs,:,:), resp)
    gsn_text_ndc( wks, label(h+2),rest@vpXF-.015, rest@vpYF + 0.01, restextl)
    text = gsn_add_text(wks,plot(h+5),labeli(hs),lon({lonmin+1}),lat({latmax-5.9}),restx)
    overlay(plot(h+5), plotp)
    plotl(h+5) = gsn_add_polyline( wks, plot(h+5), (/0.,-179.9/), (/0.,0./), polyresl)
    plotl2(h+5) = gsn_add_polyline( wks, plot(h+5), (/0.,179.9/), (/0.,0./), polyresl)

  end do
  reslb1@vpWidthF = 0.45
  gsn_labelbar_ndc( wks, dimsizes(col1), lab1, 0.15, rest@vpYF-rest@vpHeightF-0.024*2-0.01, reslb1)
  gsn_text_ndc( wks, "Surface temperature changes (~S~o~N~C)", .605, rest@vpYF-rest@vpHeightF-0.024*2-0.017, restext2)

  lat = umm_sel&lat
  print(lat(0::4))

 
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


  draw( wks)
  frame( wks)
  delete( plot)
end



