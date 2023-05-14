;Last update: May 14, 2023: Adapt for CIRES job interview 
;Last update: Aug 11, 2015
pro save_netcdf,scheme,site,num_time,num_lev,$
                sprof_lev_site,temp_site,dens_site,pp_site,ghtp_site,uu_site,vv_site,ww_site, $
                qc_site,qi_site,qv_site,qr_site,qs_site,qg_site, $
                xnor_site,xnog_site,dbz_site,smo4_site,qnrv_site,res_site,lamr_site,lamg_site, $
                xDs_site,xDs_m_site,rhos_site,rvt_new_site

filenameout='saveddata/'+scheme+'_'+site+'_TH_20051231_5min_dbz_qn_res_rhos_newRvt_BB.nc'
cdfid=ncdf_create(filenameout,/clobber)
xdimid=ncdf_dimdef(cdfid,'Time',num_time)
ydimid=ncdf_dimdef(cdfid,'Height',num_lev)
ncdf_attput, cdfid, 'creation_date', systime(), /global

vid_sprof_lev=ncdf_vardef(cdfid,'sprof_lev',[ydimid],/float)
ncdf_attput, cdfid, vid_sprof_lev, 'description',site+' S-profiler vertical levels (0.25 km interval)'
vid_temp=ncdf_vardef(cdfid,'temp',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_temp, 'description','Temperature (K)'
vid_dens=ncdf_vardef(cdfid,'dens',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_dens, 'description','Air density (kg/m^3)'
vid_pp=ncdf_vardef(cdfid,'pp',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_pp, 'description','Pressure (mb)'
vid_ghtp=ncdf_vardef(cdfid,'ghtp',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_ghtp, 'description','Geopotential height (km)'
vid_uu=ncdf_vardef(cdfid,'uu',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_uu, 'description','x-component wind (m/s)'
vid_vv=ncdf_vardef(cdfid,'vv',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_vv, 'description','y-component wind (m/s)'
vid_ww=ncdf_vardef(cdfid,'ww',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_ww, 'description','z-component wind (m/s)'
vid_qc=ncdf_vardef(cdfid,'qc',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_qc, 'description','Cloud water mixing ratio (kg/kg)'
vid_qi=ncdf_vardef(cdfid,'qi',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_qi, 'description','Cloud ice mixing ratio (kg/kg)'
vid_qv=ncdf_vardef(cdfid,'qv',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_qv, 'description','Water vapor mixing ratio (kg/kg)'
vid_qr=ncdf_vardef(cdfid,'qr',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_qr, 'description','Rain mixing ratio (kg/kg)'
vid_qs=ncdf_vardef(cdfid,'qs',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_qs, 'description','Snow mixing ratio (kg/kg)'
vid_qg=ncdf_vardef(cdfid,'qg',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_qg, 'description','Graupel mixing ratio (kg/kg)'
vid_xnor=ncdf_vardef(cdfid,'xnor',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_xnor, 'description','Intercept parameter, rain (m^-4)'
vid_xnog=ncdf_vardef(cdfid,'xnog',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_xnog, 'description','Intercept parameter, graupel (m^-4)'
vid_dbz=ncdf_vardef(cdfid,'dbz',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_dbz, 'description','Reflectivity (Rayleigh, dBZ)'
vid_smo4=ncdf_vardef(cdfid,'smo4',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_smo4, 'description','4th moment of snow particle diameter (m^4/m^3)'
vid_rvt_new=ncdf_vardef(cdfid,'rvt_new',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_rvt_new, 'description','New Reflectivity-weighted terminal velocity,dielectric factor and new combination method, total (m/s)'
vid_qnrv=ncdf_vardef(cdfid,'qnrv',[xdimid,ydimid],/float) ;use the same name in Thompson scheme as in Morrison scheme 
ncdf_attput, cdfid, vid_qnrv, 'description','Number concentration, rain, (m^-3)'
vid_res=ncdf_vardef(cdfid,'res',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_res, 'description','Effective radius, re, snow (micron)'
vid_lamr=ncdf_vardef(cdfid,'lamr',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_lamr, 'description','lamda, rain (1/m)'
vid_lamg=ncdf_vardef(cdfid,'lamg',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_lamg, 'description','lamda, graupel (1/m)'
vid_xDs=ncdf_vardef(cdfid,'xDs',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_xDs, 'description','Effective diameter M3/M2, xDs, snow (m)'
vid_xDs_m=ncdf_vardef(cdfid,'xDs_m',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_xDs_m, 'description','Effective diameter M3/M2, xDs_m, snow (m), greater than 200 micron'
vid_rhos=ncdf_vardef(cdfid,'rhos',[xdimid,ydimid],/float)
ncdf_attput, cdfid, vid_rhos, 'description','Variable Snow Density, rhos (kg/m^3)'

ncdf_control, cdfid, /endef

ncdf_varput, cdfid, vid_sprof_lev, sprof_lev_site
ncdf_varput, cdfid, vid_temp, temp_site
ncdf_varput, cdfid, vid_dens, dens_site
ncdf_varput, cdfid, vid_pp, pp_site
ncdf_varput, cdfid, vid_ghtp, ghtp_site
ncdf_varput, cdfid, vid_uu, uu_site
ncdf_varput, cdfid, vid_vv, vv_site
ncdf_varput, cdfid, vid_ww, ww_site
ncdf_varput, cdfid, vid_qc, qc_site
ncdf_varput, cdfid, vid_qi, qi_site
ncdf_varput, cdfid, vid_qv, qv_site
ncdf_varput, cdfid, vid_qr, qr_site
ncdf_varput, cdfid, vid_qs, qs_site
ncdf_varput, cdfid, vid_qg, qg_site
ncdf_varput, cdfid, vid_xnor, xnor_site
ncdf_varput, cdfid, vid_xnog, xnog_site
ncdf_varput, cdfid, vid_dbz, dbz_site
ncdf_varput, cdfid, vid_smo4, smo4_site
ncdf_varput, cdfid, vid_rvt_new, rvt_new_site
ncdf_varput, cdfid, vid_qnrv, qnrv_site
ncdf_varput, cdfid, vid_res, res_site
ncdf_varput, cdfid, vid_lamr, lamr_site
ncdf_varput, cdfid, vid_lamg, lamg_site
ncdf_varput, cdfid, vid_xDs, xDs_site
ncdf_varput, cdfid, vid_xDs_m, xDs_m_site
ncdf_varput, cdfid, vid_rhos, rhos_site

ncdf_close, cdfid

end ;end of pro save_netcdf

pro comp_th_data2prof_netcdf_Thom_BB,dbz,dbz_r,dbz_s,dbz_g
;Aug 11, 2015: This pro is modified from "comp_th_data2prof_netcdf_Thom" that was used for Han et al. 2013)
;In this pro, I used the updated reflectivity calculation to account for the bright band.
;It is for the paper to be submitted in 2015, comparing THOM with HUCM SBM.

;Give the name of scheme:
;schemename=strarr(2)
;schemename=['gsfc','wsm6']

schemename=strarr(1)
;schemename=['wsm6']
;schemename=['gsfc']
;schemename=['morr']
schemename=['thom']

;for ischeme=0,1 do begin
for ischeme=0,0 do begin
  ;Use "for-loop" to go through the time steps:
  ;Calculate for 5-min wrfout stored in $ARCHIVE:
  num_time=289 ;5-min data
  openr,10,'/nfs3m/archive/sfa_cache05/users/g13/mhan/HMT_WRFV3.1/case_thom_rst5min/wrfoutd4/filenamelist'
  path='/nfs3m/archive/sfa_cache05/users/g13/mhan/HMT_WRFV3.1/case_thom_rst5min/wrfoutd4/'

  for i=0,num_time-1 do begin
    filename=''
    readf, 10, filename
    filename=path+filename
    print,filename
    readncfile, filename,scheme=schemename[ischeme],times,dx,lon,lat,ter,ghtp,ght,rho,temp,pp,nx,ny,nz,nt,$
                qv,qc,qi,qr,qs,qg,uu,vv,ww,qni=qni,qnr=qnr,qns=qns,qng=qng,swe,snowh

    print,times
  
    case 1 of
      ;(schemename[ischeme] eq 'gsfc'): dbz_vt_calc_gsfc, qr,qs,qg,pp,temp,qv,nx,ny,nz,nt,dbz,vt,vtr,vts,vtg,r_vt,r_vtr,r_vts,r_vtg,rcp_lamr,rcp_lams,rcp_lamg,dens
      ;(schemename[ischeme] eq 'wsm6'): dbz_vt_calc_wsm6, qr,qs,qg,pp,temp,qv,nx,ny,nz,nt,dbz,vt,vtr,vts,vtg,r_vt,r_vtr,r_vts,r_vtg,rcp_lamr,rcp_lams,rcp_lamg,dens,xnos
      ;(schemename[ischeme] eq 'gsfc'): dbz_vt_calc_gsfc_newRvt, qr,qs,qg,pp,temp,qv,nx,ny,nz,nt,dbz,vt,vtr,vts,vtg,r_vt,r_vtr,r_vts,r_vtg,rcp_lamr,rcp_lams,rcp_lamg,dens,r_vt_new
      ;(schemename[ischeme] eq 'wsm6'): dbz_vt_calc_wsm6_newRvt, qr,qs,qg,pp,temp,qv,nx,ny,nz,nt,dbz,vt,vtr,vts,vtg,r_vt,r_vtr,r_vts,r_vtg,rcp_lamr,rcp_lams,rcp_lamg,dens,xnos,r_vt_new
      ;(schemename[ischeme] eq 'morr'): dbz_calc_morr, qr,qs,qg,qnr,qns,qng,pp,temp,qv,nx,ny,nz,nt,dbz,dbzl,n0r,n0s,n0g,xnor,xnos,xnog,lamr,lams,lamg,dens
                                ;In Morrison scheme, n0r, n0s, and n0g are intercepts with unit of /(Kg m). 
                                ; n0r * dens (air density) = xnor in the unit of /(m^4).
      ;(schemename[ischeme] eq 'thom'): dbz_calc_thom, qr,qs,qg,pp,temp,qv,nx,ny,nz,nt,dbz,rcp_lamr,rcp_lamg,dens, $
      ;                                qnr,nr,smo4,res,xnor,xnog,lamr,lamg,dbz_r,dbz_s,dbz_g
                                ;In Thompson scheme, qnr is number concentration per unit mass (unit is /Kg) 
                                ;                    nr (=qnr*dens) is number concentration per unit volumn (unit is /(m^3))
                                ;smo4 is the forth moment of snow  ...... proportion to radar reflectivity in Thompson Scheme!!!
                                ;xnor and xnog are in the unit of /(m^-4).

      (schemename[ischeme] eq 'thom'): dbz_vt_calc_thom_BB, qr,qs,qg,pp,temp,qv,nx,ny,nz,nt,dbz,rcp_lamr,rcp_lamg,dens, $
                                       qnr,nr,smo4,res,rhos,xDs,xDs_m,xnor,xnog,lamr,lamg,dbz_r,dbz_s,dbz_g,r_vt_new
    endcase

    qnrv=nr          ;name the number concentration per unit volumn (unit is /(m^3)) to the same name as output for Morrison 
                 ;Looks a bit crazy here in the naming methods for number number concentrations!!
    print,'max dbz=',max(dbz)
    print,'max dbz_r=',max(dbz_r),'  max dbz_s=',max(dbz_s),'  max dbz_g=',max(dbz_g)
    print,'max rhos=',max(rhos),'  min rhos=',min(rhos)

    print,'max r_vt_new=',max(r_vt_new), ' min r_vt_new=',min(r_vt_new)

    ;Give the x and y index (from the 4th domain with 1.3 km resolution) for the profiler sites:
    ix_ATA=191
    iy_ATA=293
    ix_CZC=34
    iy_CZC=249
    ix_CFC=179
    iy_CFC=283

    if (i eq 0)  then begin
      ;y-coordinate
      ;I assume the vertical levels (H_MSL) are from 0.0 to 20 km MSL with 0.25 km interval
      ;The Sprof_lev_ATA is above the local terrain height ter[ix,iy]
      H_MSL=findgen(81)/4. ;0-20 km above MSL with 0.25 km intervals
      ter=ter/1000.        ;convert m to km
      Sprof_lev_ATA=H_MSL(where(H_MSL ge ter[ix_ATA,iy_ATA]))
      Sprof_lev_CZC=H_MSL(where(H_MSL ge ter[ix_CZC,iy_CZC]))
      Sprof_lev_CFC=H_MSL(where(H_MSL ge ter[ix_CFC,iy_CFC]))
      dim_ATA=size(Sprof_lev_ATA)
      num_lev_ATA=dim_ATA[1]
      dim_CZC=size(Sprof_lev_CZC)
      num_lev_CZC=dim_CZC[1]
      dim_CFC=size(Sprof_lev_CFC)
      num_lev_CFC=dim_CFC[1]
      print,num_lev_ATA,num_lev_CZC,num_lev_CFC
      ;print,Sprof_lev_ATA,Sprof_lev_CZC,Sprof_lev_CFC

      ;Initialize variables
      temp_ATA=fltarr(num_time,num_lev_ATA)
      dens_ATA=fltarr(num_time,num_lev_ATA)
      pp_ATA=fltarr(num_time,num_lev_ATA)
      ghtp_ATA=fltarr(num_time,num_lev_ATA)
      uu_ATA=fltarr(num_time,num_lev_ATA)
      vv_ATA=fltarr(num_time,num_lev_ATA)
      ww_ATA=fltarr(num_time,num_lev_ATA)
      qc_ATA=fltarr(num_time,num_lev_ATA)
      qi_ATA=fltarr(num_time,num_lev_ATA)
      qv_ATA=fltarr(num_time,num_lev_ATA)
      qr_ATA=fltarr(num_time,num_lev_ATA)
      qs_ATA=fltarr(num_time,num_lev_ATA)
      qg_ATA=fltarr(num_time,num_lev_ATA)
      rcp_lamr_ATA=fltarr(num_time,num_lev_ATA)
      rcp_lams_ATA=fltarr(num_time,num_lev_ATA)
      rcp_lamg_ATA=fltarr(num_time,num_lev_ATA)
      xnor_ATA=fltarr(num_time,num_lev_ATA)
      xnog_ATA=fltarr(num_time,num_lev_ATA)
      dbz_ATA=fltarr(num_time,num_lev_ATA)
      smo4_ATA=fltarr(num_time,num_lev_ATA)
      rvt_new_ATA=fltarr(num_time,num_lev_ATA)
      qnrv_ATA=fltarr(num_time,num_lev_ATA)
      res_ATA=fltarr(num_time,num_lev_ATA)
      lamr_ATA=fltarr(num_time,num_lev_ATA)
      lamg_ATA=fltarr(num_time,num_lev_ATA)
      xDs_ATA=fltarr(num_time,num_lev_ATA)
      xDs_m_ATA=fltarr(num_time,num_lev_ATA)
      rhos_ATA=fltarr(num_time,num_lev_ATA)

      ......
    endif
  

    ;Interpolate 3D array to the time height (TH) 2D array at a specific site (grid point)
    ;ATA:
    temp_ATA[i,*]=interpol(temp[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    dens_ATA[i,*]=interpol(dens[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    pp_ATA[i,*]=interpol(pp[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    ghtp_ATA[i,*]=interpol(ghtp[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    uu_ATA[i,*]=interpol(uu[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    vv_ATA[i,*]=interpol(vv[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    ww_ATA[i,*]=interpol(ww[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    qc_ATA[i,*]=interpol(qc[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    qi_ATA[i,*]=interpol(qi[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    qv_ATA[i,*]=interpol(qv[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    qr_ATA[i,*]=interpol(qr[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    qs_ATA[i,*]=interpol(qs[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    qg_ATA[i,*]=interpol(qg[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    if (schemename[ischeme] eq 'thom') then begin
      xnor_ATA[i,*]=interpol(xnor[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
      xnog_ATA[i,*]=interpol(xnog[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
      smo4_ATA[i,*]=interpol(smo4[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    endif
    dbz_ATA[i,*]=interpol(dbz[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    rvt_new_ATA[i,*]=interpol(r_vt_new[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    qnrv_ATA[i,*]=interpol(qnrv[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    res_ATA[i,*]=interpol(res[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    lamr_ATA[i,*]=interpol(lamr[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    lamg_ATA[i,*]=interpol(lamg[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    xDs_ATA[i,*]=interpol(xDs[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    xDs_m_ATA[i,*]=interpol(xDs_m[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)
    rhos_ATA[i,*]=interpol(rhos[ix_ATA,iy_ATA,*],ghtp[ix_ATA,iy_ATA,*],Sprof_lev_ATA)

    ......
  endfor ; endfor of num_time

close, 10

;SAVE time-height data on each site to NetCDF file:
;For Thompson:
  save_netcdf,schemename(ischeme),'ATA',num_time,num_lev_ATA,$
              Sprof_lev_ATA,temp_ATA,dens_ATA,pp_ATA,ghtp_ATA,uu_ATA,vv_ATA,ww_ATA, $
              qc_ATA,qi_ATA,qv_ATA,qr_ATA,qs_ATA,qg_ATA, $
              xnor_ATA,xnog_ATA,dbz_ATA,smo4_ATA,qnrv_ATA,res_ATA,lamr_ATA,lamg_ATA, $
              xDs_ATA,xDs_m_ATA,rhos_ATA,rvt_new_ATA

  ......
endfor  ; endfor of scheme

end  ; end of the pro comp_th_data2prof_netcdf_Thom_BB
  
