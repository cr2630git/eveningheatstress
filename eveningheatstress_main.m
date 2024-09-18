icloud='~/Library/Mobile Documents/com~apple~CloudDocs/';
figloc=strcat(icloud,'General_Academics/Research/Evening_Heat_Stress/');
era5loc='/Volumes/ExternalDriveD/ERA5_WRFauxiliary/';
era5landloc='/Volumes/ExternalDriveD/ERA5Land_Mideast/';
pwsdataloc=strcat(icloud,'General_Academics/Research/LCLUC_Tuholske_2023-26/Data_Sources/wunderground_stns/');

saveloc='/Volumes/ExternalDriveD/Various_LCLUC_arrays/';


dostartuptasks=0; %10 min; runs all "required on start-up" tasks

reloaddata=0; %8 min; required on start-up
readstndataanadplotmultistntimeseries=0; %1 min; need to check settings when running for the first time (see loop itself for details)
    makeplot=0;
readpwsdata=0; %3 min
domultivartimeseries=0; %10 sec; requires having already run wrfrunanalysis
era5singlelevelreading_main=0; %30 min; once computed, heavy time-consuming elements can largely be reloaded
era5pressurelevelreading_main=0; %15 min
    alsodouandv=0; %sometimes only need to do T and q
era5landreading_main=0; %30 min
era5_twcalc_part1=0; %10 min for single level/all hours, same for pressure levels/4-hourly
    if era5_twcalc_part1==1;datasettodo='era5_pl';end %'era5_sl', 'era5_pl', or 'era5land'
era5_twcalc_part2=0; %4 min
readmoreradiosondelevels=0; %4 hr
reloadrharmradiosondes=0; %1.5 min
diurnalcompositesetup=0; %90 min to redo everything
era5p95diurnalcomposites=0; %4 min; prep for fig 3
makefig3_revised=0; %4 min when reloading arrays, 30 sec otherwise; creates Fig 3
makefig2_asboxplots=1; %10 sec; creates Fig 2
calcspreads_era5=0; %6 min
    minfrac=0.5; %default (for main Fig 1) is 0.5; for supplemental figure, 0.4; ALSO APPLIES TO CALCSPREADS_POINTDATA
getvariousstats=0; %2 min; required on start-up
calcspreads_pointdata=0; %3 sec; required on start-up
makefig1_revised=0; %45 sec; creates Fig 1
    if makefig1_revised==1;datasettodo='era5';end
calcs_diurnalsequence=0; %30 min to restart all (i.e. when changing hotdaysbasedon)
    hotdaysbasedon='abudhabi'; %'abudhabi' [default] or 'liwaoasis' (from besti_era5 array)
finaldiurnalseqprep=0; %5 sec; required on start-up
map_diurnalsequence=0; %25 sec; creates Fig 4
getvertprofiles=0;
    vertprofileprep=0; 
        dorharmprep=0; %2 min
makefig5_revised=0; %15 sec for plotting alone; creates Fig 5
    doera5prep=0; %15 min; requires having done era5pressurelevelreading_main
readutciandwbgtdata=0; %1 hour; UTCI and WBGT were calculated in Jupyter notebook utci_wbgt_full.ipynb
readlurompsheatindex=0; %8 min
calcheatindex=0; %15 min
getaddlindexstats=0; %25 min; for HI, UTCI, WBGT, Lu+Romps HI
calcspreads_addlindices=0; %40 min
heatindexcomparison=0; %15 sec
makefig1_suppversion_addlindices=0; %1 min
suppfigwinddirs=0; %1 min 20 sec
    moresuppcalcs=1;
    plotwindshiftmap=0;
    plotabudhabivsdubaidetail=1;
    plotdiurnalqrange=0;
   

if dostartuptasks==1;setup_nctoolbox;reloaddata=1;getvariousstats=1;calcspreads_pointdata=1;finaldiurnalseqprep=1;end

f=load(strcat(saveloc,'lcluchotdaysoutput_',hotdaysbasedon));
tarray_twhotdays_mean=f.tarray_twhotdays_mean;qarray_twhotdays_mean=f.qarray_twhotdays_mean;
uarray_twhotdays_mean=f.uarray_twhotdays_mean;varray_twhotdays_mean=f.varray_twhotdays_mean;
tarray_thotdays_mean=f.tarray_thotdays_mean;qarray_thotdays_mean=f.qarray_thotdays_mean;
uarray_thotdays_mean=f.uarray_thotdays_mean;varray_thotdays_mean=f.varray_thotdays_mean;



%Lat/lon bounds of maps that will be produced
vsmmap_w=50.24;vsmmap_n=27.52;vsmmap_e=57.5;vsmmap_s=22;
smmap_w=50.24;smmap_n=29;smmap_e=60;smmap_s=20;
lgmap_w=39.24;lgmap_n=34;lgmap_e=63;lgmap_s=15;

exist zoomedin;if ans==0;zoomedin=0;end %default
if zoomedin==1
    wb=vsmmap_w;nb=vsmmap_n;eb=vsmmap_e;sb=vsmmap_s;
else
    wb=smmap_w;nb=smmap_n;eb=smmap_e;sb=smmap_s;
end

%Settings for ERA5 data
monlens=[31;28;31;30;31;30;31;31;30;31;30;31];
climodaylen=92; %days in JJA
startyr=2000;stopyr=2020;nyr=stopyr-startyr+1;

%Any specific heat wave to focus on?
hw_year=2020;hw_firstdoy=229;hw_lastdoy=235;numhwdays=hw_lastdoy-hw_firstdoy+1;

%PG/Arabia subregions
regnames={'sPG';'wPG';'inland'};
regnames_full={'Southern Persian Gulf';'Western Persian Gulf';'Inland Saudi Arabia'};


fname=strcat(era5loc,'era5_singlelevel_abudhabi_',num2str(startyr),'06.nc');
lat1d=ncread(fname,'latitude');lon1d=ncread(fname,'longitude');
[lat2d_era5,lon2d_era5]=latlon_2dfrom1d(lat1d,lon1d);
lat_era5=double(lat2d_era5');lon_era5=double(lon2d_era5');
latnative_era5=double(lat2d_era5);lonnative_era5=double(lon2d_era5);
era5latsz=length(lat1d);era5lonsz=length(lon1d);

%Higher-res lat/lon arrays
[Xorig,Yorig]=meshgrid(1:1:161);
[Xq,Yq]=meshgrid(1:0.1:161);
lat_era5_10xres=interp2(Xorig,Yorig,lat_era5,Xq,Yq);
lon_era5_10xres=interp2(Xorig,Yorig,lon_era5,Xq,Yq);

splabels={'a)';'b)';'c)';'d)';'e)';'f)';'g)';'h)';'i)'};


%Also get closest ERA5 pt for several points of interest
%1. Abu Dhabi Airport (24.44N, 54.65E)
%2. Persian Gulf waters 50 km northwest of Abu Dhabi (24.80, 54.05)
%3. rural coastline west of city (23.98, 53.95)
%4. Bateen Liwa inland oasis (23.13, 53.81)
%5. Dubai airport (25.24, 55.39)
latsofint=[24.44;24.80;23.98;23.13;25.24];lonsofint=[54.65;54.05;53.95;53.81;55.39];
clear besti_era5;clear bestj_era5;
for loop=1:size(latsofint,1)
    bestlatdiffsofar=1000;bestlondiffsofar=1000;
    for i=1:size(lat_era5,1)
        for j=1:size(lat_era5,2)
            thislat=lat_era5(i,j);thislon=lon_era5(i,j);
            latdiff=abs(thislat-latsofint(loop));londiff=abs(thislon-lonsofint(loop));
            if latdiff<=bestlatdiffsofar
                besti_era5(loop)=i;
                bestlatdiffsofar=latdiff;
            end
            if londiff<=bestlondiffsofar
                bestj_era5(loop)=j;
                bestlondiffsofar=londiff;
            end
        end
    end
end
loccolors=[colors('medium green');colors('somewhat dark blue');colors('sirkes_orange');colors('somewhat dark red')];
newregcolors=[colors('somewhat dark blue');colors('medium green');colors('somewhat dark red')];

%Where to define based on 
if strcmp(hotdaysbasedon,'abudhabi')
    myi=besti_era5(1);myj=bestj_era5(1);mylat=latsofint(1);mylon=lonsofint(1);
elseif strcmp(hotdaysbasedon,'liwaoasis')
    myi=besti_era5(4);myj=bestj_era5(4);mylat=latsofint(4);mylon=lonsofint(4);
else
    disp('Need help in determining what to do!');return;
end

abudhabiurbanlatbounds=[24.46;24.54;24.46;24.50;24.44;24.24;24.31;24.38;24.46];
abudhabiurbanlonbounds=[54.30;54.38;54.45;54.61;54.77;54.62;54.44;54.46;54.30];

dubaisharjahurbanlatbounds=[24.98;25.15;25.29;25.49;25.31;25.15;25.13;24.96;24.93;24.98];
dubaisharjahurbanlonbounds=[55.01;55.13;55.27;55.49;55.55;55.42;55.28;55.13;55.06;55.01];

%20-km-wide strip along coast
sPGcoastallatbounds=[24.45;24.60;25.90;25.65;25.45;25.20;24.90;24.65;24.50;24.28;24.09;23.87;23.62;24.45];
sPGcoastallonbounds=[51.60;54.00;55.75;55.75;55.65;55.55;55.25;54.90;54.83;54.72;54.50;54.03;52.03;51.60];

tzadj=3; %from utc to lst, add this -- 4 time zones away, but data starts at an index of 1 whereas hours start at 00


if reloaddata==1
    %Reload hourly data from ERA5, as UTC
    %Reminder: contains hourly data for 92 days per year (JJA) for 2000-2020
    %Load very large arrays only as needed
    tmp=load(strcat(saveloc,'twera5.mat'));tw_era5=tmp.tw_era5;clear tmp;
    %tmp=load(strcat(saveloc,'tera5.mat'));t_era5=tmp.t_era5;clear tmp;
    %tmp=load(strcat(saveloc,'qera5.mat'));q_era5=tmp.q_era5;clear tmp;
    %tmp=load(strcat(saveloc,'uera5.mat'));u_era5=tmp.u_era5;clear tmp;
    %tmp=load(strcat(saveloc,'vera5.mat'));v_era5=tmp.v_era5;clear tmp;
    %tmp=load(strcat(saveloc,'psfcera5.mat'));psfc_era5=tmp.psfc_era5;clear tmp; %only needed for calculating Tw 
    %tmp=load(strcat(saveloc,'pblhera5.mat'));pblh_era5=tmp.pblh_era5;clear tmp;
    %tmp=load(strcat(saveloc,'hiera5.mat'));hi_era5=tmp.hi_era5;clear tmp;
    %tmp=load(strcat(saveloc,'utciera5.mat'));utci_era5=tmp.utci_era5;clear tmp;

    exist nlat;if ans==0;nlat=size(tw_era5,1);nlon=size(tw_era5,2);nhrs=size(tw_era5,3);ndays=nhrs/24;end
    tw2m_days=reshape(tw_era5,[nlat nlon 24 climodaylen*nyr]);

    tmp=load(strcat(saveloc,'tera5.mat'));t_era5=tmp.t_era5;clear tmp;
        t2m_days=reshape(t_era5,[nlat nlon 24 climodaylen*nyr]);clear t_era5;

    tmp=load(strcat(saveloc,'qera5.mat'));q_era5=tmp.q_era5;clear tmp;
        q2m_days=reshape(q_era5,[nlat nlon 24 climodaylen*nyr]);clear q_era5;

    tmp=load(strcat(saveloc,'hotdaysstarts.mat'));
    tw95_gridcell=tmp.tw95_gridcell;t95_gridcell=tmp.t95_gridcell;
    q90_gridcell=tmp.q90_gridcell;q75_gridcell=tmp.q75_gridcell;
    t90_gridcell=tmp.t90_gridcell;t75_gridcell=tmp.t75_gridcell;
    %tw95whereabove=tmp.tw95whereabove;t95whereabove=tmp.t95whereabove;


    %Reload ERA5-Land data only as needed, because it's large and Matlab likes to crash
    %Also, reduce size of spatial domain to help with this issue
    %tmp=load(strcat(saveloc,'twera5land.mat'));tw_era5land=tmp.tw_era5land;clear tmp;
    %tmp=load(strcat(saveloc,'tera5land.mat'));t_era5land=tmp.t_era5land;t_era5land=t_era5land(41:161,61:261,:);clear tmp;
    %tmp=load(strcat(saveloc,'qera5land.mat'));q_era5land=tmp.q_era5land;q_era5land=q_era5land(41:161,61:261,:);clear tmp;
    %tmp=load(strcat(saveloc,'uera5land.mat'));u_era5land=tmp.u_era5land;clear tmp;
    %tmp=load(strcat(saveloc,'vera5land.mat'));v_era5land=tmp.v_era5land;clear tmp;
    e5llat1d=ncread(strcat(era5landloc,'t_2020.nc'),'latitude');
    e5llon1d=ncread(strcat(era5landloc,'t_2020.nc'),'longitude');
    [lat_era5land,lon_era5land]=latlon_2dfrom1d(e5llat1d,e5llon1d);lat_era5land=lat_era5land';lon_era5land=lon_era5land';
    %Because of domain-size trimming, reduce lat and lon arrays here
    lat_era5land=double(lat_era5land(41:161,61:261));lon_era5land=double(lon_era5land(41:161,61:261));e5llatsz=size(lat_era5land,1);e5llonsz=size(lat_era5land,2);

    %Reload JJA diurnal means from ERA5
    tmp=load(strcat(saveloc,'diurnalmeansera5.mat'));
    tw2m_diurnalmean=tmp.tw2m_diurnalmean;t2m_diurnalmean=tmp.t2m_diurnalmean;q2m_diurnalmean=tmp.q2m_diurnalmean;
    pblh_diurnalmean=tmp.pblh_diurnalmean;
    tw2m_diurnalmax=tmp.tw2m_diurnalmax;t2m_diurnalmax=tmp.t2m_diurnalmax;
    winddir10_diurnalmode=tmp.winddir10_diurnalmode;windspd10_diurnalmean=tmp.windspd10_diurnalmean;
    hi2m_diurnalmean=tmp.hi2m_diurnalmean;hi2m_diurnalmax=tmp.hi2m_diurnalmax;
    utci2m_diurnalmean=tmp.utci2m_diurnalmean;utci2m_diurnalmax=tmp.utci2m_diurnalmax;
    wbgt2m_diurnalmean=tmp.wbgt2m_diurnalmean;wbgt2m_diurnalmax=tmp.wbgt2m_diurnalmax;
    loadthese=0;
    if loadthese==1 %large, so don't load by default
    tmp=load(strcat(saveloc,'diurnalmeansp95_pblh_era5.mat'));pblhp95arr=tmp.pblhp95arr;
    tmp=load(strcat(saveloc,'diurnalmeansp95_tw_era5.mat'));tw2mp95arr=tmp.tw2mp95arr;
    tmp=load(strcat(saveloc,'diurnalmeansp95_t_era5.mat'));t2mp95arr=tmp.t2mp95arr;
    tmp=load(strcat(saveloc,'diurnalmeansp95_q_era5.mat'));q2mp95arr=tmp.q2mp95arr;
    tmp=load(strcat(saveloc,'diurnalmeansp95_winddir_era5.mat'));winddir10p95arr=tmp.winddir10p95arr;
    tmp=load(strcat(saveloc,'diurnalmeansp95_windspd_era5.mat'));windspd10p95arr=tmp.windspd10p95arr;
    tmp=load(strcat(saveloc,'diurnalmeansp95_hi_era5.mat'));hi2mp95arr=tmp.hi2mp95arr;
    tmp=load(strcat(saveloc,'diurnalmeansp95_utci_era5.mat'));utci2mp95arr=tmp.utci2mp95arr;
    tmp=load(strcat(saveloc,'diurnalmeansp95_wbgt_era5.mat'));wbgt2mp95arr=tmp.wbgt2mp95arr;
    end

    %Reload central variables for ERA5
    tmp=load(strcat(saveloc,'centralvarsera5.mat'));
    centralhours_local_era5=tmp.centralhours_local_era5;
    centralwinddirs_era5=tmp.centralwinddirs_era5;
    centralpblhs_era5=tmp.centralpblhs_era5;
    centralhours_hi_local_era5=tmp.centralhours_hi_local_era5;
    centralhours_utci_local_era5=tmp.centralhours_utci_local_era5;
    centralhours_wbgt_local_era5=tmp.centralhours_wbgt_local_era5;

    %Various other ERA5 stats
    tmp=load(strcat(saveloc,'addlstats.mat'));
    hi_era5_alltimemax=tmp.hi_era5_alltimemax;hi_era5_p95=tmp.hi_era5_p95;
    utci_era5_alltimemax=tmp.utci_era5_alltimemax;utci_era5_p95=tmp.utci_era5_p95;
    wbgt_era5_alltimemax=tmp.wbgt_era5_alltimemax;wbgt_era5_p95=tmp.wbgt_era5_p95;

    %ERA5 vertical profiles
    tmp=load(strcat(saveloc,'hotdaysprofiles_new.mat'));
    twprofile_twhotdays_era5_mean=tmp.twprofile_twhotdays_era5_mean;twprofile_twhotdays_era5_stdev=tmp.twprofile_twhotdays_era5_stdev;
    tprofile_twhotdays_era5_mean=tmp.tprofile_twhotdays_era5_mean;tprofile_twhotdays_era5_stdev=tmp.tprofile_twhotdays_era5_stdev;
    qprofile_twhotdays_era5_mean=tmp.qprofile_twhotdays_era5_mean;qprofile_twhotdays_era5_stdev=tmp.qprofile_twhotdays_era5_stdev;
    twprofile_alldays_era5_mean=tmp.twprofile_alldays_era5_mean;twprofile_alldays_era5_stdev=tmp.twprofile_alldays_era5_stdev;
    tprofile_alldays_era5_mean=tmp.tprofile_alldays_era5_mean;tprofile_alldays_era5_stdev=tmp.tprofile_alldays_era5_stdev;
    qprofile_alldays_era5_mean=tmp.qprofile_alldays_era5_mean;qprofile_alldays_era5_stdev=tmp.qprofile_alldays_era5_stdev;

    tmp=load(strcat(saveloc,'lsmaskera5.mat'));lsmask_era5=tmp.lsmask_era5;
    lsmask_era5_10xres=interp2(Xorig,Yorig,lsmask_era5,Xq,Yq);

    %RHARM vertical profiles
    tmp=load(strcat(saveloc,'hotdaysprofiles_rharm.mat'));
    twprofile_alldays_rharm_byloc=tmp.twprofile_alldays_rharm_byloc;twprofile_twhotdays_rharm_byloc=tmp.twprofile_twhotdays_rharm_byloc;

    %Reload stn data
    %Reminder: contains hourly data for 365 days per year for 2000-2020 (time zone: UTC)
    f=load(strcat(saveloc,'savedstndata_inland.mat'));
    subdailyt_inland_tmp=f.subdailyt_all;subdailytd_inland_tmp=f.subdailytd_all;subdailytw_inland_tmp=f.subdailytw_all;
    subdailytimes_inland_tmp=f.subdailytimes_all;subdailywinddir_inland_tmp=f.subdailywinddir_all;subdailywindspd_inland_tmp=f.subdailywindspd_all;
    stninfo_inland=f.stninfo_all;
    f=load(strcat(saveloc,'savedstndata_sPG.mat'));
    subdailyt_sPG_tmp=f.subdailyt_all;subdailytd_sPG_tmp=f.subdailytd_all;subdailytw_sPG_tmp=f.subdailytw_all;
    subdailytimes_sPG_tmp=f.subdailytimes_all;subdailywinddir_sPG_tmp=f.subdailywinddir_all;subdailywindspd_sPG_tmp=f.subdailywindspd_all;
    stninfo_sPG=f.stninfo_all;
    f=load(strcat(saveloc,'savedstndata_wPG.mat'));
    subdailyt_wPG_tmp=f.subdailyt_all;subdailytd_wPG_tmp=f.subdailytd_all;subdailytw_wPG_tmp=f.subdailytw_all;
    subdailytimes_wPG_tmp=f.subdailytimes_all;subdailywinddir_wPG_tmp=f.subdailywinddir_all;subdailywindspd_wPG_tmp=f.subdailywindspd_all;
    stninfo_wPG=f.stninfo_all;

    %Condense so stn data timespan corresponds directly to that of ERA5 
    jjadays=repmat([zeros(151,1);[152:243]';zeros(122,1)],[nyr 1]);jjadays(jjadays>0)=1;
    jjadays=permute(repmat(jjadays,[1 size(subdailyt_inland_tmp,1) size(subdailyt_inland_tmp,3)]),[2 1 3]); %should now match dims of subdailyt_inland, etc
    test=jjadays==1;
    clear subdailytw_inland;clear subdailytw_sPG;clear subdailytw_wPG;
    clear subdailyt_inland;clear subdailyt_sPG;clear subdailyt_wPG;
    clear subdailytd_inland;clear subdailytd_sPG;clear subdailytd_wPG;
    clear subdailywinddir_inland;clear subdailywinddir_sPG;clear subdailywinddir_wPG;
    clear subdailywindspd_inland;clear subdailywindspd_sPG;clear subdailywindspd_wPG;
    %There's surely a more efficient way, but, after some wasted effort, this is what I understand and it works
    for d1=1:max([size(subdailyt_inland_tmp,1);size(subdailyt_sPG_tmp,1);size(subdailyt_wPG_tmp,1)])
        newd2=0;
        for d2=1:size(test,2)
            if test(d1,d2,1)~=0
                newd2=newd2+1;
            end
            for d3=1:size(test,3)
                if test(d1,d2,d3)~=0
                    if d1<=size(subdailyt_inland_tmp,1)
                        subdailytw_inland(d1,newd2,d3)=subdailytw_inland_tmp(d1,d2,d3);
                        subdailyt_inland(d1,newd2,d3)=subdailyt_inland_tmp(d1,d2,d3);
                        subdailytd_inland(d1,newd2,d3)=subdailytd_inland_tmp(d1,d2,d3);
                        subdailywinddir_inland(d1,newd2,d3)=subdailywinddir_inland_tmp(d1,d2,d3);
                        subdailywindspd_inland(d1,newd2,d3)=subdailywindspd_inland_tmp(d1,d2,d3);
                    end
                    if d1<=size(subdailyt_sPG_tmp,1)
                        subdailytw_sPG(d1,newd2,d3)=subdailytw_sPG_tmp(d1,d2,d3);
                        subdailyt_sPG(d1,newd2,d3)=subdailyt_sPG_tmp(d1,d2,d3);
                        subdailytd_sPG(d1,newd2,d3)=subdailytd_sPG_tmp(d1,d2,d3);
                        subdailywinddir_sPG(d1,newd2,d3)=subdailywinddir_sPG_tmp(d1,d2,d3);
                        subdailywindspd_sPG(d1,newd2,d3)=subdailywindspd_sPG_tmp(d1,d2,d3);
                    end
                    if d1<=size(subdailyt_wPG_tmp,1)
                        subdailytw_wPG(d1,newd2,d3)=subdailytw_wPG_tmp(d1,d2,d3);
                        subdailyt_wPG(d1,newd2,d3)=subdailyt_wPG_tmp(d1,d2,d3);
                        subdailytd_wPG(d1,newd2,d3)=subdailytd_wPG_tmp(d1,d2,d3);
                        subdailywinddir_wPG(d1,newd2,d3)=subdailywinddir_wPG_tmp(d1,d2,d3);
                        subdailywindspd_wPG(d1,newd2,d3)=subdailywindspd_wPG_tmp(d1,d2,d3);
                    end
                end
            end
        end
    end
    invalid=subdailywinddir_inland<0;subdailywinddir_inland(invalid)=NaN;
    invalid=subdailywinddir_sPG<0;subdailywinddir_sPG(invalid)=NaN;
    invalid=subdailywinddir_wPG<0;subdailywinddir_wPG(invalid)=NaN;
    invalid=subdailywindspd_inland<0;subdailywindspd_inland(invalid)=NaN;
    invalid=subdailywindspd_sPG<0;subdailywindspd_sPG(invalid)=NaN;
    invalid=subdailywindspd_wPG<0;subdailywindspd_wPG(invalid)=NaN;

    %Also reload PWS [personal weather station] data
    f=load(strcat(saveloc,'savedpwsdata.mat'));
    subdailytw_pws_hrs=f.subdailytw_pws_hrs;pwsstninfo=f.pwsstninfo;
    subdailyt_pws_hrs=f.subdailyt_pws_hrs;subdailytd_pws_hrs=f.subdailytd_pws_hrs;
    subdailywinddir_pws_hrs=f.subdailywinddir_pws_hrs;subdailywindspd_pws_hrs=f.subdailywindspd_pws_hrs;
    subdailypsfc_pws_hrs=f.subdailypsfc_pws_hrs;
    %Reminder: downloaded time zone is LST -- convert to UTC (04 LST = 00 UTC)
    subdailytw_pws_hrs_old=subdailytw_pws_hrs;
    subdailytw_pws_hrs(:,:,1:20)=subdailytw_pws_hrs_old(:,:,5:24);subdailytw_pws_hrs(:,:,21:24)=subdailytw_pws_hrs_old(:,:,1:4);
    subdailyt_pws_hrs_old=subdailyt_pws_hrs;
    subdailyt_pws_hrs(:,:,1:20)=subdailyt_pws_hrs_old(:,:,5:24);subdailyt_pws_hrs(:,:,21:24)=subdailyt_pws_hrs_old(:,:,1:4);
    subdailytd_pws_hrs_old=subdailytd_pws_hrs;
    subdailytd_pws_hrs(:,:,1:20)=subdailytd_pws_hrs_old(:,:,5:24);subdailytd_pws_hrs(:,:,21:24)=subdailytd_pws_hrs_old(:,:,1:4);

    f=readtable(strcat(pwsdataloc,'station_info.csv'));
    pwsstnnames=f.Var1;pwsstnlats=f.Var2;pwsstnlons=f.Var3;
    numpws=length(pwsstnnames);

    %Now, separate ERA5 gridpoints and stations (both official and PWSs) into 3 regions: 
        %(1) coastal waters (<50 km to mainland UAE)
        %(2) coastal land (<50 km to water)
        %(3) inland (UAE, w of 55.8E, >50 km to water)

    uaeshapefile=shaperead(strcat(icloud,'General_Academics/Research/KeyFiles/CountryShapefiles/UnitedArabEmirates_adm0.shp'));
        uaelats=uaeshapefile.Y;uaelons=uaeshapefile.X;
        %SW corner is in wrong location selon latest maps -- adjust
        for i=1:length(uaelats);if uaelats(i)-23<0.1 && uaelons(i)-52<0.1;uaelats(i)=22.94;uaelons(i)=52.58;end;end

    %ERA5 gridpt regionalization (1 min)
    %Also include a higher-res version for pleasingness of map (takes about 1h30 to recalculate)
    myregs_era5=NaN.*ones(size(lat_era5));
    for i=1:size(lat_era5,1)
        for j=1:size(lon_era5,2)
            mylat=lat_era5(i,j);mylon=lon_era5(i,j);lsmaskval=lsmask_era5(i,j);
            myregs_era5(i,j)=separateintouaeregions(mylat,mylon,lat_era5,lon_era5,lsmask_era5,lsmaskval,uaelats,uaelons);
        end
    end
    recalcthis=0;
    if recalcthis==1
        myregs_era5_10xres=NaN.*ones(size(lat_era5_10xres));
        for i=1:size(lat_era5_10xres,1)
            for j=1:size(lon_era5_10xres,2)
                mylat=lat_era5_10xres(i,j);mylon=lon_era5_10xres(i,j);lsmaskval=lsmask_era5_10xres(i,j);
                myregs_era5_10xres(i,j)=separateintouaeregions(mylat,mylon,lat_era5_10xres,lon_era5_10xres,lsmask_era5_10xres,lsmaskval,uaelats,uaelons);
            end
            if rem(i,20)==0;disp(i);disp(clock);end
        end
        myregs_era5_10xres_old=myregs_era5_10xres;
        %Edit to call 'coast' (region 2) points that are marked as land but not included in UAE shapefile polygon
        test=myregs_era5_10xres;
        for i=1:size(lat_era5_10xres,1)
            for j=1:size(lon_era5_10xres,2)
                if isnan(test(i,j)) && sum(~isnan(myregs_era5_10xres(1:i-1,j)))>=2 && sum(~isnan(myregs_era5_10xres(i+1:end,j)))>=2 && ...
                        max(myregs_era5_10xres(i-10:i+10,j))<=2 && min(myregs_era5_10xres(i:end,j))<=2 && ~(i<=805 && j<=775)
                    test(i,j)=2;
                end
            end
        end
        myregs_era5_10xres=test;clear test;
        save(strcat(saveloc,'regsera5'),'myregs_era5','myregs_era5_10xres','myregs_era5_10xres_old');
    else
        f=load(strcat(saveloc,'regsera5'));myregs_era5=f.myregs_era5;
        myregs_era5_10xres=f.myregs_era5_10xres;myregs_era5_10xres_old=f.myregs_era5_10xres_old;
    end

    %Get percentiles for each region, for both gridpts and stations
    approach=2; %1 -- everything together, then take percentiles; 2 -- get percentiles for each gridpt/stn, then spatially average

    %First, for ERA5
    for vl=1:3 %vars tw, t, q
        if vl==1;days=tw2m_days;elseif vl==2;days=t2m_days;elseif vl==3;days=q2m_days;end
        gridptregc=zeros(3,1);allgridpt={};ijsaver=cell(3,1);
        for reg=1:3
            for i=1:size(lat_era5,1)
                for j=1:size(lon_era5,2)
                    if myregs_era5(i,j)==reg
                        gridptregc(reg)=gridptregc(reg)+1;
                        ijsaver{reg}(gridptregc(reg),1)=i;ijsaver{reg}(gridptregc(reg),2)=j;
                        allgridpt{reg}(gridptregc(reg),:,:)=squeeze(days(i,j,:,:))';
                    end
                end
            end
        end
        gridptp10byregandhr=NaN.*ones(3,24);gridptp50byregandhr=NaN.*ones(3,24);gridptp90byregandhr=NaN.*ones(3,24);
        p10byreg={};p50byreg={};p90byreg={};
        for reg=1:3
            numgridpts=gridptregc(reg);
            if numgridpts~=0
                if approach==1
                    tmp=reshape(allgridpt{reg},[numgridpts*92*21 24]); 
                    gridptp10byregandhr(reg,:)=quantile(tmp,0.1);gridptp50byregandhr(reg,:)=quantile(tmp,0.5);
                        gridptp90byregandhr(reg,:)=quantile(tmp,0.9);
                else
                    clear p10;clear p50;clear p90;
                    for gridpt=1:size(allgridpt{reg},1)
                        p10(gridpt,:)=quantile(squeeze(allgridpt{reg}(gridpt,:,:)),0.1);
                        p50(gridpt,:)=quantile(squeeze(allgridpt{reg}(gridpt,:,:)),0.5);
                        p90(gridpt,:)=quantile(squeeze(allgridpt{reg}(gridpt,:,:)),0.9);
                    end
                    gridptp10byregandhr(reg,:)=mean(p10,'omitnan');
                    gridptp50byregandhr(reg,:)=mean(p50,'omitnan');
                    gridptp90byregandhr(reg,:)=mean(p90,'omitnan');
                end
            end
            p10byreg{reg}=p10;p50byreg{reg}=p50;p90byreg{reg}=p90;
        end

        if vl==1
            gridpttwp10byregandhr=gridptp10byregandhr;gridpttwp50byregandhr=gridptp50byregandhr;gridpttwp90byregandhr=gridptp90byregandhr;
            p10twbygridpt=p10byreg;p50twbygridpt=p50byreg;p90twbygridpt=p90byreg;
            allgridpttw=allgridpt;
        elseif vl==2
            gridpttp10byregandhr=gridptp10byregandhr;gridpttp50byregandhr=gridptp50byregandhr;gridpttp90byregandhr=gridptp90byregandhr;
            p10tbygridpt=p10byreg;p50tbygridpt=p50byreg;p90tbygridpt=p90byreg;
            allgridptt=allgridpt;
        elseif vl==3
            gridptqp10byregandhr=gridptp10byregandhr;gridptqp50byregandhr=gridptp50byregandhr;gridptqp90byregandhr=gridptp90byregandhr;
            p10qbygridpt=p10byreg;p50qbygridpt=p50byreg;p90qbygridpt=p90byreg;
            allgridptq=allgridpt;
        end
    end

    %Station regionalization
    pwsstnlatlons=[pwsstninfo{2} pwsstninfo{3} zeros(24,1)];
    stnstocheck=cat(1,stninfo_sPG,pwsstnlatlons); %all stns in UAE
    myregs_stns=[];regstnc=zeros(3,1);
    allstntw={};allstnt={};allstnq={};
    for stn=1:length(stnstocheck)
        mylat=stnstocheck(stn,1);mylon=stnstocheck(stn,2);
        reg=separateintouaeregions(mylat,mylon,lat_era5,lon_era5,lsmask_era5,1,uaelats,uaelons);

        %Find previously saved data for this station
        if stn<=4 %official -- in UTC 
            stninfo_here=stninfo_sPG;
            thistw=cat(1,squeeze(subdailytw_sPG(stn,:,:)),NaN.*ones(92*2,24)); %PWSs have 2 add'l years (2021-22), so fill NaNs to match
            thist=cat(1,squeeze(subdailyt_sPG(stn,:,:)),NaN.*ones(92*2,24));
            thisq=calcqfromTd(cat(1,squeeze(subdailytd_sPG(stn,:,:)),NaN.*ones(92*2,24)));
        else %PWS -- also in UTC
            thistw=squeeze(subdailytw_pws_hrs(stn-4,:,:));%thistw=cat(2,thistw_tmp(:,5:24),thistw_tmp(:,1:4));
            thist=squeeze(subdailyt_pws_hrs(stn-4,:,:));
            thisq=calcqfromTd(squeeze(subdailytd_pws_hrs(stn-4,:,:)));
        end

        if ~isnan(reg)
            regstnc(reg)=regstnc(reg)+1;
            myregs_stns{reg}(regstnc(reg),1)=mylat;myregs_stns{reg}(regstnc(reg),2)=mylon;
            
            allstntw{reg}(regstnc(reg),:,:)=thistw;
            allstnt{reg}(regstnc(reg),:,:)=thist;
            allstnq{reg}(regstnc(reg),:,:)=thisq;
        end
    end

    %Pool all stations in each region together and then get JJA percentiles
    for vl=1:3 %Tw, T, q
        if vl==1;allstn=allstntw;elseif vl==2;allstn=allstnt;elseif vl==3;allstn=allstnq;end

        stnp10byregandhr=NaN.*ones(3,24);stnp50byregandhr=NaN.*ones(3,24);stnp90byregandhr=NaN.*ones(3,24);clear samplesize;goodstns={};
        for reg=1:3 %water, coast, inland
            numstns=regstnc(reg);goodgridptc=0;
            if numstns~=0
                if approach==1
                    tmp=reshape(allstn{reg},[numstns*92*23 24]);
                    stnp10byregandhr(reg,:)=quantile(tmp,0.1);stnp50byregandhr(reg,:)=quantile(tmp,0.5);stnp90byregandhr(reg,:)=quantile(tmp,0.9);
                else
                    clear p10;clear p50;clear p90;
                    for gridpt=1:size(allstn{reg},1)
                        samplesize(gridpt,:)=squeeze(sum(~isnan(allstn{reg}(gridpt,:,:))));
    
                        if mean(samplesize(gridpt,2))>=100 %min sample size, in days
                            goodgridptc=goodgridptc+1;
                            %disp(gridpt); %troubleshooting
                            p10(goodgridptc,:)=quantile(squeeze(allstn{reg}(gridpt,:,:)),0.1);
                            p50(goodgridptc,:)=quantile(squeeze(allstn{reg}(gridpt,:,:)),0.5);
                            p90(goodgridptc,:)=quantile(squeeze(allstn{reg}(gridpt,:,:)),0.9);
                            goodstns{reg}(goodgridptc,1)=myregs_stns{reg}(gridpt,1);goodstns{reg}(goodgridptc,2)=myregs_stns{reg}(gridpt,2);
                        end
                    end
                    stnp10byregandhr(reg,:)=mean(p10,'omitnan');
                    stnp50byregandhr(reg,:)=mean(p50,'omitnan');
                    stnp90byregandhr(reg,:)=mean(p90,'omitnan');
                end
            end
        end

        if vl==1
            stntwp10byregandhr=stnp10byregandhr;stntwp50byregandhr=stnp50byregandhr;stntwp90byregandhr=stnp90byregandhr;
        elseif vl==2
            stntp10byregandhr=stnp10byregandhr;stntp50byregandhr=stnp50byregandhr;stntp90byregandhr=stnp90byregandhr;
        elseif vl==3
            stnqp10byregandhr=stnp10byregandhr;stnqp50byregandhr=stnp50byregandhr;stnqp90byregandhr=stnp90byregandhr;
        end
    end

    %Create version of ERA5 gridpt array that's weighted to match station
    %locations in each region
    exist t2m_days;
    if ans==0
        tmp=load(strcat(saveloc,'tera5.mat'));t_era5=tmp.t_era5;clear tmp;
        t2m_days=reshape(t_era5,[nlat nlon 24 climodaylen*nyr]);clear t_era5;
    end

    exist q2m_days;
    if ans==0
        tmp=load(strcat(saveloc,'qera5.mat'));q_era5=tmp.q_era5;clear tmp;
        q2m_days=reshape(q_era5,[nlat nlon 24 climodaylen*nyr]);clear q_era5;
    end

    for vl=1:3
        if vl==1;days=tw2m_days;elseif vl==2;days=t2m_days;elseif vl==3;days=q2m_days;end
        for reg=1:3 %water, coast, inland
            p10_here=NaN.*ones(length(goodstns{reg}),24);p10_here_alt=NaN.*ones(length(goodstns{reg}),24);
            p50_here=NaN.*ones(length(goodstns{reg}),24);p50_here_alt=NaN.*ones(length(goodstns{reg}),24);
            p90_here=NaN.*ones(length(goodstns{reg}),24);p90_here_alt=NaN.*ones(length(goodstns{reg}),24);
            for stn=1:length(goodstns{reg})
                stnlat=goodstns{reg}(stn,1);stnlon=goodstns{reg}(stn,2);
                %Get closest 4 ERA5 gridpts
                outarr=interpolate2dlatlonarray(stnlat,stnlon,lat_era5,lon_era5,'-180-180');
    
                for pct=1:3
                    if pct==1;pctile=0.1;elseif pct==2;pctile=0.5;elseif pct==3;pctile=0.9;end
                    p=outarr(1,4).*quantile(squeeze(days(outarr(1,1),outarr(1,2),:,:))',pctile)+...
                        outarr(2,4).*quantile(squeeze(days(outarr(2,1),outarr(2,2),:,:))',pctile)+...
                        outarr(3,4).*quantile(squeeze(days(outarr(3,1),outarr(3,2),:,:))',pctile)+...
                        outarr(4,4).*quantile(squeeze(days(outarr(4,1),outarr(4,2),:,:))',pctile);

                    %Also try using just single closest ERA5 gridpt
                    p_single=quantile(squeeze(days(outarr(1,1),outarr(1,2),:,:))',pctile);

                    if pct==1
                        p10_here(stn,:)=p;p10_here_alt(stn,:)=p_single;
                    elseif pct==2
                        p50_here(stn,:)=p;p50_here_alt(stn,:)=p_single;
                    elseif pct==3
                        p90_here(stn,:)=p;p90_here_alt(stn,:)=p_single;
                    end
                end
            end
            arr10_w(reg,:)=mean(p10_here,'omitnan');arr50_w(reg,:)=mean(p50_here,'omitnan');arr90_w(reg,:)=mean(p90_here,'omitnan');
            arr10_w_alt(reg,:)=mean(p10_here_alt,'omitnan');arr50_w_alt(reg,:)=mean(p50_here_alt,'omitnan');arr90_w_alt(reg,:)=mean(p90_here_alt,'omitnan');
        end
        if vl==1
            gridpttwp10byregandhr_w=arr10_w;gridpttwp50byregandhr_w=arr50_w;gridpttwp90byregandhr_w=arr90_w;
            gridpttwp10byregandhr_w_alt=arr10_w_alt;gridpttwp50byregandhr_w_alt=arr50_w_alt;gridpttwp90byregandhr_w_alt=arr90_w_alt;
        elseif vl==2
            gridpttp10byregandhr_w=arr10_w;gridpttp50byregandhr_w=arr50_w;gridpttp90byregandhr_w=arr90_w;
            gridpttp10byregandhr_w_alt=arr10_w_alt;gridpttp50byregandhr_w_alt=arr50_w_alt;gridpttp90byregandhr_w_alt=arr90_w_alt;
        elseif vl==3
            gridptqp10byregandhr_w=arr10_w;gridptqp50byregandhr_w=arr50_w;gridptqp90byregandhr_w=arr90_w;
            gridptqp10byregandhr_w_alt=arr10_w_alt;gridptqp50byregandhr_w_alt=arr50_w_alt;gridptqp90byregandhr_w_alt=arr90_w_alt;
        end
    end
    clear days;

    f=load(strcat(saveloc,'era5diurnalmeans'));
    pblhp95_diurnalmean=f.pblhp95_diurnalmean;t2mp95_diurnalmean=f.t2mp95_diurnalmean;tw2mp95_diurnalmean=f.tw2mp95_diurnalmean;
    q2mp95_diurnalmean=f.q2mp95_diurnalmean;winddir10p95_diurnalmode=f.winddir10p95_diurnalmode;windspd10p95_diurnalmean=f.windspd10p95_diurnalmean;
    hi2mp95_diurnalmean=f.hi2mp95_diurnalmean;utci2mp95_diurnalmean=f.utci2mp95_diurnalmean;wbgt2mp95_diurnalmean=f.wbgt2mp95_diurnalmean;
end


if readstndataanadplotmultistntimeseries==1
    %Tw timeseries for each station
    %Stations found/checked using completestndataanalysis.m script
    citynames_sPG={'Abu Dhabi AE';'Dubai AE';'Ras Al Khaimah AE';'Sharjah AE'};
    citynms_sPG={'abudhabi';'dubai';'rasalkhaimah';'sharjah'};
    citycodes_sPG=[1760;1757;1756;1758];
    citynames_wPG={'Doha QA';'Manama BH';'Kuwait KW'};
    citynms_wPG={'doha';'manama';'kuwait'};
    citycodes_wPG=[1755;1754;1705];
    citynames_GO={'Sohar OM';'Jiwani PK'}; %on second thought, omit these
    citynms_go={'sohar';'jiwani'};
    citycodes_GO=[1762;1783];
    citynames_inland={'Rafha SA';'Al Qaisumah SA';'Hail SA';'Buraydah SA';'Al Hofuf SA';'Riyadh SA';'Al Jouf SA'};
    citynms_inland={'rafha';'alqaisumah';'hail';'buraydah';'alhofuf';'riyadh';'aljouf'};
    citycodes_inland=[1692;1693;1696;1698;1700;1703;1691];

    if makeplot==1;figure(84);clf;end
    for reg=1:length(regnames)
        if makeplot==1;subplot(size(regnames,1),1,reg);end
        citynames=eval(['citynames_' regnames{reg}]);citycodes=eval(['citycodes_' regnames{reg}]);citynms=eval(['citynms_' regnames{reg}]);
        outputstnlist=citycodes;
        %mycolors=varycolor(size(citynames,1));
        if makeplot==1;mycb=colormaps('classy rainbow','more','pale');cbsz=size(mycb,1);end

        %Necessary settings for completestndataanalysis (also note that it requires data from ExternalDriveZ):
        %lookforstns=0;
        %maxmissing=25;
        %year1=startyr;year2=stopyr;doy1=1;doy2=365;
        %gett=1;gettd=1;computewetbulb=1;getwinddir=1;getwindspd=1;
        %daystogoout=0;tosave=0;
        %completestndataanalysis; %approx. 20 sec per stn obtained; output of interest is subdailytw_all whose first dimension should be stn number

        if makeplot==1
            for city=1:size(citynames,1)
                thiscityhwdata=squeeze(subdailytw_all(city,hw_firstdoy:hw_lastdoy,:))';
                cbidx=round(cbsz/(size(citynames,1)*2)+(city-1)*cbsz/size(citynames,1));
                plot(reshape(thiscityhwdata,[numhwdays*24 1]),'color',mycb(cbidx,:),'linewidth',1.5);
                %maketitle(citynames{city});
        
                ylim([10 35]);xlim([1 numhwdays*24]);
                set(gca,'xtick',1:24:numhwdays*24,'xticklabel',hw_firstdoy:hw_lastdoy);
                set(gca,'fontweight','bold','fontname','arial','fontsize',11);
                %Add dashed line corresponding to all-time max for Washington DC, to really put an exclamation mark on it
                hold on;
                if city==size(citynames,1)
                    if reg==size(regnames,1);xlabel('Day of Year in 2020','fontweight','bold','fontname','arial','fontsize',12);end
                    title(regnames_full{reg},'fontweight','bold','fontname','arial','fontsize',14);
    
                    myy=30.66.*ones(1,numhwdays*24);
                    myx=1:numhwdays*24;
                    plot(myx,myy,'b','linewidth',1.5,'linestyle','--');
                end
                %CHANGE TEXT TO BE DATA UNITS IN ORDER TO HAVE IT APPEAR USING EXPORT_FIG
                %t=text(1.01,1.1-0.15*city,citynames{city},'units','normalized');set(t,'color',mycb(cbidx,:),'fontweight','bold','fontname','arial','fontsize',11);
            end
        end

        %Manually recreate stninfo array because it wasn't saved before
        %if reg==1 %s PG
        %    stninfo_all=[24.43 54.65 27;25.26 55.36 10;25.61 55.94 31;25.33 55.52 34];
        %elseif reg==2 %w PG
        %    stninfo_all=[25.26 51.57 11;26.27 50.63 2;29.23 47.97 63];
        %elseif reg==3 %inland
        %    stninfo_all=[29.63 43.49 449;28.34 46.13 358;27.44 41.69 1015;26.3 43.77 648;25.29 49.49 179;24.71 46.73 635;29.79 40.1 689];
        %end

        %Save stn data for convenience
        save(strcat(saveloc,'savedstndata_',regnames{reg},'.mat'),'subdailytw_all','subdailytimes_all','subdailyt_all','subdailytd_all','subdailywinddir_all','subdailywindspd_all','stninfo_all');
        %...also save as csv for sharing with others
        for c=1:length(citynms)
            writematrix(squeeze(subdailytw_all(c,:,:)),strcat(saveloc,'stndata_',citynms{c},'_tw.txt'));
            writematrix(squeeze(subdailytimes_all(c,:,:)),strcat(saveloc,'stndata_',citynms{c},'_times.txt'));
            writematrix(squeeze(subdailyt_all(c,:,:)),strcat(saveloc,'stndata_',citynms{c},'_t.txt'));
            writematrix(squeeze(subdailytd_all(c,:,:)),strcat(saveloc,'stndata_',citynms{c},'_td.txt'));
            writematrix(squeeze(subdailywinddir_all(c,:,:)),strcat(saveloc,'stndata_',citynms{c},'_winddir.txt'));
            writematrix(squeeze(subdailywindspd_all(c,:,:)),strcat(saveloc,'stndata_',citynms{c},'_windspd.txt'));
            writematrix(squeeze(stninfo(c,:,:)),strcat(saveloc,'stndata_',citynms{c},'_metadata.txt'));
        end
    end
    %set(gcf,'color','w');figname=strcat(figloc,'twtimeseries');
    %curpart=1;highqualityfiguresetup;curpart=2;highqualityfiguresetup;
end



if readpwsdata==1
    nyr_pws=nyr+2; %so total length is 2000-2022; do this because many stations only have data for the last couple years

    subdailytimes_pws=NaN.*ones(nyr_pws*climodaylen,24*12);
    subdailyt_pws=NaN.*ones(length(pwsstnnames),nyr_pws*climodaylen,24*12); %dims are stns | days | hourly/subhourly [allows for 5-min intervals]
    subdailytd_pws=NaN.*ones(length(pwsstnnames),nyr_pws*climodaylen,24*12);
    subdailywinddir_pws=NaN.*ones(length(pwsstnnames),nyr_pws*climodaylen,24*12);
    subdailywindspd_pws=NaN.*ones(length(pwsstnnames),nyr_pws*climodaylen,24*12);
    subdailypsfc_pws=NaN.*ones(length(pwsstnnames),nyr_pws*climodaylen,24*12);
   
    yeararr_pws=[];for y=2000:2022;yeararr_pws=[yeararr_pws;repmat(y,[climodaylen 1])];end
    doyarr_pws=repmat([152:152+climodaylen-1]',[nyr+2 1]);
    hourarr_pws=[repmat(0,[12 1]);repmat(1,[12 1]);repmat(2,[12 1]);repmat(3,[12 1]);repmat(4,[12 1]);repmat(5,[12 1]);repmat(6,[12 1]);repmat(7,[12 1]);...
        repmat(8,[12 1]);repmat(9,[12 1]);repmat(10,[12 1]);repmat(11,[12 1]);repmat(12,[12 1]);repmat(13,[12 1]);repmat(14,[12 1]);repmat(15,[12 1]);...
        repmat(16,[12 1]);repmat(17,[12 1]);repmat(18,[12 1]);repmat(19,[12 1]);repmat(20,[12 1]);repmat(21,[12 1]);repmat(22,[12 1]);repmat(23,[12 1])];
    minutearr_pws=[repmat(linspace(0,55,12)',[24 1])];

    for s=1:length(pwsstnnames)
        if strcmp(pwsstnnames{s},"I90581516") || strcmp(pwsstnnames{s},"IABUDH16") || strcmp(pwsstnnames{s},"IABUDH36") || strcmp(pwsstnnames{s},"IABUDHAB100") ||...
            strcmp(pwsstnnames{s},"IABUDHAB116") || strcmp(pwsstnnames{s},"IABUDHAB141") || strcmp(pwsstnnames{s},"IALHASA3") ||...
            strcmp(pwsstnnames{s},"IDUBAI108")
            firstyear=2015;
        else
            firstyear=2020;
        end
        for yr=firstyear:2022
            mytable=readtable(strcat(pwsdataloc,pwsstnnames{s},'_',num2str(yr),'.csv'));
            if size(mytable,1)>=3 %some data exists
                for row=2:size(mytable,1)
                    tmp=char(mytable.Date(row));
                    if strcmp(tmp(5),'/') %e.g. of form 2020/01/01
                        thismon=str2num(tmp(6:7));thisdom=str2num(tmp(9:10));
                    else %e.g. of form 01/01/2020
                        thismon=str2num(tmp(1:2));thisdom=str2num(tmp(4:5));
                    end

                    if thismon>=6 && thismon<=8 %%JJA only
                        thisdoyrel=DatetoDOY(thismon,thisdom,2021)-151; %so Jun 1 is always in position 1
                        dayidx=(yr-2000)*climodaylen+thisdoyrel;

                        %Convert time to 24-hour time (00:00-23:59)
                        tmp=char(mytable.Time(row));
                        if strcmp(tmp(2),':') %e.g. of form 1:00 AM
                            curhour=str2num(tmp(1));thisminute=str2num(tmp(3:4));curampm=tmp(6:7);
                        else %e.g. of form 01:00 AM
                            curhour=str2num(tmp(1:2));thisminute=str2num(tmp(4:5));curampm=tmp(7:8);
                        end
                        if curhour==12 && strcmp(curampm,'AM')
                            thishour=0;
                        elseif curhour==12 && strcmp(curampm,'PM')
                            thishour=12;
                        elseif strcmp(curampm,'PM')
                            thishour=curhour+12;
                        else
                            thishour=curhour;
                        end
                        %Where in the day does this 5-minute interval fall?
                        subhourlyidx=thishour*12+round2(thisminute/5,1,'floor')+1;

                        %Register temperature, dewpoint, wind, pressure recordings at this time
                        if ~iscell(mytable.Temperature_C(row))
                            subdailyt_pws(s,dayidx,subhourlyidx)=mytable.Temperature_C(row);
                        end

                        if ~iscell(mytable.Dew_Point_C(row))
                            subdailytd_pws(s,dayidx,subhourlyidx)=mytable.Dew_Point_C(row);
                        end

                        if ~iscell(mytable.Speed_kmh(row))
                            subdailywindspd_pws(s,dayidx,subhourlyidx)=mytable.Speed_kmh(row);
                            %if subdailywindspd_pws(s,dayidx,subhourlyidx)~=0;disp('333');return;end
                        end

                        if ~(isnan(subdailywindspd_pws(s,dayidx,subhourlyidx)) || subdailywindspd_pws(s,dayidx,subhourlyidx)==0)
                            %if speed is NaN or 0, wind dir must be too
                            if ~strcmp(mytable.Wind{row},'NA')
                                subdailywinddir_pws(s,dayidx,subhourlyidx)=winddirfromtextdescrip(char(mytable.Wind(row)));
                            end
                        end

                        if ~iscell(mytable.Pressure_hPa(row))
                            subdailypsfc_pws(s,dayidx,subhourlyidx)=mytable.Pressure_hPa(row);
                        end
                    end
                end
            end
        end
        disp(s);disp(clock);
    end
    pwsstninfo={pwsstnnames pwsstnlats pwsstnlons};

    invalid=subdailyt_pws==0;subdailyt_pws(invalid)=NaN;
    invalid=subdailytd_pws==0;subdailytd_pws(invalid)=NaN;
    invalid=subdailypsfc_pws==0;subdailypsfc_pws(invalid)=NaN;

    subdailytw_pws=calcwbt_daviesjones(subdailyt_pws,subdailypsfc_pws.*100,calcqfromTd(subdailytd_pws)./1000);

    %Save for convenience
    save(strcat(saveloc,'savedpwsdata.mat'),'subdailytw_pws','subdailyt_pws','subdailytd_pws','subdailywinddir_pws','subdailywindspd_pws','subdailypsfc_pws','pwsstninfo','-append');


    %For comparison purposes, also create a version that consists of hourly means of this same data
    subdailytw_pws_hrs=NaN.*ones(size(subdailytw_pws,1),size(subdailytw_pws,2),size(subdailytw_pws,3)/12);
    subdailyt_pws_hrs=NaN.*ones(size(subdailytw_pws,1),size(subdailytw_pws,2),size(subdailytw_pws,3)/12);
    subdailytd_pws_hrs=NaN.*ones(size(subdailytw_pws,1),size(subdailytw_pws,2),size(subdailytw_pws,3)/12);
    subdailywinddir_pws_hrs=NaN.*ones(size(subdailytw_pws,1),size(subdailytw_pws,2),size(subdailytw_pws,3)/12);
    subdailywindspd_pws_hrs=NaN.*ones(size(subdailytw_pws,1),size(subdailytw_pws,2),size(subdailytw_pws,3)/12);
    subdailypsfc_pws_hrs=NaN.*ones(size(subdailytw_pws,1),size(subdailytw_pws,2),size(subdailytw_pws,3)/12);
    for hr=1:24
        int1=hr*12-11;int2=hr*12; %all 5-min intervals in this hour
        subdailytw_pws_hrs(:,:,hr)=mean(subdailytw_pws(:,:,int1:int2),3,'omitnan');
        subdailyt_pws_hrs(:,:,hr)=mean(subdailyt_pws(:,:,int1:int2),3,'omitnan');
        subdailytd_pws_hrs(:,:,hr)=mean(subdailytd_pws(:,:,int1:int2),3,'omitnan');
        subdailywinddir_pws_hrs(:,:,hr)=mean(subdailywinddir_pws(:,:,int1:int2),3,'omitnan');
        subdailywindspd_pws_hrs(:,:,hr)=mean(subdailywindspd_pws(:,:,int1:int2),3,'omitnan');
        subdailypsfc_pws_hrs(:,:,hr)=mean(subdailypsfc_pws(:,:,int1:int2),3,'omitnan');
    end
    invalid=subdailywinddir_pws_hrs<0;subdailywinddir_pws_hrs(invalid)=NaN;
    invalid=subdailywindspd_pws_hrs<0;subdailywindspd_pws_hrs(invalid)=NaN;

    save(strcat(saveloc,'savedpwsdata.mat'),'subdailytw_pws_hrs','subdailyt_pws_hrs','subdailytd_pws_hrs',...
        'subdailywinddir_pws_hrs','subdailywindspd_pws_hrs','subdailypsfc_pws_hrs','-append');
end

if era5singlelevelreading_main==1
    %Read all single-level data for JJA 2000-2020
    t2m=NaN.*ones(era5lonsz,era5latsz,nyr*climodaylen*24);q2m=NaN.*ones(era5lonsz,era5latsz,nyr*climodaylen*24);psfc=NaN.*ones(era5lonsz,era5latsz,nyr*climodaylen*24);
    u10=NaN.*ones(era5lonsz,era5latsz,nyr*climodaylen*24);v10=NaN.*ones(era5lonsz,era5latsz,nyr*climodaylen*24);pblh=NaN.*ones(era5lonsz,era5latsz,nyr*climodaylen*24);
    hrstart=1;

    for y=startyr:stopyr
        ystr=num2str(y);
        for mon=6:8
            if mon<=9;zm='0';else;zm='';end
            thislen=monlens(mon);
            fname=strcat(era5loc,'era5_singlelevel_abudhabi_',ystr,zm,num2str(mon),'.nc');
            psfc(:,:,hrstart:hrstart+thislen*24-1)=ncread(fname,'sp');
            t2m(:,:,hrstart:hrstart+thislen*24-1)=ncread(fname,'t2m')-273.15;
            q2m(:,:,hrstart:hrstart+thislen*24-1)=calcqfromTd_dynamicP(ncread(fname,'d2m')-273.15,ncread(fname,'sp'));
            u10(:,:,hrstart:hrstart+thislen*24-1)=ncread(fname,'u10');v10(:,:,hrstart:hrstart+thislen*24-1)=ncread(fname,'v10');
            pblh(:,:,hrstart:hrstart+thislen*24-1)=ncread(fname,'blh');
            hrstart=hrstart+thislen*24;
        end
    end

    %Rotate everything to normal orientation (1 min for 20 years)
    t_era5=permute(t2m,[2 1 3]);clear t2m;q_era5=permute(q2m,[2 1 3]);clear q2m;
    u_era5=permute(u10,[2 1 3]);clear u10;v_era5=permute(v10,[2 1 3]);clear v10;
    psfc_era5=permute(psfc,[2 1 3]);clear psfc;pblh_era5=permute(pblh,[2 1 3]);clear pblh;

    %Save again (20 min, but saves lots of time later in case of any crashes)
    u_era5=round2(u_era5,0.01);v_era5=round2(v_era5,0.01);pblh_era5=round2(pblh_era5,0.1); %extra precision not needed, and reduces array sizes pleasantly
    save(strcat(saveloc,'tera5.mat'),'t_era5','-v7.3');save(strcat(saveloc,'qera5.mat'),'q_era5','-v7.3');
    save(strcat(saveloc,'uera5.mat'),'u_era5','-v7.3');save(strcat(saveloc,'vera5.mat'),'v_era5','-v7.3');
    save(strcat(saveloc,'psfcera5.mat'),'psfc_era5','-v7.3');save(strcat(saveloc,'pblhera5.mat'),'pblh_era5','-v7.3');
end

if era5pressurelevelreading_main==1
    %Read all pressure-level data for JJA 2000-2020
    %6 levels obtained are 500, 700, 850, 925, 975, 1000 (last one should approximate the single-level i.e. near-surface data)
    %To avoid overloading Matlab, keep only data every 4 hours (i.e. 6 time points per day)
    %Default time points are (UTC): 00, 04, 08, 12, 16, 20
    tlevs=NaN.*ones(era5lonsz,era5latsz,6,nyr*climodaylen*6);qlevs=NaN.*ones(era5lonsz,era5latsz,6,nyr*climodaylen*6);
    if alsodouandv==1;ulevs=NaN.*ones(era5lonsz,era5latsz,6,nyr*climodaylen*6);vlevs=NaN.*ones(era5lonsz,era5latsz,6,nyr*climodaylen*6);end
    era5plevs=[500;700;850;925;975;1000];
    hrstart=1;everynthhour=4;

    %UTC-to-LST conversion notes
    %hr 1 = 00 UTC == 04 LST
    %hr 21 = 20 UTC = 00 LST
    %hr=1:4:21
    %relhr=round2(hr/4,1,'ceil');
    %localrelhr=relhr+1;
    %if localrelhr==7;localrelhr=1;end 

    for y=startyr:stopyr
        ystr=num2str(y);
        for mon=6:8
            if mon<=9;zm='0';else;zm='';end
            thislen=monlens(mon);
            fname=strcat(era5loc,'era5_pressurelevel_abudhabi_',ystr,zm,num2str(mon),'.nc');
            startvec=[1 1 1 1];countvec=[Inf Inf Inf thislen*(24/everynthhour)];stridevec=[1 1 1 everynthhour];
            tlevs(:,:,:,hrstart:hrstart+thislen*(24/everynthhour)-1)=ncread(fname,'t',startvec,countvec,stridevec)-273.15; %C
            qlevs(:,:,:,hrstart:hrstart+thislen*(24/everynthhour)-1)=ncread(fname,'q',startvec,countvec,stridevec).*1000; %g/kg
            if alsodouandv==1
            ulevs(:,:,:,hrstart:hrstart+thislen*(24/everynthhour)-1)=ncread(fname,'u',startvec,countvec,stridevec); %m/s
            vlevs(:,:,:,hrstart:hrstart+thislen*(24/everynthhour)-1)=ncread(fname,'v',startvec,countvec,stridevec); %m/s
            end

            hrstart=hrstart+thislen*(24/everynthhour);
        end
        if rem(y,5)==0;disp(y);disp(clock);end
    end

    %Now, note that first element of dim 4 is hr 1 = 00 UTC = 04 LST 
    %and last element is 00 LST

    %Rotate everything to normal spatial orientation (1 min for 20 years)
    tlevs_era5=permute(tlevs,[2 1 3 4]);clear tlevs;qlevs_era5=permute(qlevs,[2 1 3 4]);clear qlevs;
    if alsodouandv==1
    ulevs_era5=permute(ulevs,[2 1 3 4]);clear ulevs;vlevs_era5=permute(vlevs,[2 1 3 4]);clear vlevs;
    ulevs_era5=round2(ulevs_era5,0.01);vlevs_era5=round2(vlevs_era5,0.01); %extra precision not needed, and reduces array sizes pleasantly
    end

    %save(strcat(saveloc,'tlevs_era5'),'tlevs_era5');save(strcat(saveloc,'qlevs_era5'),'qlevs_era5');
end

if era5landreading_main==1
    %Read all data for JJA 2000-2020
    %Reminder: data was downloaded for MJJAS

    %Code here is arranged in a seemingly inefficient way but one that's
    %necessary to keep Matlab from crashing with the permute operation
    t_era5land=zeros(e5llatsz,e5llonsz,nyr*climodaylen*24);hrstart=1;
    for y=startyr:stopyr
        ystr=num2str(y);
        tmpvals=ncread(strcat(era5landloc,'t_',ystr,'.nc'),'t2m')-273.15;t_era5land(:,:,hrstart:hrstart+climodaylen*24-1)=permute(tmpvals(:,:,31*24+1:(31+climodaylen)*24),[2 1 3]);
        hrstart=hrstart+climodaylen*24;
    end
    save(strcat(saveloc,'tera5land.mat'),'t_era5land','-v7.3');

    q_era5land=zeros(e5llatsz,e5llonsz,nyr*climodaylen*24);hrstart=1;
    for y=startyr:stopyr
        ystr=num2str(y);
        tmpvals=calcqfromTd(ncread(strcat(era5landloc,'td_',ystr,'.nc'),'d2m')-273.15);
        q_era5land(:,:,hrstart:hrstart+climodaylen*24-1)=permute(tmpvals(:,:,31*24+1:(31+climodaylen)*24),[2 1 3]);
        hrstart=hrstart+climodaylen*24;
    end
    save(strcat(saveloc,'qera5land.mat'),'q_era5land','-v7.3');

    u_era5land=zeros(e5llatsz,e5llonsz,nyr*climodaylen*24);hrstart=1;
    for y=startyr:stopyr
        ystr=num2str(y);
        tmpvals=ncread(strcat(era5landloc,'u_',ystr,'.nc'),'u10');u_era5land(:,:,hrstart:hrstart+climodaylen*24-1)=permute(tmpvals(:,:,31*24+1:(31+climodaylen)*24),[2 1 3]);
        hrstart=hrstart+climodaylen*24;
    end
    save(strcat(saveloc,'uera5land.mat'),'u_era5land','-v7.3');clear u_era5land; %reload it only as necessary

    v_era5land=zeros(e5llatsz,e5llonsz,nyr*climodaylen*24);hrstart=1;
    for y=startyr:stopyr
        ystr=num2str(y);
        tmpvals=ncread(strcat(era5landloc,'v_',ystr,'.nc'),'v10');v_era5land(:,:,hrstart:hrstart+climodaylen*24-1)=permute(tmpvals(:,:,31*24+1:(31+climodaylen)*24),[2 1 3]);
        hrstart=hrstart+climodaylen*24;
    end
    save(strcat(saveloc,'vera5land.mat'),'v_era5land','-v7.3');clear v_era5land; %reload it only as necessary
end


if era5_twcalc_part1==1  
    %tw2m=NaN.*ones(size(t2m));
    %for k=1:size(t2m,3)
    %    tw2m(:,:,k)=calcwbt_daviesjones(t2m(:,:,k),psfc(:,:,k),q2m(:,:,k)./1000); %10 min PER YEAR
    %end

    %It's about 100x faster to calculate Tw in Python and then read it back in
    %Would be great to do this directly in Matlab, e.g. via pyrunfile, but
    %that is throwing out errors that seem difficult to remedy

    %Note: T in C, q in g/kg, psfc in Pa

    if strcmp(datasettodo,'era5_sl')
        arrt=t_era5;arrq=q_era5;arrp=psfc_era5;
    elseif strcmp(datasettodo,'era5_pl')
        arrt=tlevs_era5;arrq=qlevs_era5;arrp=zeros(size(arrt));
        for lev=1:size(arrp,3)
            arrp(:,:,lev,:)=100.*repmat(era5plevs(lev),[size(arrp,1) size(arrp,2) size(arrp,4)]);
        end
    elseif strcmp(datasettodo,'era5land')
        exist psfc_era5land;
        if ans==0
            psfc_era5land=NaN.*ones(e5llatsz,e5llonsz,size(psfc_era5,3));
            for d=1:size(psfc_era5,3)
                psfc_era5land(:,:,d)=interp2(lon_era5,lat_era5,squeeze(psfc_era5(:,:,d)),lon_era5land,lat_era5land);
            end
        end
        arrt=t_era5land;arrq=q_era5land;arrp=psfc_era5land;
    end
    calctw_matlabpythonhelper_part1; %then, need to run Python script manually as described in that file
end

if era5_twcalc_part2==1   
    myfname='twlevs'; %file prefix to use in saving
    calctw_matlabpythonhelper_part2; %requires having done calctw_matlabpythonhelper_part1 already

    %Output is tw_era5
end


if readmoreradiosondelevels==1
    rharmdir='/Volumes/ExternalDriveD/RHARM_radiosondes/';
    extsavedir='/Volumes/ExternalDriveF/Humidity_Trends/';
    savefilename='rharmradiosondearray_extralowlevels.mat';

    goallevs=[950;975;990];
    vnames={'station_name';'radiosonde_code';'sensor_model';'report_timestamp';'actual_time';'report_id';'longitude';'latitude';'height_of_station_above_sea_level';...
    'air_pressure';'air_temperature';'air_temperature_total_uncertainty';'relative_humidity';'relative_humidity_total_uncertainty';'wind_speed';'wind_from_direction';...
    'eastward_wind_component';'northward_wind_component';'geopotential_height'};
    
    numyrs=51;startyr=1970; %1970-2020 as a legacy to match what was done before in obs_humidity_trends_analysis
    x=1;save(strcat(extsavedir,savefilename),'x'); %to be able to later append whatever I want to this file

    readrharmradiosondedata(rharmdir,extsavedir,savefilename,goallevs,vnames,numyrs,1980,2019,startyr);

    %dims of standard output (e.g. tdata_rharm) are stn | year | month | day of month | hour of day (00 or 12 UTC) | pressure level (same order as in goallevs)
end


if reloadrharmradiosondes==1
    extsavedir='/Volumes/ExternalDriveF/Humidity_Trends/';
    f=load(strcat(extsavedir,'rharmradiosondearray.mat'));
    tdata_rharm=f.tdata;rhdata_rharm=f.rhdata;latlist_rharm=f.latlist;lonlist_rharm=f.lonlist;
    windspddata_rharm=f.windspddata;winddirdata_rharm=f.winddirdata;
    f=load(strcat(extsavedir,'rharmradiosondearray_extralowlevels.mat'));
    tdata_rharm_ell=f.tdata_rharm;rhdata_rharm_ell=f.rhdata_rharm;latlist_rharm_ell=f.latlist_rharm;lonlist_rharm_ell=f.lonlist_rharm;
    windspddata_rharm_ell=f.windspddata_rharm;winddirdata_rharm_ell=f.winddirdata_rharm;

    %Harmonize pressure levels
    %Investigation reveals extra-low-levels version has acquired an
        %additional 3 stations, but can be truncated to exactly match original
    %Pressure levels are now 200, 300, 500, 700, 850, 925, 950, 975, 990, 1000
    tdata_rharm_final=cat(6,tdata_rharm(:,:,:,:,:,1:6),tdata_rharm_ell(1:695,:,:,:,:,:),tdata_rharm(:,:,:,:,:,7));
    rhdata_rharm_final=cat(6,rhdata_rharm(:,:,:,:,:,1:6),rhdata_rharm_ell(1:695,:,:,:,:,:),rhdata_rharm(:,:,:,:,:,7));
    windspddata_rharm_final=cat(6,windspddata_rharm(:,:,:,:,:,1:6),windspddata_rharm_ell(1:695,:,:,:,:,:),windspddata_rharm(:,:,:,:,:,7));
    winddirdata_rharm_final=cat(6,winddirdata_rharm(:,:,:,:,:,1:6),winddirdata_rharm_ell(1:695,:,:,:,:,:),winddirdata_rharm(:,:,:,:,:,7));
    latlist_rharm_final=latlist_rharm;lonlist_rharm_final=lonlist_rharm;

    clear tdata_rharm;clear rhdata_rharm;clear windspddata_rharm;clear winddirdata_rharm;
    clear tdata_rharm_ell;clear rhdata_rharm_ell;clear windspddata_rharm_ell;clear winddirdata_rharm_ell;


    c=0;clear t_rharm;clear rh_rharm;
    for i=1:size(latlist_rharm,1)
        if latlist_rharm(i)>=lgmap_s && latlist_rharm(i)<=lgmap_n && lonlist_rharm(i)>=lgmap_w && lonlist_rharm(i)<=lgmap_e
            c=c+1;
            t_rharm(c,:,:,:,:,:)=tdata_rharm_final(i,:,:,:,:,:); %C
            rh_rharm(c,:,:,:,:,:)=rhdata_rharm_final(i,:,:,:,:,:); %percent
            windspd_rharm(c,:,:,:,:,:)=windspddata_rharm_final(i,:,:,:,:,:);
            winddir_rharm(c,:,:,:,:,:)=winddirdata_rharm_final(i,:,:,:,:,:);
            stnlats_rharm(c)=latlist_rharm(i);
            stnlons_rharm(c)=lonlist_rharm(i);
        end
    end
    pdata=permute(repmat([200;300;500;700;850;925;950;975;990;1000].*100,...
        [1 size(t_rharm,1) size(t_rharm,2) size(t_rharm,3) size(t_rharm,4) size(t_rharm,5)]),[2 3 4 5 6 1]);
    tw_rharm=calcwbt_daviesjones(t_rharm,pdata,rh_rharm,1);
    for stn=1:size(t_rharm,1)
        frac00avail(stn)=sum(sum(sum(~isnan(tw_rharm(stn,:,:,:,1,5)))))./(40*size(t_rharm,3)*30); 
        frac12avail(stn)=sum(sum(sum(~isnan(tw_rharm(stn,:,:,:,2,5)))))./(40*size(t_rharm,3)*30);
    end

    tw_rharm_jja20002020_00utc=squeeze(tw_rharm(:,31:51,6:8,:,1,:)); %dims are stns, years, months, days, then final dim is vert level -- ending with 1000
    tw_rharm_jja20002020_12utc=squeeze(tw_rharm(:,31:51,6:8,:,2,:));
    for stn=1:size(t_rharm,1)
        for lev=1:size(tw_rharm,6)
            tw_rharm_00utc_p95(stn,lev)=quantile(reshape(tw_rharm_jja20002020_00utc(stn,:,:,:,lev),[nyr*3*31 1]),0.95);
            tw_rharm_12utc_p95(stn,lev)=quantile(reshape(tw_rharm_jja20002020_12utc(stn,:,:,:,lev),[nyr*3*31 1]),0.95);
        end
    end
    invalid=tw_rharm_00utc_p95<=0;tw_rharm_00utc_p95(invalid)=NaN;
    invalid=tw_rharm_12utc_p95<=0;tw_rharm_12utc_p95(invalid)=NaN;
end





%For several variables: diurnal composite TS at Abu Dhabi for all days in JJA for all years
if diurnalcompositesetup==1
    exist nlat;if ans==0;nlat=size(tw_era5,1);nlon=size(tw_era5,2);nhrs=size(tw_era5,3);end

    %Get diurnal means (4 min for 20 years)
    %these are all based on UTC so far
    redodiurnalmeans=0;

    if redodiurnalmeans==1;tw2m_diurnalmean=squeeze(mean(tw2m_days,4));tw2m_diurnalmax=squeeze(max(tw2m_days,[],3));end

    tmp=load(strcat(saveloc,'pblhera5.mat'));pblh_era5=tmp.pblh_era5;clear tmp;
    pblh_days=reshape(pblh_era5,[nlat nlon 24 climodaylen*nyr]);clear pblh_era5;
    if redodiurnalmeans==1;pblh_diurnalmean=squeeze(mean(pblh_days,4));end

    exist t2m_days;if ans==0;exist t_era5;if ans==0;tmp=load(strcat(saveloc,'tera5.mat'));t_era5=tmp.t_era5;clear tmp;end;...
            t2m_days=reshape(t_era5,[nlat nlon 24 climodaylen*nyr]);clear t_era5;end
    if redodiurnalmeans==1;t2m_diurnalmean=squeeze(mean(t2m_days,4));t2m_diurnalmax=squeeze(max(t2m_days,[],3));end

    exist q2m_days;if ans==0;exist q_era5;if ans==0;tmp=load(strcat(saveloc,'qera5.mat'));q_era5=tmp.q_era5;clear tmp;end;q2m_days=reshape(q_era5,[nlat nlon 24 climodaylen*nyr]);clear q_era5;end
    if redodiurnalmeans==1;q2m_diurnalmean=squeeze(mean(q2m_days,4));end

    exist u_era5;if ans==0;tmp=load(strcat(saveloc,'uera5.mat'));u_era5=tmp.u_era5;clear tmp;tmp=load(strcat(saveloc,'vera5.mat'));v_era5=tmp.v_era5;clear tmp;end
    winddir10=winddirfromuandv(u_era5,v_era5);
    winddir10_rounded=round2(winddir10,30);clear winddir10; %necessary to get meaningful wind-dir modes
    winddir10_days=reshape(winddir10_rounded,[nlat nlon 24 climodaylen*nyr]);clear winddir10_rounded;
    if redodiurnalmeans==1;winddir10_diurnalmode=squeeze(mode(winddir10_days,4));end

    windspd10=sqrt(u_era5.^2+v_era5.^2);clear u_era5;clear v_era5;
    windspd10_days=reshape(windspd10,[nlat nlon 24 climodaylen*nyr]);clear windspd10;
    if redodiurnalmeans==1;windspd10_diurnalmean=squeeze(mean(windspd10_days,4));end

    if redodiurnalmeans==1;hi2m_diurnalmean=squeeze(mean(hi_days,4));hi2m_diurnalmax=squeeze(max(hi_days,[],3));end
    if redodiurnalmeans==1;utci2m_diurnalmean=squeeze(mean(utci_days,4,'omitnan'));utci2m_diurnalmax=squeeze(max(utci_days,[],3));end
    if redodiurnalmeans==1;wbgt2m_diurnalmean=squeeze(mean(wbgt_days,4,'omitnan'));wbgt2m_diurnalmax=squeeze(max(wbgt_days,[],3));end

    if redodiurnalmeans==1
        save(strcat(saveloc,'diurnalmeansera5.mat'),'tw2m_diurnalmean','pblh_diurnalmean','t2m_diurnalmean','q2m_diurnalmean','winddir10_diurnalmode','windspd10_diurnalmean',...
            'tw2m_diurnalmax','t2m_diurnalmax','hi2m_diurnalmean','hi2m_diurnalmax','utci2m_diurnalmean','utci2m_diurnalmax','wbgt2m_diurnalmean','wbgt2m_diurnalmax');
    end



    %Repeat diurnal-mean calculation but this time for Tw-defined hot days (any hour exceeding the local overall p95) only
    %Code is organized as it is for speed, to avoid having many big arrays unnecessarily early in the calculation
    %Runtime: about 1 hour
    pblhp95arr=cell(nlat,nlon);t2mp95arr=cell(nlat,nlon);tw2mp95arr=cell(nlat,nlon);q2mp95arr=cell(nlat,nlon);
    winddir10p95arr=cell(nlat,nlon);windspd10p95arr=cell(nlat,nlon);daysincl=zeros(nlat,nlon);
    hi2mp95arr=cell(nlat,nlon);utci2mp95arr=cell(nlat,nlon);wbgt2mp95arr=cell(nlat,nlon);

    for i=1:nlat;for j=1:nlon;c=0;for day=1:ndays;if sum(tw95whereabove(i,j,:,day))>=1;c=c+1;pblhp95arr{i,j}(c,:)=pblh_days(i,j,:,day);end;end;daysincl(i,j)=c;end;end
    save(strcat(saveloc,'diurnalmeansp95_pblh_era5.mat'),'pblhp95arr');clear pblhp95arr; %clear to speed up remainder of calculation

    for i=1:nlat;for j=1:nlon;c=0;for day=1:ndays;if sum(tw95whereabove(i,j,:,day))>=1;c=c+1;tw2mp95arr{i,j}(c,:)=tw2m_days(i,j,:,day);end;end;end;end
    save(strcat(saveloc,'diurnalmeansp95_tw_era5.mat'),'tw2mp95arr');clear tw2mp95arr;

    for i=1:nlat;for j=1:nlon;c=0;for day=1:ndays;if sum(tw95whereabove(i,j,:,day))>=1;c=c+1;t2mp95arr{i,j}(c,:)=t2m_days(i,j,:,day);end;end;end;end
    save(strcat(saveloc,'diurnalmeansp95_t_era5.mat'),'t2mp95arr');clear t2mp95arr;

    for i=1:nlat;for j=1:nlon;c=0;for day=1:ndays;if sum(tw95whereabove(i,j,:,day))>=1;c=c+1;q2mp95arr{i,j}(c,:)=q2m_days(i,j,:,day);end;end;end;end
    save(strcat(saveloc,'diurnalmeansp95_q_era5.mat'),'q2mp95arr');clear q2mp95arr;

    for i=1:nlat;for j=1:nlon;c=0;for day=1:ndays;if sum(tw95whereabove(i,j,:,day))>=1;c=c+1;winddir10p95arr{i,j}(c,:)=winddir10_days(i,j,:,day);end;end;end;end
    save(strcat(saveloc,'diurnalmeansp95_winddir_era5.mat'),'winddir10p95arr');clear winddir10p95arr;

    for i=1:nlat;for j=1:nlon;c=0;for day=1:ndays;if sum(tw95whereabove(i,j,:,day))>=1;c=c+1;windspd10p95arr{i,j}(c,:)=windspd10_days(i,j,:,day);end;end;end;end
    save(strcat(saveloc,'diurnalmeansp95_windspd_era5.mat'),'windspd10p95arr');clear windspd10p95arr;

    for i=1:nlat;for j=1:nlon;c=0;for day=1:ndays;if sum(tw95whereabove(i,j,:,day))>=1;c=c+1;hi2mp95arr{i,j}(c,:)=hi_days(i,j,:,day);end;end;end;end
    save(strcat(saveloc,'diurnalmeansp95_hi_era5.mat'),'hi2mp95arr');clear hi2mp95arr;

    for i=1:nlat;for j=1:nlon;c=0;for day=1:ndays;if sum(tw95whereabove(i,j,:,day))>=1;c=c+1;utci2mp95arr{i,j}(c,:)=utci_days(i,j,:,day);end;end;end;end
    save(strcat(saveloc,'diurnalmeansp95_utci_era5.mat'),'utci2mp95arr');clear utci2mp95arr;

    for i=1:nlat;for j=1:nlon;c=0;for day=1:ndays;if sum(tw95whereabove(i,j,:,day))>=1;c=c+1;wbgt2mp95arr{i,j}(c,:)=wbgt_days(i,j,:,day);end;end;end;end
    save(strcat(saveloc,'diurnalmeansp95_wbgt_era5.mat'),'wbgt2mp95arr');clear wbgt2mp95arr;
end


if era5p95diurnalcomposites==1
    exist pblhp95arr;if ans==0;f=load(strcat(saveloc,'diurnalmeansp95_pblh_era5.mat'));pblhp95arr=f.pblhp95arr;end
    exist tw2mp95arr;if ans==0;f=load(strcat(saveloc,'diurnalmeansp95_tw_era5.mat'));tw2mp95arr=f.tw2mp95arr;end
    exist t2mp95arr;if ans==0;f=load(strcat(saveloc,'diurnalmeansp95_t_era5.mat'));t2mp95arr=f.t2mp95arr;end
    exist q2mp95arr;if ans==0;f=load(strcat(saveloc,'diurnalmeansp95_q_era5.mat'));q2mp95arr=f.q2mp95arr;end
    exist winddir10p95arr;if ans==0;f=load(strcat(saveloc,'diurnalmeansp95_winddir_era5.mat'));winddir10p95arr=f.winddir10p95arr;end
    exist windspd10p95arr;if ans==0;f=load(strcat(saveloc,'diurnalmeansp95_windspd_era5.mat'));windspd10p95arr=f.windspd10p95arr;end
    exist hi2mp95arr;if ans==0;f=load(strcat(saveloc,'diurnalmeansp95_hi_era5.mat'));hi2mp95arr=f.hi2mp95arr;end
    exist utci2mp95arr;if ans==0;f=load(strcat(saveloc,'diurnalmeansp95_utci_era5.mat'));utci2mp95arr=f.utci2mp95arr;end
    exist wbgt2mp95arr;if ans==0;f=load(strcat(saveloc,'diurnalmeansp95_wbgt_era5.mat'));wbgt2mp95arr=f.wbgt2mp95arr;end

    pblhp95_diurnalmean=NaN.*ones(nlat,nlon,24);t2mp95_diurnalmean=NaN.*ones(nlat,nlon,24);tw2mp95_diurnalmean=NaN.*ones(nlat,nlon,24);
    q2mp95_diurnalmean=NaN.*ones(nlat,nlon,24);winddir10p95_diurnalmode=NaN.*ones(nlat,nlon,24);windspd10p95_diurnalmean=NaN.*ones(nlat,nlon,24);
    hi2mp95_diurnalmean=NaN.*ones(nlat,nlon,24);utci2mp95_diurnalmean=NaN.*ones(nlat,nlon,24);wbgt2mp95_diurnalmean=NaN.*ones(nlat,nlon,24);
    for i=1:nlat
        for j=1:nlon
            pblhp95_diurnalmean(i,j,:)=mean(pblhp95arr{i,j});
            t2mp95_diurnalmean(i,j,:)=mean(t2mp95arr{i,j});
            tw2mp95_diurnalmean(i,j,:)=mean(tw2mp95arr{i,j});
            q2mp95_diurnalmean(i,j,:)=mean(q2mp95arr{i,j});
            winddir10p95_diurnalmode(i,j,:)=mode(winddir10p95arr{i,j});
            windspd10p95_diurnalmean(i,j,:)=mean(windspd10p95arr{i,j});
            hi2mp95_diurnalmean(i,j,:)=mean(hi2mp95arr{i,j});
            utci2mp95_diurnalmean(i,j,:)=mean(utci2mp95arr{i,j},'omitnan');
            wbgt2mp95_diurnalmean(i,j,:)=mean(wbgt2mp95arr{i,j},'omitnan');
        end
    end

    northwind=winddir10p95_diurnalmode==360;winddir10p95_diurnalmode(northwind)=0;
    northwind=winddir10_diurnalmode==360;winddir10_diurnalmode(northwind)=0;

    %Further, for smoothness, convert wind dirs to 3-hour rolling means
    winddir10p95_diurnalmode_3hrrolling=NaN.*ones(size(winddir10p95_diurnalmode));
    winddir10_diurnalmode_3hrrolling=NaN.*ones(size(winddir10_diurnalmode));
    for hr=1:24
        if hr==1
            tmp=cat(3,winddir10p95_diurnalmode(:,:,24),winddir10p95_diurnalmode(:,:,1),winddir10p95_diurnalmode(:,:,2));
            winddir10p95_diurnalmode_3hrrolling(:,:,hr)=mean(tmp,3);
            tmp=cat(3,winddir10_diurnalmode(:,:,24),winddir10_diurnalmode(:,:,1),winddir10_diurnalmode(:,:,2));
            winddir10_diurnalmode_3hrrolling(:,:,hr)=mean(tmp,3);
        elseif hr==24
            tmp=cat(3,winddir10p95_diurnalmode(:,:,23),winddir10p95_diurnalmode(:,:,24),winddir10p95_diurnalmode(:,:,1));
            winddir10p95_diurnalmode_3hrrolling(:,:,hr)=mean(tmp,3);
            tmp=cat(3,winddir10_diurnalmode(:,:,23),winddir10_diurnalmode(:,:,24),winddir10_diurnalmode(:,:,1));
            winddir10_diurnalmode_3hrrolling(:,:,hr)=mean(tmp,3);
        else
            winddir10p95_diurnalmode_3hrrolling(:,:,hr)=mean(winddir10p95_diurnalmode(:,:,hr-1:hr+1),3);
            winddir10_diurnalmode_3hrrolling(:,:,hr)=mean(winddir10_diurnalmode(:,:,hr-1:hr+1),3);
        end
    end

    save(strcat(saveloc,'era5diurnalmeans'),'pblhp95_diurnalmean','t2mp95_diurnalmean','tw2mp95_diurnalmean','q2mp95_diurnalmean',...
        'winddir10p95_diurnalmode','windspd10p95_diurnalmean','hi2mp95_diurnalmean','utci2mp95_diurnalmean','wbgt2mp95_diurnalmean');
end


if makefig3_revised==1
    %Diurnal composite for several points near Abu Dhabi
    %As above, 1: Abu Dhabi city, 2: Persian Gulf waters, 3: rural coast, 4: inland oasis
    %LST is UTC+4, so shift diurnal-mean arrays to reflect this
    %Row 1: Tw, UTCI
    %Row 2: HI, PBLH
    %Row 3: T, q
    %Row 4: wind dir, wind spd
    figure(777);clf;
    xticks=1:4:24;xtl={'00','04','08','12','16','20'};


    %Get means for each of the 3 new regions in myregs_era5: coastal water, coastal land, inland
    arrsall={tw2m_diurnalmean;utci2m_diurnalmean;wbgt2m_diurnalmean;hi2m_diurnalmean;...
        t2m_diurnalmean;q2m_diurnalmean;pblh_diurnalmean;windspd10_diurnalmean;winddir10_diurnalmode};
    arrs95={tw2mp95_diurnalmean;utci2mp95_diurnalmean;wbgt2mp95_diurnalmean;hi2mp95_diurnalmean;...
        t2mp95_diurnalmean;q2mp95_diurnalmean;pblhp95_diurnalmean;windspd10p95_diurnalmean;winddir10p95_diurnalmode};

    for idx=1:length(arrsall)
        arrall=arrsall{idx};arrallholder=cell(3,1);arrallregmean=NaN.*ones(3,24);
        arr95=arrs95{idx};arr95holder=cell(3,1);arr95regmean=NaN.*ones(3,24);
        if idx==length(arrsall) %wind dir
            uallholder=cell(3,1);vallholder=cell(3,1);u95holder=cell(3,1);v95holder=cell(3,1);
        end
        c=zeros(3,1); %number of new regions
        for i=1:nlon
            for j=1:nlat
                r=NaN;
                if myregs_era5(i,j)==1 %coastal water
                    r=1;
                elseif myregs_era5(i,j)==2 %coastal land
                    r=2;
                elseif myregs_era5(i,j)==3 %inland
                    r=3;
                end
                if ~isnan(r)
                    c(r)=c(r)+1;
                    if idx<length(arrsall) %everything but wind dir
                        arrallholder{r}(c(r),:)=squeeze(arrall(i,j,:));
                        arr95holder{r}(c(r),:)=squeeze(arr95(i,j,:));
                    else %wind dir
                        for hr=1:24
                            [ucompall(hr),vcompall(hr)]=uandvfromwinddirandspeed(arrall(i,j,hr),windspd10_diurnalmean(i,j,hr));
                            [ucomp95(hr),vcomp95(hr)]=uandvfromwinddirandspeed(arr95(i,j,hr),windspd10p95_diurnalmean(i,j,hr));
                        end
                        uallholder{r}(c(r),:)=ucompall;u95holder{r}(c(r),:)=ucomp95;
                        vallholder{r}(c(r),:)=vcompall;v95holder{r}(c(r),:)=vcomp95;
                    end
                end
            end
        end
        if idx<length(arrsall) %for everything but wind dir
            for r=1:3
                arrallregmean(r,:)=mean(arrallholder{r});
                arrallregmeanminus1stdev(r,:)=mean(arrallholder{r})-std(arrallholder{r});
                arrallregmeanplus1stdev(r,:)=mean(arrallholder{r})+std(arrallholder{r});

                arr95regmean(r,:)=mean(arr95holder{r});
                arr95regmeanminus1stdev(r,:)=mean(arr95holder{r})+std(arr95holder{r});
                arr95regmeanplus1stdev(r,:)=mean(arr95holder{r})-std(arr95holder{r});
            end
            arrsall_mean{idx}=arrallregmean;arrsall_meanminus1stdev{idx}=arrallregmeanminus1stdev;arrsall_meanplus1stdev{idx}=arrallregmeanplus1stdev;
            arrs95_mean{idx}=arr95regmean;arrs95_meanminus1stdev{idx}=arr95regmeanminus1stdev;arrs95_meanplus1stdev{idx}=arr95regmeanplus1stdev;
        else %wind dir
            for r=1:3
                uallregmean(r,:)=mean(uallholder{r});vallregmean(r,:)=mean(vallholder{r});
                uallregmeanminus1stdev(r,:)=mean(uallholder{r})-std(uallholder{r});uallregmeanplus1stdev(r,:)=mean(uallholder{r})+std(uallholder{r});
                vallregmeanminus1stdev(r,:)=mean(vallholder{r})-std(vallholder{r});vallregmeanplus1stdev(r,:)=mean(vallholder{r})+std(vallholder{r});

                u95regmean(r,:)=mean(u95holder{r});v95regmean(r,:)=mean(v95holder{r});
                u95regmeanminus1stdev(r,:)=mean(u95holder{r})-std(u95holder{r});u95regmeanplus1stdev(r,:)=mean(u95holder{r})+std(u95holder{r});
                v95regmeanminus1stdev(r,:)=mean(v95holder{r})-std(v95holder{r});v95regmeanplus1stdev(r,:)=mean(v95holder{r})+std(v95holder{r});
            end
        end
    end

    
    titles={'Tw';'UTCI';'WBGT';'Heat Index';'T';'q';'PBL height'};
    ylabs={strcat(char(176),'C');strcat(char(176),'C');strcat(char(176),'C');strcat(char(176),'C');...
        strcat(char(176),'C');'g/kg';'m';'';'m/s'};
    lefts=[0.03;0.35;0.67;0.03;0.35;0.67;0.03;0.35;0.67];bottoms=[0.68;0.68;0.68;0.35;0.35;0.35;0.02;0.02;0.02];wwidth=0.32;hheight=0.32;

    %For each, plot 00 LST as first element
    for loop=1:3 %regions
        cr=newregcolors(loop,:);
        difffrommax=max(cr)-cr;cr_pale=cr+0.4*difffrommax;

        for sp=1:7
            subplot(3,3,sp);
            arrall_mean=[squeeze(arrsall_mean{sp}(loop,21:24))';squeeze(arrsall_mean{sp}(loop,1:20))'];
                arrall_minus1stdev=[squeeze(arrsall_meanminus1stdev{sp}(loop,21:24))';squeeze(arrsall_meanminus1stdev{sp}(loop,1:20))'];
                arrall_plus1stdev=[squeeze(arrsall_meanplus1stdev{sp}(loop,21:24))';squeeze(arrsall_meanplus1stdev{sp}(loop,1:20))'];
            arr95_mean=[squeeze(arrs95_mean{sp}(loop,21:24))';squeeze(arrs95_mean{sp}(loop,1:20))'];
                arr95_minus1stdev=[squeeze(arrs95_meanminus1stdev{sp}(loop,21:24))';squeeze(arrs95_meanminus1stdev{sp}(loop,1:20))'];
                arr95_plus1stdev=[squeeze(arrs95_meanplus1stdev{sp}(loop,21:24))';squeeze(arrs95_meanplus1stdev{sp}(loop,1:20))'];

            %Include error bars only at 00, 04, 08, etc. LST
            arrall_neg=arrall_minus1stdev-arrall_mean;arrall_pos=arrall_plus1stdev-arrall_mean;
            arr95_neg=arr95_minus1stdev-arr95_mean;arr95_pos=arr95_plus1stdev-arr95_mean;
            for i=1:24
                if rem(i,4)~=1;arrall_neg(i)=NaN;arrall_pos(i)=NaN;end
                if rem(i,4)~=2;arr95_neg(i)=NaN;arr95_pos(i)=NaN;end
            end
                
            %plot(arrall_mean,'linestyle',':','linewidth',1.7,'color',cr);hold on;
            %plot(arr95_mean,'linestyle','-','linewidth',1.5,'color',cr);
            p=errorbar(1:24,arrall_mean,arrall_neg,arrall_pos,'linestyle',':','linewidth',1.7,'color',cr_pale,'CapSize',0);hold on;p.Bar.LineStyle='dotted';
            p=errorbar(1:24,arr95_mean,arr95_neg,arr95_pos,'linestyle','-','linewidth',1.7,'color',cr,'CapSize',0);
            
            xlim([1 24]);maketitle(titles{sp});set(gca,'fontweight','bold');set(gca,'xtick',xticks,'xticklabel',xtl);
            ylabel(ylabs{sp},'fontweight','bold','fontname','arial','fontsize',11);
            if sp==1;ylim([20 32]);elseif sp==2;ylim([24 55]);elseif sp==3;ylim([23 35]);elseif sp==4;ylim([30 52]);elseif sp==5;ylim([28 45]);...
            elseif sp==6;ylim([8 28]);elseif sp==7;ylim([0 2800]);end
            t=text(-0.05,0.55,splabels{sp},'units','normalized');set(t,'fontweight','bold','fontname','arial','fontsize',12);
            if sp==7;xlabel('Time (LST)','fontweight','bold','fontname','arial','fontsize',12);end
        end
        
        %Wind direction
        sp=8;subplot(3,3,sp);hold on;
            %For clarity, plot wind vectors
            the_yorder=[1;2;3];

            uall_lst=[uallregmean(loop,21:24)';uallregmean(loop,1:20)'];u95_lst=[u95regmean(loop,21:24)';u95regmean(loop,1:20)'];
            vall_lst=[vallregmean(loop,21:24)';vallregmean(loop,1:20)'];v95_lst=[v95regmean(loop,21:24)';v95regmean(loop,1:20)'];

            for hr=1:4:24
                xbase=(hr-0.5)/23;ybase=1.05-0.25*the_yorder(loop);

                ucomp=uall_lst(hr);vcomp=vall_lst(hr);
                ucomp_norm=ucomp./(8.*sqrt(ucomp.^2+vcomp.^2));vcomp_norm=vcomp./(8.*sqrt(ucomp.^2+vcomp.^2)); %arbritrary shrinkage to look good
                anArrow=annotation('arrow','color',cr_pale,'linewidth',1.5,'linestyle',':','headlength',7,'headwidth',7,'headstyle','deltoid');anArrow.Parent=gca;
                %This part is a bit hacky, but c'est la vie
                %In a 3x2 grid, by default each subplot is about 3cm tall by 8cm wide on the screen
                %Therefore, the x component needs to be multiplied by 3/8 so that e.g. 315 deg truly looks like a NW wind, 337.5 is NNW, etc
                anArrow.Position=[xbase,ybase,ucomp_norm*3/8,vcomp_norm];

                ucomp=u95_lst(hr);vcomp=v95_lst(hr);
                ucomp_norm=ucomp./(8.*sqrt(ucomp.^2+vcomp.^2));vcomp_norm=vcomp./(8.*sqrt(ucomp.^2+vcomp.^2));
                anArrow=annotation('arrow','color',cr,'linewidth',1.5,'linestyle','-','headlength',7,'headwidth',7,'headstyle','deltoid');anArrow.Parent=gca;
                anArrow.Position=[xbase,ybase-0.05,ucomp_norm*3/8,vcomp_norm];
            end
            maketitle('Wind Dir');set(gca,'fontweight','bold');set(gca,'xtick',(xticks-1)./23,'xticklabel',xtl);
            set(gca,'ytick',0.3:0.25:0.8,'yticklabel',{'Inland';'Coast';'Gulf'});ax=gca;ax.Clipping='off';box on;
            t=text(-0.05,0.55,splabels{sp},'units','normalized');set(t,'fontweight','bold','fontname','arial','fontsize',12);
            xlabel('Time (LST)','fontweight','bold','fontname','arial','fontsize',12);

        %Wind speed
        sp=9;subplot(3,3,sp);hold on;
            arrall_mean=[squeeze(arrsall_mean{sp-1}(loop,21:24))';squeeze(arrsall_mean{sp-1}(loop,1:20))'];
                arrall_minus1stdev=[squeeze(arrsall_meanminus1stdev{sp-1}(loop,21:24))';squeeze(arrsall_meanminus1stdev{sp-1}(loop,1:20))'];
                arrall_plus1stdev=[squeeze(arrsall_meanplus1stdev{sp-1}(loop,21:24))';squeeze(arrsall_meanplus1stdev{sp-1}(loop,1:20))'];
            arr95_mean=[squeeze(arrs95_mean{sp-1}(loop,21:24))';squeeze(arrs95_mean{sp-1}(loop,1:20))'];
                arr95_minus1stdev=[squeeze(arrs95_meanminus1stdev{sp-1}(loop,21:24))';squeeze(arrs95_meanminus1stdev{sp-1}(loop,1:20))'];
                arr95_plus1stdev=[squeeze(arrs95_meanplus1stdev{sp-1}(loop,21:24))';squeeze(arrs95_meanplus1stdev{sp-1}(loop,1:20))'];

            %Include error bars only at 00, 04, 08, etc. LST
            arrall_neg=arrall_minus1stdev-arrall_mean;arrall_pos=arrall_plus1stdev-arrall_mean;
            arr95_neg=arr95_minus1stdev-arr95_mean;arr95_pos=arr95_plus1stdev-arr95_mean;
            for i=1:24
                if rem(i,4)~=1;arrall_neg(i)=NaN;arrall_pos(i)=NaN;end
                if rem(i,4)~=2;arr95_neg(i)=NaN;arr95_pos(i)=NaN;end
            end

            %plot(arrall_mean,'linestyle',':','linewidth',1.7,'color',cr);hold on;
            %plot(arr95_mean,'linestyle','-','linewidth',1.5,'color',cr);
            p=errorbar(1:24,arrall_mean,arrall_neg,arrall_pos,'linestyle',':','linewidth',1.7,'color',cr_pale,'CapSize',0);hold on;p.Bar.LineStyle='dotted';
            p=errorbar(1:24,arr95_mean,arr95_neg,arr95_pos,'linestyle','-','linewidth',1.7,'color',cr,'CapSize',0);

            xlim([1 24]);maketitle('Wind Speed');set(gca,'fontweight','bold');set(gca,'xtick',xticks,'xticklabel',xtl);
            ylabel(ylabs{sp},'fontweight','bold','fontname','arial','fontsize',11);
            ylim([1.7 6.2]);box on;
            t=text(-0.05,0.55,splabels{sp},'units','normalized');set(t,'fontweight','bold','fontname','arial','fontsize',12);
            xlabel('Time (LST)','fontweight','bold','fontname','arial','fontsize',12);
    end
    %Add a few final legend elements, then save
    a=annotation('line',[0.27 0.35],[0.03 0.03]);set(a,'linewidth',2,'color',newregcolors(1,:));
    t=text(-0.786,-0.19,'Gulf','units','normalized');set(t,'fontweight','bold','fontname','arial','fontsize',12,'color',newregcolors(1,:));
    a=annotation('line',[0.46 0.54],[0.03 0.03]);set(a,'linewidth',2,'color',newregcolors(2,:));
    t=text(-0.343,-0.19,'Coast','units','normalized');set(t,'fontweight','bold','fontname','arial','fontsize',12,'color',newregcolors(2,:));
    a=annotation('line',[0.65 0.73],[0.03 0.03]);set(a,'linewidth',2,'color',newregcolors(3,:));
    t=text(0.10,-0.19,'Inland','units','normalized');set(t,'fontweight','bold','fontname','arial','fontsize',12,'color',newregcolors(3,:));
    set(gcf,'color','w');
    curpart=1;highqualityfiguresetup;figname='fig3_rev';curpart=2;highqualityfiguresetup;


    makefigr1=0; %default 0, as this is just for reviewer response
    if makefigr1==1
        tall=arrsall_mean{5};tdall=calcTdfromq(arrsall_mean{6});
        rhall=calcrhfromTandTd(tall,tdall);
        t95=arrs95_mean{5};td95=calcTdfromq(arrs95_mean{6});
        rh95=calcrhfromTandTd(t95,td95);

        figure(45);clf;hold on;
        for loop=1:3
            cr=newregcolors(loop,:);
            plot([squeeze(rhall(loop,21:24))';squeeze(rhall(loop,1:20))'],'linestyle',':','linewidth',1.7,'color',cr);hold on;
            plot([squeeze(rh95(loop,21:24))';squeeze(rh95(loop,1:20))'],'linestyle','-','linewidth',1.5,'color',cr);
            xlim([1 24]);
            title('RH','fontweight','bold','fontname','arial','fontsize',14);
            set(gca,'fontweight','bold');set(gca,'xtick',xticks,'xticklabel',xtl);
            ylabel('Percent','fontweight','bold','fontname','arial','fontsize',12);
        end
        xlabel('Time (LST)','fontweight','bold','fontname','arial','fontsize',12);
        set(gcf,'color','w');
        curpart=1;highqualityfiguresetup;figname='figr1_rh';curpart=2;highqualityfiguresetup;
    end
end




if makefig2_asboxplots==1
    figure(1002);clf;c=0;figname='fig2_boxplots';
    lefts=[0.05 0.52;0.05 0.52;0.05 0.52];bottoms=[0.69 0.69;0.385 0.385;0.08 0.08];wwidth=0.42;hheight=0.26;
    for vl=1:3
        if vl==1
            allstn=allstntw;allgridpt=allgridpttw;ylab=strcat('Tw (',char(176),'C)');
        elseif vl==2
            allstn=allstnt;allgridpt=allgridptt;ylab=strcat('T (',char(176),'C)');
        elseif vl==3
            allstn=allstnq;allgridpt=allgridptq;ylab='q (g/kg)';
        end
        for reg=2:3
            subplot(3,2,vl*2-2+reg-1);hold on;c=c+1;
            stndatabyhr_tmp=reshape(allstn{reg},[size(allstn{reg},1)*size(allstn{reg},2) size(allstn{reg},3)]);
            gridptdatabyhr_tmp=reshape(allgridpt{reg},[size(allgridpt{reg},1)*size(allgridpt{reg},2) size(allgridpt{reg},3)]);

            %Convert from UTC to LST
            stndatabyhr=cat(2,stndatabyhr_tmp(:,21:24),stndatabyhr_tmp(:,1:20));
            gridptdatabyhr=cat(2,gridptdatabyhr_tmp(:,21:24),gridptdatabyhr_tmp(:,1:20));


            data=[];groups=[];stnsz=size(stndatabyhr,1);gridptsz=size(gridptdatabyhr,1);
            for hr=1:2:24
                data=cat(1,data,stndatabyhr(:,hr));data=cat(1,data,gridptdatabyhr(:,hr));
                groups=cat(1,groups,(hr-1)*2+ones(stnsz,1));groups=cat(1,groups,(hr-1)*2+1+ones(gridptsz,1));
            end
            positions=[1 1.7 3 3.7 5 5.7 7 7.7 9 9.7 11 11.7 13 13.7 15 15.7 17 17.7 19 19.7 21 21.7 23 23.7];
            bh=boxplot(data,groups,'positions',positions,'symbol','','width',0.4); %with outliers removed

            set(bh,'linewidth',1.5);
            ax=gca;
            for i=25:168
                if rem(i,2)==0
                    ax.Children.Children(i).Color=colors('sirkes_dark gold'); %stations
                else
                    ax.Children.Children(i).Color='k'; %ERA5 gridpts
                end
                ax.Children.Children(i).LineWidth=1;
            end
            %Delete whisker endcaps for station data, to help distinguish groups
            allTags=get(bh,'tag');isCap=contains(allTags,'Adjacent');considerdeleting=bh(isCap);
            %for i=1:4:length(considerdeleting);delete(considerdeleting(i));end
            %for i=2:4:length(considerdeleting);delete(considerdeleting(i));end
            for i=1:length(considerdeleting);delete(considerdeleting(i));end

            %Fill boxes with color for better visual impression
            colorarr=[colors('gray');colors('sirkes_dark gold')];colorarr=repmat(colorarr,[12 1]);
            h=findobj(gca,'Tag','Box');
            for j=1:length(h);patch(get(h(j),'XData'),get(h(j),'YData'),colorarr(j,:),'FaceAlpha',0.5);end

            %set(gca,'ylim',[round2(min(data),1,'floor') round2(max(data),1,'ceil')]);
            if vl==1;set(gca,'ylim',[11 34]);elseif vl==2;set(gca,'ylim',[19 50]);elseif vl==3;set(gca,'ylim',[0 32]);end
            set(gca,'xtick',[1:4:24],'xticklabel',{'00';'04';'08';'12';'16';'20'});
            set(gca,'fontweight','bold','fontname','arial','fontsize',11);
            
            t=text(-0.045,0.53,splabels{c},'units','normalized');set(t,'fontweight','bold','fontname','arial','fontsize',12);
            if vl==3;xlabel('Time (LST)','fontweight','bold','fontname','arial','fontsize',12);end
            if reg==2;ylabel(ylab,'fontweight','bold','fontname','arial','fontsize',12);end
            if reg==2;titletext='Coast';elseif reg==3;titletext='Inland';end
            if vl==1;title(titletext,'fontweight','bold','fontname','arial','fontsize',12);end

            otherapproach=0;
            if otherapproach==1
                clear A;A{1}=stndatabyhr;A{2}=gridptdatabyhr;
                boxplotGroup(A,'whisker',3);
    
                %To improve cleanliness of plot, remove boxplots for hours 01, 03, 05 LST -- #2, 4, 6, etc.
                ax=gca;
                for bp=1:2 %each boxplot group
                    bh=ax.Children(bp).Children;
                    for hrtoremove=2:2:24
                        idxtoremove=25-hrtoremove;
                        for elem=1:7
                            bh(24*(elem-1)+idxtoremove).Color='w';
                            if elem==1;bh(24*(elem-1)+idxtoremove).Marker='none';end
                        end
                    end
                    for hr=1:24;bh(hr).MarkerEdgeColor='k';bh(hr).Size=4;end %instead of red and 6
                end
            end

            set(gca,'position',[lefts(vl,reg-1) bottoms(vl,reg-1) wwidth hheight]);
        end
    end
    set(gcf,'color','w');
    curpart=1;highqualityfiguresetup;curpart=2;highqualityfiguresetup;
end


if calcspreads_era5==1
    %Refine hour-of-max-exceedance-probability figure by only plotting 'central hours': points where
    %>=50% of the p95 values occur within a 6-hour period centered on that
    %hour (in other words, distribution is reasonably peaky)
    %Runtime 1.5 min
    centralhours_era5=NaN.*ones(nlon,nlat);hrofmaxp95prob_tw_era5=NaN.*ones(nlon,nlat);
    for i=1:nlon
        for j=1:nlat
            %95th pctile across all hours at this gridcell
            for hr=1:24;datathishr=squeeze(tw2m_days(i,j,hr,:));probp95_tw_era5(hr)=sum(datathishr>quantile(tw_era5(i,j,:),0.95))/size(tw2m_days,4);end
            [~,hrofmaxp95prob_tw_era5(i,j)]=max(probp95_tw_era5);
            %Do at least 50% of p95 values occur within 3 of the peak?
            hrmax=hrofmaxp95prob_tw_era5(i,j);
            if hrmax>=22
                tosum=[probp95_tw_era5(hrmax-3:24) probp95_tw_era5(1:hrmax-21)];
            elseif hrmax<=3
                tosum=[probp95_tw_era5(hrmax+21) probp95_tw_era5(1:hrmax+3)];
            else
                tosum=[probp95_tw_era5(hrmax-3:hrmax+3)];
            end
            sum_within3=sum(tosum);
            sum_all=sum(probp95_tw_era5);
            frac_within3=round2(sum_within3/sum_all,0.05);
            if frac_within3>=minfrac
                centralhours_era5(i,j)=hrmax;
            end
        end
    end
    a=centralhours_era5+tzadj;a(a>=24)=a(a>=24)-24;centralhours_local_era5=a;
    %figure(10);clf;imagescnan(centralhours_local');colorbar;

    %Analogously, define 'central wind dirs' associated with p95 Tw -- at
    %least 50% within a range of 90 deg
    %Note on these requirements:
        %-requiring 50% of extreme Tw's to fall within 25% of possible values: provides quick glimpse across region
        %-requiring 50% of extreme Tw's to fall within 25% of actual obs: also separates mean from extreme conditions
    centralwinddirs_era5=NaN.*ones(nlon,nlat);
    exist winddir10_days;
    if ans==0
        exist u_era5;if ans==0;tmp=load(strcat(saveloc,'uera5.mat'));u_era5=tmp.u_era5;clear tmp;tmp=load(strcat(saveloc,'vera5.mat'));v_era5=tmp.v_era5;clear tmp;end
        winddir10=winddirfromuandv(u_era5,v_era5);clear u_era5;clear v_era5;
        winddir10_rounded=round2(winddir10,30);clear winddir10; %necessary to get meaningful wind-dir modes
        winddir10_days=reshape(winddir10_rounded,[nlat nlon 24 climodaylen*nyr]);clear winddir10_rounded;
    end
    for i=1:nlon
        for j=1:nlat
            thisp95=quantile(tw_era5(i,j,:),0.95);c=0;
            for hr=1:24
                for day=1:size(tw2m_days,4)
                    if tw2m_days(i,j,hr,day)>thisp95
                        c=c+1;
                        winddirvec(c)=winddir10_days(i,j,hr,day);
                    end
                end
            end
            peakdir=mode(winddirvec);
            frac_within45deg=round2(sum(abs(winddirvec-peakdir)<=45)/c,0.05);
            if frac_within45deg>=minfrac
                centralwinddirs_era5(i,j)=peakdir;
            end
        end
    end

    %Repeat for PBL height (only slightly interesting; wind spd not at all)
    dothis=0;
    if dothis==1
        centralpblhs_era5=NaN.*ones(nlon,nlat);

        exist pblh_days;
        if ans==0
            tmp=load(strcat(saveloc,'pblhera5.mat'));pblh_era5=tmp.pblh_era5;clear tmp;
            pblh_days=reshape(pblh_era5,[nlat nlon 24 climodaylen*nyr]);clear pblh_era5;
        end
        for i=1:nlon
            for j=1:nlat
                thisp95=quantile(tw_era5(i,j,:),0.95);c=0;
                for hr=1:24
                    for day=1:climodaylen
                        if tw2m_days(i,j,hr,day)>thisp95
                            c=c+1;
                            pblhvec(c)=pblh_days(i,j,hr,day);
                        end
                    end
                end
                peakpblh=mean(pblhvec);
    
                allpblhs_thispt=squeeze(pblh_era5(i,j,:));continueon=1;
                for tryout=10:10:5000
                    if continueon==1
                        fracpblhswithin=sum(abs(allpblhs_thispt-peakpblh)<=tryout)/nhrs;
                        if fracpblhswithin>=0.25;continueon=0;disttouse=tryout;end
                    end
                end
                frac_quarter=sum(abs(pblhvec-peakpblh)<=disttouse)/c;
                if frac_quarter>=minfrac
                    centralpblhs_era5(i,j)=peakpblh;
                end
            end
        end
        clear pblh_days;
    end

    if minfrac==0.5 %default
        save(strcat(saveloc,'centralvarsera5.mat'),'centralhours_local_era5','centralwinddirs_era5','centralpblhs_era5');
    elseif minfrac==0.4 %other typical choice
        save(strcat(saveloc,'centralvarsera5_minfrac0pt4.mat'),'centralhours_local_era5','centralwinddirs_era5','centralpblhs_era5');
    end
end



if domultivartimeseries==1
    %Figure 2: Timeseries of various variables at Abu Dhabi Intl Airport,
    %from various sources incl WRF

    %black: weather station
    %green: ERA5
    %orange: WRF no-urban run
    %red: WRF latest run

    dayofsummer=81;hourofday=11;
    hheight=(hw_year-startyr)*climodaylen*24+1+(dayofsummer-1)*24+hourofday-1;

    firstday=365*(hw_year-startyr)+hw_firstdoy;lastday=365*(hw_year-startyr)+hw_lastdoy;

    figure(85);clf;
    subplot(3,2,1);plot(reshape(squeeze(subdailytw_all(7,firstday:lastday,:))',[numhwdays*24 1]),'k','linewidth',1.5);maketitle('Wet-Bulb Temp.');hold on;
        xlim([1 numhwdays*24]);set(gca,'xtick',1:24:numhwdays*24,'xticklabel',firstday:lastday);set(gca,'fontweight','bold','fontname','arial','fontsize',11);
        %Add WRF output
        thists=cat(2,NaN.*ones(1,96),tw_nourban(1,:,besti,bestj),NaN.*ones(1,48))';plot(thists,'color',colors('orange'),'linewidth',1.5); %no-urban run
        thists=cat(2,NaN.*ones(1,96),tw_withdamp(1,:,besti,bestj),NaN.*ones(1,48))';plot(thists,'color',colors('red'),'linewidth',1.5); %latest run
        %Add ERA5
        thists=squeeze(tw_era5(besti_era5,bestj_era5,1825:1825+167));plot(thists,'color',colors('medium green'),'linewidth',1.5);

    %PBLH
    subplot(3,2,2);
    %WRF
    thists=cat(2,NaN.*ones(1,96),pblh_nourban(1,:,besti,bestj),NaN.*ones(1,48))';plot(thists,'color',colors('orange'),'linewidth',1.5); %no-urban run
        maketitle('PBL Height');hold on;
        xlim([1 numhwdays*24]);set(gca,'xtick',1:24:numhwdays*24,'xticklabel',firstday:lastday);set(gca,'fontweight','bold','fontname','arial','fontsize',11);  
    thists=cat(2,NaN.*ones(1,96),pblh_withdamp(1,:,besti,bestj),NaN.*ones(1,48))';plot(thists,'color',colors('red'),'linewidth',1.5); %latest run
    %ERA5
    thists=squeeze(pblhfinal_era5(besti_era5,bestj_era5,1825:1825+167));plot(thists,'color',colors('medium green'),'linewidth',1.5);
    

    %T
    subplot(3,2,3);plot(reshape(squeeze(subdailyt_all(7,firstday:lastday,:))',[numhwdays*24 1]),'k','linewidth',1.5);maketitle('Temperature');hold on;
        xlim([1 numhwdays*24]);set(gca,'xtick',1:24:numhwdays*24,'xticklabel',firstday:lastday);set(gca,'fontweight','bold','fontname','arial','fontsize',11);
        thists=cat(2,NaN.*ones(1,96),t2_nourban(1,:,besti,bestj),NaN.*ones(1,48))';plot(thists,'color',colors('orange'),'linewidth',1.5);
        thists=cat(2,NaN.*ones(1,96),t2_withdamp(1,:,besti,bestj),NaN.*ones(1,48))';plot(thists,'color',colors('red'),'linewidth',1.5); %latest run
        thists=squeeze(tfinal_era5(besti_era5,bestj_era5,1825:1825+167));plot(thists,'color',colors('medium green'),'linewidth',1.5);

    subplot(3,2,4);plot(reshape(squeeze(subdailytd_all(7,firstday:lastday,:))',[numhwdays*24 1]),'k','linewidth',1.5);maketitle('Dewpoint Temp.');hold on;
        xlim([1 numhwdays*24]);set(gca,'xtick',1:24:numhwdays*24,'xticklabel',firstday:lastday);set(gca,'fontweight','bold','fontname','arial','fontsize',11);
        thists=cat(2,NaN.*ones(1,96),calcTdfromq(q2_nourban(1,:,besti,bestj).*1000),NaN.*ones(1,48))';plot(thists,'color',colors('orange'),'linewidth',1.5);
        thists=cat(2,NaN.*ones(1,96),calcTdfromq(q2_withdamp(1,:,besti,bestj).*1000),NaN.*ones(1,48))';plot(thists,'color',colors('red'),'linewidth',1.5);
        thists=squeeze(calcTdfromq(q_era5(besti_era5,bestj_era5,1825:1825+167)));plot(thists,'color',colors('medium green'),'linewidth',1.5);

    subplot(3,2,5);plot(reshape(squeeze(subdailywinddir_all(7,firstday:lastday,:))',[numhwdays*24 1]),'k','linewidth',1.5);maketitle('Wind Direction');hold on;
        xlim([1 numhwdays*24]);set(gca,'xtick',1:24:numhwdays*24,'xticklabel',firstday:lastday);set(gca,'fontweight','bold','fontname','arial','fontsize',11);
        xlabel('Day of Year in 2020','fontweight','bold','fontname','arial','fontsize',12);ylim([0 360]);
        thists=cat(2,NaN.*ones(1,96),winddirfromuandv(u10_nourban(1,:,besti,bestj),v10_nourban(1,:,besti,bestj)),NaN.*ones(1,48))';plot(thists,'color',colors('orange'),'linewidth',1.5);
        thists=cat(2,NaN.*ones(1,96),winddirfromuandv(u10_withdamp(1,:,besti,bestj),v10_withdamp(1,:,besti,bestj)),NaN.*ones(1,48))';plot(thists,'color',colors('red'),'linewidth',1.5);
        thists=squeeze(winddirfromuandv(squeeze(u_era5(besti_era5,bestj_era5,1825:1825+167)),squeeze(v_era5(besti_era5,bestj_era5,1825:1825+167))));plot(thists,'color',colors('medium green'),'linewidth',1.5);

    subplot(3,2,6);plot(reshape(squeeze(subdailywindspd_all(7,firstday:lastday,:))',[numhwdays*24 1]),'k','linewidth',1.5);maketitle('Wind Speed');hold on;
        xlim([1 numhwdays*24]);set(gca,'xtick',1:24:numhwdays*24,'xticklabel',firstday:lastday);set(gca,'fontweight','bold','fontname','arial','fontsize',11);
        xlabel('Day of Year in 2020','fontweight','bold','fontname','arial','fontsize',12);
        thists=cat(2,NaN.*ones(1,96),sqrt(u10_nourban(1,:,besti,bestj).^2+v10_nourban(1,:,besti,bestj).^2),NaN.*ones(1,48))';plot(thists,'color',colors('orange'),'linewidth',1.5);
        thists=cat(2,NaN.*ones(1,96),sqrt(u10_withdamp(1,:,besti,bestj).^2+v10_withdamp(1,:,besti,bestj).^2),NaN.*ones(1,48))';plot(thists,'color',colors('red'),'linewidth',1.5);
        thists=squeeze(sqrt(u_era5(besti_era5,bestj_era5,1825:1825+167).^2+v_era5(besti_era5,bestj_era5,1825:1825+167).^2));plot(thists,'color',colors('medium green'),'linewidth',1.5);

    set(gcf,'color','w');
    curpart=1;highqualityfiguresetup;figname='abudhabitimeseries';curpart=2;highqualityfiguresetup;
    

    %tw1d=reshape(squeeze(subdailytw_all(7,firstday:lastday,:))',[numeventdays*24 1]);
    %twtendency=diff(tw1d);
    %winddir1d=reshape(squeeze(subdailywinddir_all(7,firstday:lastday,:))',[numeventdays*24 1]);
    %figure(86);clf;
    %scatter(winddir1d,tw1d);
end





if getvariousstats==1
    %p95 and all-time max at each point or station
    %this is computed directly from *hourly* data

    %ERA5
    tw_era5_alltimemax=NaN.*ones(nlat,nlon);tw_era5_p95=NaN.*ones(nlat,nlon);
    for i=1:nlat
        for j=1:nlon
            tw_era5_alltimemax(i,j)=max(tw_era5(i,j,:));
            tw_era5_p95(i,j)=quantile(tw_era5(i,j,:),0.95);
            maxversusp95_era5(i,j)=tw_era5_alltimemax(i,j)-tw_era5_p95(i,j);
        end
    end

    %ERA5-Land
    %tw_era5l_alltimemax=NaN.*ones(e5llatsz,e5llonsz);tw_era5l_p95=NaN.*ones(e5llatsz,e5llonsz);
    %for i=1:e5llatsz
    %    for j=1:e5llonsz
    %        tw_era5l_alltimemax(i,j)=max(tw_era5land(i,j,:));
    %        tw_era5l_p95(i,j)=quantile(tw_era5land(i,j,:),0.95);
    %        maxversusp95_era5l(i,j)=tw_era5l_alltimemax(i,j)-tw_era5l_p95(i,j);
    %    end
    %end

    %Weather stations
    tw_inlandstns_p95=quantile(reshape(subdailytw_inland,[size(subdailytw_inland,1) size(subdailytw_inland,2)*size(subdailytw_inland,3)]),0.95,2);
    tw_sPGstns_p95=quantile(reshape(subdailytw_sPG,[size(subdailytw_sPG,1) size(subdailytw_sPG,2)*size(subdailytw_sPG,3)]),0.95,2);
    tw_wPGstns_p95=quantile(reshape(subdailytw_wPG,[size(subdailytw_wPG,1) size(subdailytw_wPG,2)*size(subdailytw_wPG,3)]),0.95,2);

    %PWSs -- exactly 12 of 24 meet min-sample-size-of-100 criterion
    tw_pws_p95=NaN.*ones(numpws,1);pwsgoodstns=zeros(numpws,1);
    for s=1:numpws
        samplesize=squeeze(sum(~isnan(subdailytw_pws_hrs(s,:,12))));
        if samplesize>=100 %minimum sample size
            tw_pws_p95(s)=quantile(reshape(squeeze(subdailytw_pws_hrs(s,:,:)),[size(subdailytw_pws_hrs,2)*24 1]),0.95);
            pwsgoodstns(s)=1;
        end
    end
    

    %ERA5-Land interpolated to ERA5 grid
    %tw_era5l_p95_era5res=interp2(lon_era5land,lat_era5land,tw_era5l_p95,lon_era5,lat_era5);

    %Compare!
    %e5l=tw_era5l_p95_era5res(57:105,29:109);
    %e5=tw_era5_p95(57:105,29:109);
    %imsq(e5l-e5);
    %ERA5-Land is significantly (1-2C) colder than ERA5, for some reason...
end


if calcspreads_pointdata==1
    %First, for weather stations -- again, time zone is UTC
    for reg=1:size(regnames,1)
        r=regnames{reg};
        twp95vec=eval(['tw_' r 'stns_p95;']);subdailytw=eval(['subdailytw_' r ';']);subdailyt=eval(['subdailyt_' r ';']);

        %Wind directions
        subdailywinddir=eval(['subdailywinddir_' r ';']);
        centralwinddirs=NaN.*ones(length(twp95vec),1);
        for stn=1:length(twp95vec)
            thisp95(stn)=twp95vec(stn);c=0;clear winddirvec;
            for hr=1:24
                for day=1:size(subdailytw,2)
                    if subdailytw(stn,day,hr)>thisp95(stn)
                        c=c+1;
                        winddirvec(c)=subdailywinddir(stn,day,hr);
                    end
                end
            end
            peakdir=mode(winddirvec);
            frac_within45deg=round2(sum(abs(winddirvec-peakdir)<=45)/c,0.05);
            if frac_within45deg>=minfrac
                centralwinddirs(stn)=peakdir;
            end
        end
        eval(['centralwinddirs_' r '=centralwinddirs;']);

        %Hours of day
        clear hrofmaxp95prob_tw;clear centralhours_local;centralhours=NaN.*ones(length(twp95vec),1);
        for stn=1:length(twp95vec)
            for hr=1:24
                datathishr=squeeze(subdailytw(stn,:,hr));numnonnan_hr=sum(~isnan(datathishr));
                probp95_tw(stn,hr)=sum(datathishr>thisp95(stn))/(numnonnan_hr);
            end
            [~,hrofmaxp95prob_tw(stn)]=max(probp95_tw(stn,:));
        
            %Do at least 50% of p95 values occur within 3 hours of the peak probability, i.e. within a 6-hour window?
            hrmax=hrofmaxp95prob_tw(stn);
            if hrmax>=22
                tosum=[probp95_tw(stn,hrmax-3:24) probp95_tw(stn,1:hrmax-21)];
            elseif hrmax<=3
                tosum=[probp95_tw(stn,hrmax+21) probp95_tw(stn,1:hrmax+3)];
            else
                tosum=[probp95_tw(stn,hrmax-3:hrmax+3)];
            end
            sum_within3=sum(tosum);
            sum_all=sum(probp95_tw(stn,:));
            frac_within3=round2(sum_within3/sum_all,0.05);
            if frac_within3>=minfrac
                centralhours(stn)=hrmax;
            end
            %Convert to LST
            chlocal_tmp=centralhours(stn)+tzadj;
            chlocal_tmp(chlocal_tmp>=24)=chlocal_tmp(chlocal_tmp>=24)-24;
            centralhours_local(stn)=chlocal_tmp;
        end
        eval(['centralhours_local_' r '=centralhours_local;']);
    end


    %Second, for personal weather stations (PWSs)
    twp95vec=tw_pws_p95;subdailytw=subdailytw_pws_hrs;subdailyt=subdailyt_pws_hrs;
    
    %Wind directions
    subdailywinddir=subdailywinddir_pws_hrs;
    centralwinddirs_pws=NaN.*ones(numpws,1);
    for stn=1:length(twp95vec)
        if pwsgoodstns(stn)==1
            thisp95(stn)=twp95vec(stn);c=0;clear winddirvec;
            for hr=1:24
                for day=1:size(subdailytw,2)
                    if subdailytw(stn,day,hr)>thisp95(stn)
                        c=c+1;
                        winddirvec(c)=subdailywinddir(stn,day,hr);
                    end
                end
            end
            exist winddirvec;
            if ans==1
                peakdir=mode(winddirvec);
                frac_within45deg=round2(sum(abs(winddirvec-peakdir)<=45)/c,0.05);
                if frac_within45deg>=minfrac
                    centralwinddirs_pws(stn)=peakdir;
                end
            end
        end
    end

    %Hours of day
    hrofmaxp95prob_tw=NaN.*ones(numpws,1);centralhours_local_pws=NaN.*ones(numpws,1);centralhours=NaN.*ones(length(twp95vec),1);
    for stn=1:length(twp95vec)
        if pwsgoodstns(stn)==1
            for hr=1:24
                datathishr=squeeze(subdailytw(stn,:,hr));numnonnan_hr=sum(~isnan(datathishr));
                probp95_tw(stn,hr)=sum(datathishr>thisp95(stn))/(numnonnan_hr);
            end
            [~,hrofmaxp95prob_tw(stn)]=max(probp95_tw(stn,:));
        
            %Do at least 50% of p95 values occur within 3 hours of the peak?
            hrmax=hrofmaxp95prob_tw(stn);
            if hrmax>=22
                tosum=[probp95_tw(stn,hrmax-3:24) probp95_tw(stn,1:hrmax-21)];
            elseif hrmax<=3
                tosum=[probp95_tw(stn,hrmax+21) probp95_tw(stn,1:hrmax+3)];
            else
                tosum=[probp95_tw(stn,hrmax-3:hrmax+3)];
            end
            sum_within3=sum(tosum);
            sum_all=sum(probp95_tw(stn,:));
            frac_within3=round2(sum_within3/sum_all,0.05);
            if frac_within3>=minfrac
                centralhours(stn)=hrmax;
            end
            chlocal_tmp=centralhours(stn)+tzadj;chlocal_tmp(chlocal_tmp>=24)=chlocal_tmp(chlocal_tmp>=24)-24;
            centralhours_local_pws(stn)=chlocal_tmp;
        end
    end
end


if makefig1_revised==1
    figure(800);clf;

    %Panel a: map
    ax=subplot(10,10,100);gca=ax;
    data={lat_era5_10xres;lon_era5_10xres;myregs_era5_10xres};
    cmap=newregcolors;
    vararginnew={'mapproj';'mercator';'datatounderlay';data;'underlaycaxismin';0.5;'underlaycaxismax';3.5;'underlaystepsize';1;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'conttoplot';'all';'nonewfig';1;'omitfirstsubplotcolorbar';1;'underlaytransparency';0.5};
    datatype='custom';region={wb;nb;eb;sb};plotModelData(data,region,vararginnew,datatype);

    t=text(0.19,0.21,'UAE','units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',12);
    t=text(0.35,0.15,'Oman','units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',12);
    t=text(0.046,0.258,'Qatar','units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',12,'rotation',90);
    t=text(0.08,0.1,{'Saudi','Arabia'},'units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',12);
    t=text(0.32,0.43,'Iran','units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',12);
    t=text(0.1,0.318,{'Persian/','Arabian','Gulf'},'units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',12,'fontangle','italic');
    t=text(0.35,0.268,{'Gulf','of Oman'},'units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',12,'fontangle','italic');

    set(gca,'position',[0.16 0.52 0.47 0.45]);
    t=text(-0.05,0.51,splabels{1},'units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',12);
    %Lat/lon labels
    t=text(-0.045,0.28,'25N','units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',10);
    t=text(-0.045,0.005,'20N','units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',10);
    t=text(-0.018,0.515,'50E','units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',10);
    t=text(0.231,0.515,'55E','units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',10);
    t=text(0.48,0.515,'60E','units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',10);


    %Panel b: Tw values
    ax=subplot(10,10,100);gca=ax;
    cmin=25;cmax=31;
    if strcmp(datasettodo,'era5')
        data={lat_era5;lon_era5;tw_era5_p95};
    elseif strcmp(datasettodo,'era5land')
        data={lat_era5land;lon_era5land;tw_era5l_p95};
    end
    cmap=colormaps('wbt','more','not');
    vararginnew={'mapproj';'mercator';'datatounderlay';data;'underlaycaxismin';cmin;'underlaycaxismax';cmax;'underlaystepsize';0.5;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'conttoplot';'all';'nonewfig';1;'omitfirstsubplotcolorbar';0;...
        'colorbarfontsize';12};
    datatype='custom';region={wb;nb;eb;sb};plotModelData(data,region,vararginnew,datatype);
    hold on;

    %Overlay PWSs
    for stn=1:numpws
        thisp95val=tw_pws_p95(stn);
        if ~isnan(thisp95val)
            if thisp95val<cmin
                valcolor=cmap(1,:);
            elseif thisp95val>cmax
                valcolor=cmap(end,:);
            else
                valcolor=cmap(round(size(cmap,1)*(thisp95val-cmin)/(cmax-cmin)),:);
            end
            geoshow(pwsstnlats(stn),pwsstnlons(stn),'DisplayType','Point','Marker','d','MarkerFaceColor',valcolor,'MarkerEdgeColor','k','MarkerSize',8,'linewidth',1.5);
        end
    end

    %Overlay stn data
    %Put it on top because it's more reliable
    for reg=1:3
        r=regnames{reg};stninfo=eval(['stninfo_' r]);
        for stn=1:size(stninfo,1)
            thisp95val=eval(['tw_' r 'stns_p95(stn);']);
            if ~isnan(thisp95val)
                if thisp95val<cmin
                    valcolor=cmap(1,:);
                elseif thisp95val>cmax
                    valcolor=cmap(end,:);
                else
                    valcolor=cmap(round(size(cmap,1)*(thisp95val-cmin)/(cmax-cmin)),:);
                end
                geoshow(stninfo(stn,1),stninfo(stn,2),'DisplayType','Point','Marker','o','MarkerFaceColor',valcolor,'MarkerEdgeColor','k','MarkerSize',8,'linewidth',3);
            end
        end
    end
    maketitle('Tw p95');
    set(gca,'position',[0.52 0.52 0.47 0.45]);
    t=text(-0.05,0.51,splabels{2},'units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',12);

    
    %Panel c: central hours
    ax=subplot(10,10,100);gca=ax;
    data={lat_era5;lon_era5;centralhours_local_era5};
    cmap=colormaps('mygbm','24','not');cmin=1;cmax=24;
    vararginnew={'mapproj';'mercator';'datatounderlay';data;'underlaycaxismin';cmin;'underlaycaxismax';cmax;'underlaystepsize';1;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'variable';'wind';'conttoplot';'Asia';'nonewfig';1;'omitfirstsubplotcolorbar';0;...
        'colorbarfontsize';12;'colorbarticks';[4:4:24]};
    datatype='custom';region={wb;nb;eb;sb};plotModelData(data,region,vararginnew,datatype);

    %Overlay PWSs (diamonds)
    for stn=1:numpws
        stnhour=centralhours_local_pws(stn);
        if ~isnan(stnhour)
            if stnhour<cmin
                valcolor=cmap(1,:);
            elseif stnhour>cmax
                valcolor=cmap(end,:);
            else
                valcolor=cmap(round(size(cmap,1)*(stnhour-cmin)/(cmax-cmin)),:);
            end
            geoshow(pwsstnlats(stn),pwsstnlons(stn),'DisplayType','Point','Marker','d','MarkerFaceColor',valcolor,'MarkerEdgeColor','k','MarkerSize',8,'linewidth',1.5);
        end
    end

    %Overlay stn data (circles)
    hold on;
    for reg=1:3
        r=regnames{reg};stninfo=eval(['stninfo_' r]);
        for stn=1:size(stninfo,1)
            stnhour=eval(['centralhours_local_' r '(stn);']);
            if ~isnan(stnhour)
                if stnhour<=cmin
                    valcolor=cmap(1,:);
                elseif stnhour>=cmax
                    valcolor=cmap(end,:);
                else
                    valcolor=cmap(round(size(cmap,1)*(stnhour-cmin)/(cmax-cmin)),:);
                end
    
                geoshow(stninfo(stn,1),stninfo(stn,2),'DisplayType','Point','Marker','o','MarkerFaceColor',valcolor,'MarkerEdgeColor','k','MarkerSize',8,'linewidth',3);
            end
        end
    end
    maketitle('Hour of Peak Tw (LST)');
    set(gca,'position',[0.16 0.02 0.47 0.45]);
    t=text(-0.05,0.51,splabels{3},'units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',12);


    %Panel d: central wind directions
    ax=subplot(10,10,100);gca=ax;
    data={lat_era5;lon_era5;centralwinddirs_era5};
    cmap=colormaps('mygbm','12','not');cmin=0;cmax=360;
    vararginnew={'mapproj';'mercator';'datatounderlay';data;'underlaycaxismin';cmin;'underlaycaxismax';cmax;'underlaystepsize';30;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'variable';'wind';'conttoplot';'Asia';'nonewfig';1;'omitfirstsubplotcolorbar';0;...
        'colorbarfontsize';12;'colorbarticks';[0:60:360]};
    datatype='custom';region={wb;nb;eb;sb};plotModelData(data,region,vararginnew,datatype);

    %Overlay PWSs
    for stn=1:numpws
        stnwinddir=centralwinddirs_pws(stn);
        if ~isnan(stnwinddir)
            if stnwinddir<=cmin
                valcolor=cmap(1,:);
            elseif stnwinddir>=cmax
                valcolor=cmap(end,:);
            else
                valcolor=cmap(round(size(cmap,1)*(stnwinddir-cmin)/(cmax-cmin)),:);
            end
            geoshow(pwsstnlats(stn),pwsstnlons(stn),'DisplayType','Point','Marker','d','MarkerFaceColor',valcolor,'MarkerEdgeColor','k','MarkerSize',8,'linewidth',1.5);
        end
    end

    %Overlay stn data
    hold on;
    for reg=1:3
        r=regnames{reg};stninfo=eval(['stninfo_' r]);
        for stn=1:size(stninfo,1)
            stnwinddir=eval(['centralwinddirs_' r '(stn);']);
            if ~isnan(stnwinddir)
                if stnwinddir<=cmin
                    valcolor=cmap(1,:);
                elseif stnwinddir>=cmax
                    valcolor=cmap(end,:);
                else
                    valcolor=cmap(round(size(cmap,1)*(stnwinddir-cmin)/(cmax-cmin)),:);
                end
    
                geoshow(stninfo(stn,1),stninfo(stn,2),'DisplayType','Point','Marker','o','MarkerFaceColor',valcolor,'MarkerEdgeColor','k','MarkerSize',8,'linewidth',3);
            end
        end
    end
    maketitle('Wind Direction of Peak Tw');
    set(gca,'position',[0.52 0.02 0.47 0.45]);
    t=text(-0.05,0.51,splabels{4},'units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',12);


    if minfrac==0.5;figname='fig1_rev';elseif minfrac==0.4;figname='figsX_minfrac04';else;disp('See line 2600!');return;end
    set(gcf,'color','w');curpart=1;highqualityfiguresetup;curpart=2;highqualityfiguresetup;

    clear gca;
end


%Shaded maps of T and q for T- and Tw-defined hot days, at 4-hour intervals
if calcs_diurnalsequence==1
    nlat=size(tw_era5,1);nlon=size(tw_era5,2);nhrs=size(tw_era5,3);
    exist t2m_days;if ans==0;exist t_era5;if ans==0;tmp=load(strcat(saveloc,'tera5.mat'));t_era5=tmp.t_era5;clear tmp;end;...
            t2m_days=reshape(t_era5,[nlat nlon 24 climodaylen*nyr]);clear t_era5;end
    exist q2m_days;if ans==0;exist q_era5;if ans==0;tmp=load(strcat(saveloc,'qera5.mat'));q_era5=tmp.q_era5;clear tmp;end;...
            q2m_days=reshape(q_era5,[nlat nlon 24 climodaylen*nyr]);clear q_era5;end

    dothis=0;
    if dothis==1
    exist q75whereabove;
    if ans==0
        q75whereabove=NaN.*ones(era5lonsz,era5latsz,24,climodaylen*nyr);
        t75whereabove=NaN.*ones(era5lonsz,era5latsz,24,climodaylen*nyr);
        for i=1:era5lonsz
            for j=1:era5latsz
                %95th pctile across all hours at this gridcell
                q75_gridcell=quantile(reshape(squeeze(q2m_days(i,j,:,:)),[24*climodaylen*nyr 1]),0.75);
                %for hr=1:24;q75whereabove(i,j,hr,:)=squeeze(q2m_days(i,j,hr,:))>=q75_thisgridcell;end 
    
                t75_gridcell=quantile(reshape(squeeze(t2m_days(i,j,:,:)),[24*climodaylen*nyr 1]),0.75);
                %for hr=1:24;t75whereabove(i,j,hr,:)=squeeze(t2m_days(i,j,hr,:))>=t75_thisgridcell;end
            end
        end
    end
    end

    %Identify p95 hot days by T and by Tw
    %These are days when any hour is above the overall p95 at that gridcell
    %Runtime: 6 min
    exist tw95whereabove;
    if ans==0
        tw95_gridcell=NaN.*ones(nlon,nlat);t95_gridcell=NaN.*ones(nlon,nlat);
        q75_gridcell=NaN.*ones(nlon,nlat);q90_gridcell=NaN.*ones(nlon,nlat);
        t75_gridcell=NaN.*ones(nlon,nlat);t90_gridcell=NaN.*ones(nlon,nlat);
        for i=1:nlon
            for j=1:nlat
                tw95_gridcell(i,j)=quantile(reshape(squeeze(tw2m_days(i,j,:,:)),[24*climodaylen*nyr 1]),0.95);
                t95_gridcell(i,j)=quantile(reshape(squeeze(t2m_days(i,j,:,:)),[24*climodaylen*nyr 1]),0.95);

                q75_gridcell(i,j)=quantile(reshape(squeeze(q2m_days(i,j,:,:)),[24*climodaylen*nyr 1]),0.75);
                q90_gridcell(i,j)=quantile(reshape(squeeze(q2m_days(i,j,:,:)),[24*climodaylen*nyr 1]),0.90);

                t75_gridcell(i,j)=quantile(reshape(squeeze(t2m_days(i,j,:,:)),[24*climodaylen*nyr 1]),0.75);
                t90_gridcell(i,j)=quantile(reshape(squeeze(t2m_days(i,j,:,:)),[24*climodaylen*nyr 1]),0.90);
            end
        end

        twp95rep=repmat(tw95_gridcell,[1 1 24 climodaylen*nyr]);
        tw95whereabove=tw2m_days>=twp95rep;clear twp95rep;

        tp95rep=repmat(t95_gridcell,[1 1 24 climodaylen*nyr]);
        t95whereabove=t2m_days>=tp95rep;clear tp95rep;
    end
    save(strcat(saveloc,'hotdaysstarts.mat'),'tw95_gridcell','t95_gridcell','q75_gridcell','q90_gridcell','t75_gridcell','t90_gridcell',...
        'tw95whereabove','t95whereabove');


    %Define hot days (25 min)
    %Must be redone for each new place hot days are based upon
    c_tw=0;clear tarray_twhotdays;clear qarray_twhotdays;clear uarray_twhotdays;clear varray_twhotdays;
    c_t=0;clear tarray_thotdays;clear qarray_thotdays;clear uarray_thotdays;clear varray_thotdays;

    exist u_era5;if ans==0;tmp=load(strcat(saveloc,'uera5.mat'));u_era5=tmp.u_era5;clear tmp;u10m_days=reshape(u_era5,[nlat nlon 24 climodaylen*nyr]);clear u_era5;end
    exist v_era5;if ans==0;tmp=load(strcat(saveloc,'vera5.mat'));v_era5=tmp.v_era5;clear tmp;v10m_days=reshape(v_era5,[nlat nlon 24 climodaylen*nyr]);clear v_era5;end
    for day=1:climodaylen*nyr
        %Is this a Tw-defined hot day?
        if sum(tw95whereabove(myi,myj,:,day))>=1 %if so, get data for whole region on this day
            c_tw=c_tw+1;
            tarray_twhotdays(c_tw,:,:,:)=t2m_days(:,:,:,day); 
            qarray_twhotdays(c_tw,:,:,:)=q2m_days(:,:,:,day);
            uarray_twhotdays(c_tw,:,:,:)=u10m_days(:,:,:,day);
            varray_twhotdays(c_tw,:,:,:)=v10m_days(:,:,:,day);
        end

        %Separately: is this a T-defined hot day?
        if sum(t95whereabove(myi,myj,:,day))>=1 %if so, get data for whole region on this day
            c_t=c_t+1;
            tarray_thotdays(c_t,:,:,:)=t2m_days(:,:,:,day); 
            qarray_thotdays(c_t,:,:,:)=q2m_days(:,:,:,day);
            uarray_thotdays(c_t,:,:,:)=u10m_days(:,:,:,day);
            varray_thotdays(c_t,:,:,:)=v10m_days(:,:,:,day);
        end
    end

    tarray_twhotdays_mean=squeeze(mean(tarray_twhotdays));
    qarray_twhotdays_mean=squeeze(mean(qarray_twhotdays));
    uarray_twhotdays_mean=squeeze(mean(uarray_twhotdays));
    varray_twhotdays_mean=squeeze(mean(varray_twhotdays));

    tarray_thotdays_mean=squeeze(mean(tarray_thotdays));
    qarray_thotdays_mean=squeeze(mean(qarray_thotdays));
    uarray_thotdays_mean=squeeze(mean(uarray_thotdays));
    varray_thotdays_mean=squeeze(mean(varray_thotdays));

    save(strcat(saveloc,'lcluchotdaysoutput_',hotdaysbasedon),'tarray_twhotdays_mean','qarray_twhotdays_mean','uarray_twhotdays_mean','varray_twhotdays_mean',...
        'tarray_thotdays_mean','qarray_thotdays_mean','uarray_thotdays_mean','varray_thotdays_mean');

    clear t95whereabove;clear tw95whereabove;
end

if finaldiurnalseqprep==1
    %All based on whatever hotdaysbasedon location is currently selected!
    %Get means of arrays on hot days (10 sec)
    %Convert T and q arrays to binary versions -- one against the 75th pct, the other against the 90th
    tarray_twhotdays_abovep95=tarray_twhotdays_mean>=repmat(t95_gridcell,[1 1 24]);
    tarray_twhotdays_abovep90=tarray_twhotdays_mean>=repmat(t90_gridcell,[1 1 24]);tarray_twhotdays_abovep75=tarray_twhotdays_mean>=repmat(t75_gridcell,[1 1 24]);
    qarray_twhotdays_abovep90=qarray_twhotdays_mean>=repmat(q90_gridcell,[1 1 24]);qarray_twhotdays_abovep75=qarray_twhotdays_mean>=repmat(q75_gridcell,[1 1 24]);

    tarray_thotdays_abovep95=tarray_thotdays_mean>=repmat(t95_gridcell,[1 1 24]);
    tarray_thotdays_abovep90=tarray_thotdays_mean>=repmat(t90_gridcell,[1 1 24]);tarray_thotdays_abovep75=tarray_thotdays_mean>=repmat(t75_gridcell,[1 1 24]);
    qarray_thotdays_abovep90=qarray_thotdays_mean>=repmat(q90_gridcell,[1 1 24]);qarray_thotdays_abovep75=qarray_thotdays_mean>=repmat(q75_gridcell,[1 1 24]);

    %Define in terms of hours (20 sec)
    %Localrelhr: 12am, 4am, 8am, 12pm, 4pm, 8pm LST
    for hr=1:4:21 %UTC, so LST is +4 from this; also note hr 1 = 00 UTC
        relhr=round2(hr/4,1,'ceil');
        localrelhr=relhr+1;
        if localrelhr==7;localrelhr=1;end %hr 21 = 20 UTC = 00 LST

        tmap95_twhotdays(:,:,localrelhr)=squeeze(tarray_twhotdays_abovep95(:,:,hr));
        tmap90_twhotdays(:,:,localrelhr)=squeeze(tarray_twhotdays_abovep90(:,:,hr));
        tmap75_twhotdays(:,:,localrelhr)=squeeze(tarray_twhotdays_abovep75(:,:,hr));
        qmap90_twhotdays(:,:,localrelhr)=squeeze(qarray_twhotdays_abovep90(:,:,hr));
        qmap75_twhotdays(:,:,localrelhr)=squeeze(qarray_twhotdays_abovep75(:,:,hr));

        tmap95_thotdays(:,:,localrelhr)=squeeze(tarray_thotdays_abovep95(:,:,hr));
        tmap90_thotdays(:,:,localrelhr)=squeeze(tarray_thotdays_abovep90(:,:,hr));
        tmap75_thotdays(:,:,localrelhr)=squeeze(tarray_thotdays_abovep75(:,:,hr));
        qmap90_thotdays(:,:,localrelhr)=squeeze(qarray_thotdays_abovep90(:,:,hr));
        qmap75_thotdays(:,:,localrelhr)=squeeze(qarray_thotdays_abovep75(:,:,hr));


        umap_twhotdays(:,:,localrelhr)=squeeze(uarray_twhotdays_mean(:,:,hr));
        vmap_twhotdays(:,:,localrelhr)=squeeze(varray_twhotdays_mean(:,:,hr));
        umap_thotdays(:,:,localrelhr)=squeeze(uarray_thotdays_mean(:,:,hr));
        vmap_thotdays(:,:,localrelhr)=squeeze(varray_thotdays_mean(:,:,hr));
    end
end


if map_diurnalsequence==1
    %Make plots

    %Top row: T-defined hot days
    %Second row: Tw-defined hot days

    %Columns: 12am, 4am, 8am, 12pm, 4pm, 8pm LST

    %Color scheme:
    %greens: q>75 or >90
    %oranges: T>75 or >90
    %purples: T&q>75 or >90
    %if T and q are of a different category, use dominant one
    %e.g. T>90 but q only >75 --> shade as dark orange
    figure(11);clf;hold on;clear gca;
    figname=strcat('diurnalsequencemap_',hotdaysbasedon);
    lefts=[0.03 0.19 0.35 0.51 0.67 0.83;0.03 0.19 0.35 0.51 0.67 0.83];
    bottoms=[0.5.*ones(1,6);0.2.*ones(1,6)];
    splabels_row1={'a)';'b)';'c)';'d)';'e)';'f)'};splabels_row2={'g)';'h)';'i)';'j)';'k)';'l)'};
    timelabels={'12am';'4am';'8am';'12pm';'4pm';'8pm'};
    
    for row=1:2
        if row==1
            qarr75=qmap75_thotdays;qarr90=qmap90_thotdays;tarr75=tmap75_thotdays;tarr90=tmap90_thotdays;uarr=umap_thotdays;varr=vmap_thotdays;splabels=splabels_row1;
        elseif row==2
            qarr75=qmap75_twhotdays;qarr90=qmap90_twhotdays;tarr75=tmap75_twhotdays;tarr90=tmap90_twhotdays;uarr=umap_twhotdays;varr=vmap_twhotdays;splabels=splabels_row2;
        end
        for col=1:6 %times of day
    
            subplot(10,10,100);hold on;
    
            mymap_u=double(squeeze(uarr(:,:,col)));mymap_v=double(squeeze(varr(:,:,col)));
            q75=double(squeeze(qarr75(:,:,col)));q90=double(squeeze(qarr90(:,:,col)));
            t75=double(squeeze(tarr75(:,:,col)));t90=double(squeeze(tarr90(:,:,col)));
            mymap_q=zeros(era5lonsz,era5latsz);mymap_t=zeros(era5lonsz,era5latsz);mymap=zeros(era5lonsz,era5latsz);mymap_both=zeros(era5lonsz,era5latsz);
            mymap_q_holes=zeros(era5lonsz,era5latsz);mymap_t_holes=zeros(era5lonsz,era5latsz);
            for i=1:era5lonsz
                for j=1:era5latsz
                    if q90(i,j)==1;mymap_q(i,j)=5;elseif q75(i,j)==1;mymap_q(i,j)=4;end
                    if t90(i,j)==1;mymap_t(i,j)=2;elseif t75(i,j)==1;mymap_t(i,j)=1;end

                    mymap_q_holes(i,j)=mymap_q(i,j);mymap_t_holes(i,j)=mymap_t(i,j); %preliminarily...

                    if q90(i,j)==1 && t90(i,j)==1 %high tie
                        mymap(i,j)=8;mymap_both(i,j)=8;mymap_q_holes(i,j)=0;mymap_t_holes(i,j)=0;
                    elseif q90(i,j)==1 && t90(i,j)==0 %q wins
                        mymap(i,j)=5;mymap_t_holes(i,j)=0;
                    elseif t90(i,j)==1 && q90(i,j)==0 %t wins
                        mymap(i,j)=2;mymap_q_holes(i,j)=0;
                    elseif q75(i,j)==1 && t75(i,j)==1 %low tie
                        mymap(i,j)=7;mymap_both(i,j)=7;mymap_q_holes(i,j)=0;mymap_t_holes(i,j)=0;
                    elseif q75(i,j)==1 && t75(i,j)==0 %q wins
                        mymap(i,j)=4;mymap_t_holes(i,j)=0;
                    elseif t75(i,j)==1 && q75(i,j)==0 %t wins
                        mymap(i,j)=1;mymap_q_holes(i,j)=0;
                    end
                end
            end
    
            myax=worldmap([sb nb],[wb eb]);tightmap;set(gcf,'color','w');
            set(findall(myax,'Tag','PLabel'),'visible','off');set(findall(myax,'Tag','MLabel'),'visible','off');

            R=georefcells([min(min(lat_era5)) max(max(lat_era5))],[min(min(lon_era5)) max(max(lon_era5))],size(mymap_t),'ColumnsStartFrom','north');
            
            maxval=8;v=-0.5:maxval-0.5;
            cmap=[colormaps('whitelightorangedarkorange','more','not');colormaps('whitelightgreendarkgreen','more','not');...
                colormaps('whitelightpurpledarkpurple','more','not')];
            colormap(cmap);clim([0 maxval+0.05]);

            Z1=mymap_t_holes;
            settonan=Z1==0;Z1(settonan)=NaN;Z1(1,1)=0; %set blank areas to NaN
            contourfm(Z1,R,v,'LineColor','none');alpha(0.7);

            if row==1 && col==1
                xticks=[4/3-0.1:4/3:8];xticklabel={'T>p75';'T>p90';'q>p75';'q>p90';'T,q>p75';'T,q>p90'};
                cb=colorbar('horiz');set(cb,'position',[0.3 0.14 0.4 0.02],'xtick',xticks,'xticklabel',xticklabel,'fontweight','bold','fontsize',10);
            end
                        
            Z2=mymap_q_holes;
            settonan=Z2==0;Z2(settonan)=NaN;Z2(1,1)=0; %set blank areas to NaN to avoid overwriting previously plotted colors with white
            contourfm(Z2,R,v,'LineColor','none');alpha(0.7);

            Z3=mymap_both;
            settonan=Z3==0;Z3(settonan)=NaN;Z3(1,1)=0; %set blank areas to NaN to avoid overwriting previously plotted colors with white
            v=6.5:8.5;
            contourfm(Z3,R,v,'LineColor','none');alpha(0.7);
            
            conttoplot='Asia';co='w';cbc='k';clw=1;addborders;
            set(findall(myax,'Tag','PLabel'),'visible','off');set(findall(myax,'Tag','MLabel'),'visible','off');
            
            set(gca,'position',[lefts(row,col) bottoms(row,col) 0.15 0.25]);
    
            quivermc(lat_era5,lon_era5,mymap_u,mymap_v,'dontaddtext','maparea',((nb-sb)*(eb-wb)),'elongationfactor',9,'skipfactor',4,'linewidth',1.3);

            t=text(-0.03,0.53,splabels{col},'units','normalized');set(t,'fontsize',11,'fontweight','bold','fontname','arial');
            title(timelabels{col},'fontsize',12,'fontweight','bold','fontname','arial');

            %Show location with red pentagram
            geoshow(mylat,mylon,'DisplayType','Point','Marker','p','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',9);

            %Legend
            if row==2 && col==6
                t=annotation('arrow',[0.815 0.84],[0.147 0.147],'color','k','linewidth',1.5,'headstyle','plain','headwidth',6,'headlength',7);
                t=annotation('arrow',[0.90 0.915],[0.147 0.147],'color','k','linewidth',1.7,'headstyle','none');
                t=text(-0.065,-0.15,'5 m/s','units','normalized');set(t,'fontsize',10.5,'fontweight','bold','fontname','arial');
                t=text(0.187,-0.15,'100 km','units','normalized');set(t,'fontsize',10.5,'fontweight','bold','fontname','arial');
            end

            %Lat/lon labels
            if col==1
                toshow=strcat('25',char(176));t=text(-0.065,0.28,toshow,'units','normalized');set(t,'fontsize',10.5,'fontweight','bold','fontname','arial');
                toshow=strcat('20',char(176));t=text(-0.07,0.005,toshow,'units','normalized');set(t,'fontsize',10.5,'fontweight','bold','fontname','arial');
            end
            if row==2 && rem(col,2)==1
                toshow=strcat('50',char(176));t=text(-0.025,-0.025,toshow,'units','normalized');set(t,'fontsize',10.5,'fontweight','bold','fontname','arial');
                toshow=strcat('60',char(176));t=text(0.475,-0.025,toshow,'units','normalized');set(t,'fontsize',10.5,'fontweight','bold','fontname','arial');
            end
        end
    end
    curpart=1;highqualityfiguresetup;curpart=2;highqualityfiguresetup;
end



if getvertprofiles==1
    if vertprofileprep==1
        if dorharmprep==1
            %RHARM
            %As in p90diurnalprob figure: combine data for available stns within 0.375 deg of center point 
            data00utcbyloc=cell(4,1);data12utcbyloc=cell(4,1);
            for loc=1:4
                c=0;
                centerlat=latsofint(loc);centerlon=lonsofint(loc);
                for stntest=1:size(tw_rharm_jja20002020_00utc,1)
                    if stnlats_rharm(stntest)>=centerlat-0.375 && stnlats_rharm(stntest)<=centerlat+0.375 &&...
                            stnlons_rharm(stntest)>=centerlon-0.375 && stnlons_rharm(stntest)<=centerlon+0.375
                        c=c+1;
                        data00utcbyloc{loc}(c,:,:,:,:)=tw_rharm_jja20002020_00utc(stntest,:,:,:,:);
                        data12utcbyloc{loc}(c,:,:,:,:)=tw_rharm_jja20002020_12utc(stntest,:,:,:,:);
                        %dims here are stns, years, months, days, then final dim is vert level -- 200, 300, 500, 700, 850, 925, 950, 975, 990, 1000
                    end
                end
            end
            %Turns out that only Abu Dhabi (loc 1) has anything... better than nothing!
    
            tmp=squeeze(data00utcbyloc{1}(1,:,:,:,:));tmp=permute(tmp,[2 1 3 4]);
            tmp2=reshape(tmp,[21*3 31 10]);tmp2=permute(tmp2,[2 1 3]); %now dims are months (in order) | days of month | vert levels
            tw_rharm00utc_loc1=reshape(tmp2,[63*31 10]); %now dims are days (in chronological order across all months and years) | vert levels
    
            tmp=squeeze(data12utcbyloc{1}(1,:,:,:,:));tmp=permute(tmp,[2 1 3 4]);
            tmp2=reshape(tmp,[21*3 31 10]);tmp2=permute(tmp2,[2 1 3]);
            tw_rharm12utc_loc1=reshape(tmp2,[63*31 10]);
    
             %Eliminate bad data (which for reasons I forgot and aren't that important varies by level: -6.77204, -7.16102, or several others...)
             bad=abs(tw_rharm00utc_loc1+6.42575)<10^-5;tw_rharm00utc_loc1(bad)=NaN;bad=abs(tw_rharm12utc_loc1+6.42575)<10^-5;tw_rharm12utc_loc1(bad)=NaN;
             bad=abs(tw_rharm00utc_loc1+6.77204)<10^-5;tw_rharm00utc_loc1(bad)=NaN;bad=abs(tw_rharm12utc_loc1+6.77204)<10^-5;tw_rharm12utc_loc1(bad)=NaN;
             bad=abs(tw_rharm00utc_loc1+7.16102)<10^-5;tw_rharm00utc_loc1(bad)=NaN;bad=abs(tw_rharm12utc_loc1+7.16102)<10^-5;tw_rharm12utc_loc1(bad)=NaN;
             bad=abs(tw_rharm00utc_loc1+8.10630)<10^-5;tw_rharm00utc_loc1(bad)=NaN;bad=abs(tw_rharm12utc_loc1+8.10630)<10^-5;tw_rharm12utc_loc1(bad)=NaN;
             bad=abs(tw_rharm00utc_loc1+9.91026)<10^-5;tw_rharm00utc_loc1(bad)=NaN;bad=abs(tw_rharm12utc_loc1+9.91026)<10^-5;tw_rharm12utc_loc1(bad)=NaN;
             bad=abs(tw_rharm00utc_loc1+13.0142)<10^-4;tw_rharm00utc_loc1(bad)=NaN;bad=abs(tw_rharm12utc_loc1+13.0142)<10^-4;tw_rharm12utc_loc1(bad)=NaN;
             bad=abs(tw_rharm00utc_loc1+15.743)<10^-4;tw_rharm00utc_loc1(bad)=NaN;bad=abs(tw_rharm12utc_loc1+15.743)<10^-4;tw_rharm12utc_loc1(bad)=NaN;
    
             %Eliminate non-existent day that was implicitly included (Jun 31 of each year)
             toremove=31:climodaylen:1953;
             newi=0;goodi=0;clear goodvec;
             for i=1:1953
                if ~ismember(i,toremove)
                    goodi=goodi+1;
                    goodvec(goodi)=i;
                end
             end
             tw_rharm00utc_loc1_new=tw_rharm00utc_loc1(goodvec,:);
             tw_rharm12utc_loc1_new=tw_rharm12utc_loc1(goodvec,:);
    
             clear myarr;
             %twprofile_alldays_rharm_byloc dims are loc | lev | hour
             myarr(1,:,2)=mean(tw_rharm00utc_loc1_new(:,3:10),'omitnan'); %500, 700, 850, 925, 950, 975, 990, 1000
             myarr(1,:,5)=mean(tw_rharm12utc_loc1_new(:,3:10),'omitnan');
             bad=myarr==0;myarr(bad)=NaN;
             twprofile_alldays_rharm_byloc=myarr;
    
             c=0;clear tmp00;clear tmp12;
             for day=1:climodaylen*nyr
                 if sum(tw95whereabove(besti_era5(1),bestj_era5(1),:,day))>=1
                     c=c+1;
                     tmp00(c,:)=tw_rharm00utc_loc1_new(day,:);
                     tmp12(c,:)=tw_rharm12utc_loc1_new(day,:);
                 end
             end
             clear twprofile_twhotdays_rharm_byloc;
             twprofile_twhotdays_rharm_byloc(1,:,2)=mean(tmp00(:,3:10),'omitnan');
             twprofile_twhotdays_rharm_byloc(1,:,5)=mean(tmp12(:,3:10),'omitnan');
             invalid=twprofile_twhotdays_rharm_byloc==0;twprofile_twhotdays_rharm_byloc(invalid)=NaN;

             save(strcat(saveloc,'hotdaysprofiles_rharm.mat'),'twprofile_alldays_rharm_byloc','twprofile_twhotdays_rharm_byloc');
        end
    end
end


if makefig5_revised==1
    if doera5prep==1
        %Reminder: levels obtained here were 500, 700, 850, 925, 975, 1000
        exist twlev_days;
        if ans==0
            f=load(strcat(saveloc,'twlevs.mat'));twlevs_era5=f.tw2m;clear f;
            twlev_days=reshape(twlevs_era5,[nlat nlon 6 6 climodaylen*nyr]);clear twlevs_era5; %dims are lat | lon | lev | hours in day | days
        end
        exist tlevs_days;
        if ans==0
            tlev_days=reshape(tlevs_era5,[nlat nlon 6 6 climodaylen*nyr]);clear tlevs_era5;
            qlev_days=reshape(qlevs_era5,[nlat nlon 6 6 climodaylen*nyr]);clear qlevs_era5;
        end

        %Define new regions and get vert profiles on all days and on p95 days
        %1 min
        c=zeros(3,1);twholder=cell(3,1);tholder=cell(3,1);qholder=cell(3,1);
        for i=1:nlat
            for j=1:nlon
                r=NaN;if myregs_era5(i,j)==1;r=1;elseif myregs_era5(i,j)==2;r=2;elseif myregs_era5(i,j)==3;r=3;end
                if ~isnan(r)
                    c(r)=c(r)+1;
                    twholder{r}(c(r),:,:,:)=twlev_days(i,j,:,:,:);
                    tholder{r}(c(r),:,:,:)=tlev_days(i,j,:,:,:);
                    qholder{r}(c(r),:,:,:)=qlev_days(i,j,:,:,:);
                end
            end
        end
        %Get region means and st devs
        twprofile_alldays_era5_mean=NaN.*ones(3,6,6);tprofile_alldays_era5_mean=NaN.*ones(3,6,6);qprofile_alldays_era5_mean=NaN.*ones(3,6,6);
        twprofile_alldays_era5_stdev=NaN.*ones(3,6,6);tprofile_alldays_era5_stdev=NaN.*ones(3,6,6);qprofile_alldays_era5_stdev=NaN.*ones(3,6,6);
        for r=1:3
            twprofile_alldays_era5_mean(r,:,:)=squeeze(mean(squeeze(mean(twholder{r},'omitnan')),3,'omitnan'));
            tprofile_alldays_era5_mean(r,:,:)=squeeze(mean(squeeze(mean(tholder{r},'omitnan')),3,'omitnan'));
            qprofile_alldays_era5_mean(r,:,:)=squeeze(mean(squeeze(mean(qholder{r},'omitnan')),3,'omitnan'));

            d1=size(twholder{r},1);d4=size(twholder{r},4);
            twprofile_alldays_era5_stdev(r,:,:)=squeeze(std(reshape(permute(twholder{r},[4 1 2 3]),[d1*d4 6 6])));
            tprofile_alldays_era5_stdev(r,:,:)=squeeze(std(reshape(permute(tholder{r},[4 1 2 3]),[d1*d4 6 6])));
            qprofile_alldays_era5_stdev(r,:,:)=squeeze(std(reshape(permute(qholder{r},[4 1 2 3]),[d1*d4 6 6])));
        end
    
        %Repeat for p95 days
        exist tw95whereabove;if ans==0;f=load(strcat(saveloc,'hotdaysstarts.mat'));tw95whereabove=f.tw95whereabove;end;clear tmp;
        twprofile_twhotdays=NaN.*ones(nlat,nlon,6,6);tprofile_twhotdays=NaN.*ones(nlat,nlon,6,6);qprofile_twhotdays=NaN.*ones(nlat,nlon,6,6);
        for i=1:nlat
            for j=1:nlon
                r=NaN;if myregs_era5(i,j)==1;r=1;elseif myregs_era5(i,j)==2;r=2;elseif myregs_era5(i,j)==3;r=3;end
                if ~isnan(r) %whether to bother with this point at all
                    c=0;tmptw=zeros(1,6,6);tmpt=zeros(1,6,6);tmpq=zeros(1,6,6);
                    for day=1:climodaylen*nyr
                        if sum(tw95whereabove(i,j,:,day))>=1
                            c=c+1;
                            tmptw(c,:,:)=twlev_days(i,j,:,:,day);
                            tmpt(c,:,:)=tlev_days(i,j,:,:,day);
                            tmpq(c,:,:)=qlev_days(i,j,:,:,day);
                        end
                    end
                    twprofile_twhotdays(i,j,:,:)=squeeze(mean(tmptw,'omitnan'));
                    tprofile_twhotdays(i,j,:,:)=squeeze(mean(tmpt,'omitnan'));
                    qprofile_twhotdays(i,j,:,:)=squeeze(mean(tmpq,'omitnan'));
                end
            end
        end
        %Separate into regions
        c=zeros(3,1);twholder=cell(3,1);tholder=cell(3,1);qholder=cell(3,1);
        for i=1:nlat
            for j=1:nlon
                r=NaN;if myregs_era5(i,j)==1;r=1;elseif myregs_era5(i,j)==2;r=2;elseif myregs_era5(i,j)==3;r=3;end
                if ~isnan(r)
                    c(r)=c(r)+1;
                    twholder{r}(c(r),:,:)=twprofile_twhotdays(i,j,:,:);
                    tholder{r}(c(r),:,:)=tprofile_twhotdays(i,j,:,:);
                    qholder{r}(c(r),:,:)=qprofile_twhotdays(i,j,:,:);
                end
            end
        end
        %Get region means and st devs
        twprofile_twhotdays_era5_mean=NaN.*ones(3,6,6);tprofile_twhotdays_era5_mean=NaN.*ones(3,6,6);qprofile_twhotdays_era5_mean=NaN.*ones(3,6,6);
        twprofile_twhotdays_era5_stdev=NaN.*ones(3,6,6);tprofile_twhotdays_era5_stdev=NaN.*ones(3,6,6);qprofile_twhotdays_era5_stdev=NaN.*ones(3,6,6);
        for r=1:3
            twprofile_twhotdays_era5_mean(r,:,:)=squeeze(mean(twholder{r}));
            tprofile_twhotdays_era5_mean(r,:,:)=squeeze(mean(tholder{r}));
            qprofile_twhotdays_era5_mean(r,:,:)=squeeze(mean(qholder{r}));

            twprofile_twhotdays_era5_stdev(r,:,:)=squeeze(std(twholder{r}));
            tprofile_twhotdays_era5_stdev(r,:,:)=squeeze(std(tholder{r}));
            qprofile_twhotdays_era5_stdev(r,:,:)=squeeze(std(qholder{r}));
        end
        save(strcat(saveloc,'hotdaysprofiles_new.mat'),'twprofile_twhotdays_era5_mean','tprofile_twhotdays_era5_mean','qprofile_twhotdays_era5_mean',...
            'twprofile_alldays_era5_mean','tprofile_alldays_era5_mean','qprofile_alldays_era5_mean',...
            'twprofile_twhotdays_era5_stdev','tprofile_twhotdays_era5_stdev','qprofile_twhotdays_era5_stdev',...
            'twprofile_alldays_era5_stdev','tprofile_alldays_era5_stdev','qprofile_alldays_era5_stdev');
    end


    %Stn preparation (only Abu Dhabi -- stn #1)
    stntwarray=squeeze(subdailytw_sPG(1,:,:));
    exist tw95whereabove;
    if ans==0;tmp=load(strcat(saveloc,'hotdaysstarts.mat'));tw95whereabove=tmp.tw95whereabove;end
    c=0;clear stntw;
    for day=1:climodaylen*nyr
        if sum(tw95whereabove(besti_era5(1),bestj_era5(1),:,day))>=1
            c=c+1;
            stntw(c,:)=stntwarray(day,:);
        end
    end
    abudhabitw_twhotdays=mean(stntw,'omitnan');
    abudhabitw_alldays=mean(stntwarray,'omitnan');

    figure(700);clf;subplotc=0;
    horizspacing=0.16;
    lefts=repmat([0.845-horizspacing*5;0.845-horizspacing*4;0.845-horizspacing*3;0.845-horizspacing*2;0.845-horizspacing*1;0.845],4,1);
    bottoms=[repmat(0.75,6,1);repmat(0.51,6,1);repmat(0.27,6,1)];
    w_big=0.14;h_big=0.22;
    cmap_tw=colormaps('whitelightpurpledarkpurple','more','not');twcolor=cmap_tw(round(3*end/4),:);
    cmap_t=colormaps('whitelightorangedarkorange','more','not');tcolor=cmap_t(round(3*end/4),:);
    cmap_q=colormaps('whitelightgreendarkgreen','more','not');qcolor=cmap_q(round(3*end/4),:);
    difffrommax=max(twcolor)-twcolor;twcolor_pale=twcolor+0.4*difffrommax;
    difffrommax=max(tcolor)-tcolor;tcolor_pale=tcolor+0.4*difffrommax;
    difffrommax=max(qcolor)-qcolor;qcolor_pale=qcolor+0.4*difffrommax;
    timelabels={'12am';'4am';'8am';'12pm';'4pm';'8pm'};
    reglabels={'Gulf';'Coast';'Inland'};

    for row=1:3 %regions -- water, coast, inland
        for col=1:6 %times
            subplotc=subplotc+1;
            axes('position',[0.1 0.1 0.1 0.1]);hold on;
            if col>=2;coltopullfrom=col-1;elseif col==1;coltopullfrom=6;end %to account for UTC-LST difference

            allmeantoplot=fliplr(twprofile_alldays_era5_mean(row,2:6,coltopullfrom));
                allmean_1stdevneg=-1.*fliplr(twprofile_alldays_era5_stdev(row,2:6,coltopullfrom));
                allmean_1stdevpos=1.*fliplr(twprofile_alldays_era5_stdev(row,2:6,coltopullfrom));
            hotmeantoplot=fliplr(twprofile_twhotdays_era5_mean(row,2:6,coltopullfrom));
                hotmean_1stdevneg=-1.*fliplr(twprofile_twhotdays_era5_stdev(row,2:6,coltopullfrom));
                hotmean_1stdevpos=1.*fliplr(twprofile_twhotdays_era5_stdev(row,2:6,coltopullfrom));

            %In plotting, leave out 500 mb to better see lower-level behavior
            plot(allmeantoplot,1:5,'color',twcolor_pale,'linestyle',':','linewidth',2);hold on;
            p=errorbar(allmeantoplot,1:5,NaN,NaN,allmean_1stdevneg,allmean_1stdevpos,'linestyle',':','linewidth',1.2,'color',twcolor_pale,'CapSize',1.5);
                for i=1:5;p(i).Bar.LineStyle='dotted';end
            plot(hotmeantoplot,1:5,'color',twcolor,'linestyle','-','linewidth',2); 
            p=errorbar(hotmeantoplot,1:5,NaN,NaN,hotmean_1stdevneg,hotmean_1stdevpos,'linestyle','-','linewidth',1.7,'color',twcolor,'CapSize',2.5);

            xlim([3 32]);
            if col==1
                set(gca,'ytick',1:5,'yticklabel',{'1000';'975';'925';'850';'700'},'ylim',[1 5]);
            else
                set(gca,'ytick',1:5,'yticklabel','','ylim',[1 5]);
            end
            if row==1
                title(timelabels{col},'fontweight','bold','fontname','arial','fontsize',12);
            end
            if row==3 %last row
                set(gca,'xtick',10:10:30);
            else
                set(gca,'xtick',10:10:30,'xticklabel','');
            end
            set(gca,'ticklength',[0.02 0.01])
            set(gca,'fontweight','bold','fontname','arial','fontsize',10);
            set(gca,'position',[lefts(subplotc),bottoms(subplotc),w_big,h_big]);

            if col==1
                t=ylabel(reglabels{row},'fontweight','bold','fontname','arial','fontsize',11);
                tpos=t.Position;set(t,'Position',[tpos(1)+0.2 tpos(2) tpos(3)]);
            end

            %Add RHARM data where available (only Abu Dhabi)
            %Note: current labeled version does not have this
            %radiosonde data included, just because it doesn't really add much
            rharmhotcolor=colors('pink');rharmallcolor=colors('sky blue');
            stnhotcolor=colors('orange');stnallcolor=colors('sirkes_gold');
            correspyidx=[NaN;5;4;3;2.5;2;1.4;1];
            %4am LST
            if row==2 && col==2
                for plotlev=2:size(twprofile_twhotdays_rharm_byloc,2) %i.e. will plot 700, 850, 925, 950, 975, 990, 1000
                    scatter(twprofile_twhotdays_rharm_byloc(1,plotlev,2),correspyidx(plotlev),50,...
                        'marker','pentagram','MarkerFaceColor',rharmhotcolor,'MarkerEdgeColor',rharmhotcolor);
                    scatter(twprofile_alldays_rharm_byloc(1,plotlev,2),correspyidx(plotlev),50,...
                        'marker','pentagram','MarkerFaceColor',rharmallcolor,'MarkerEdgeColor',rharmallcolor);
                end
                %At bottom, add HadISD station
                scatter(abudhabitw_twhotdays(1),1,50,'marker','s','MarkerFaceColor',stnhotcolor,'MarkerEdgeColor',stnhotcolor);
                scatter(abudhabitw_alldays(1),1,50,'marker','s','MarkerFaceColor',stnallcolor,'MarkerEdgeColor',stnallcolor);
            end
            %4pm LST
            if row==2 && col==5
                for plotlev=2:size(twprofile_twhotdays_rharm_byloc,2)
                    scatter(twprofile_twhotdays_rharm_byloc(1,plotlev,5),correspyidx(plotlev),50,...
                        'marker','pentagram','MarkerFaceColor',rharmhotcolor,'MarkerEdgeColor',rharmhotcolor);
                    scatter(twprofile_alldays_rharm_byloc(1,plotlev,5),correspyidx(plotlev),50,...
                        'marker','pentagram','MarkerFaceColor',rharmallcolor,'MarkerEdgeColor',rharmallcolor);
                end
                %At bottom, add HadISD station
                scatter(abudhabitw_twhotdays(13),1,50,'marker','s','MarkerFaceColor',stnhotcolor,'MarkerEdgeColor',stnhotcolor);
                scatter(abudhabitw_alldays(13),1,50,'marker','s','MarkerFaceColor',stnallcolor,'MarkerEdgeColor',stnallcolor);
            end
            %Infill HadISD only for Abu Dhabi/coast:
            if row==2 && col==3 %8am
                scatter(abudhabitw_twhotdays(5),1,50,'marker','s','MarkerFaceColor',stnhotcolor,'MarkerEdgeColor',stnhotcolor);
                scatter(abudhabitw_alldays(5),1,50,'marker','s','MarkerFaceColor',stnallcolor,'MarkerEdgeColor',stnallcolor);
            elseif row==2 && col==4 %12pm
                scatter(abudhabitw_twhotdays(9),1,50,'marker','s','MarkerFaceColor',stnhotcolor,'MarkerEdgeColor',stnhotcolor);
                scatter(abudhabitw_alldays(9),1,50,'marker','s','MarkerFaceColor',stnallcolor,'MarkerEdgeColor',stnallcolor);
            elseif row==2 && col==6 %8pm
                scatter(abudhabitw_twhotdays(17),1,50,'marker','s','MarkerFaceColor',stnhotcolor,'MarkerEdgeColor',stnhotcolor);
                scatter(abudhabitw_alldays(17),1,50,'marker','s','MarkerFaceColor',stnallcolor,'MarkerEdgeColor',stnallcolor);
            elseif row==2 && col==1 %12am
                scatter(abudhabitw_twhotdays(21),1,50,'marker','s','MarkerFaceColor',stnhotcolor,'MarkerEdgeColor',stnhotcolor);
                scatter(abudhabitw_alldays(21),1,50,'marker','s','MarkerFaceColor',stnallcolor,'MarkerEdgeColor',stnallcolor);
            end

            if row==3 && col==3
                t=text(0.24,-0.12,strcat('Wet-bulb Temperature (',char(176),'C)'),'units','normalized');
                set(t,'fontname','arial','fontsize',11,'fontweight','bold');
            end

            %Now add T and q in small plots overlaid!
            %T
            axes('position',[0.1 0.1 0.1 0.1]);hold on;
            allmeantoplot=fliplr(tprofile_alldays_era5_mean(row,2:6,coltopullfrom));
                allmean_1stdevneg=-1.*fliplr(tprofile_alldays_era5_stdev(row,2:6,coltopullfrom));
                allmean_1stdevpos=1.*fliplr(tprofile_alldays_era5_stdev(row,2:6,coltopullfrom));
            hotmeantoplot=fliplr(tprofile_twhotdays_era5_mean(row,2:6,coltopullfrom));
                hotmean_1stdevneg=-1.*fliplr(tprofile_twhotdays_era5_stdev(row,2:6,coltopullfrom));
                hotmean_1stdevpos=1.*fliplr(tprofile_twhotdays_era5_stdev(row,2:6,coltopullfrom));

            plot(allmeantoplot,1:5,'color',tcolor_pale,'linestyle',':','linewidth',1.7);
            p=errorbar(allmeantoplot,1:5,NaN,NaN,allmean_1stdevneg,allmean_1stdevpos,'linestyle',':','linewidth',1.3,'color',tcolor_pale,'CapSize',0);
                for i=1:5;p(i).Bar.LineStyle='dotted';end
            plot(hotmeantoplot,1:5,'color',tcolor,'linestyle','-','linewidth',1.7); 
            p=errorbar(hotmeantoplot,1:5,NaN,NaN,hotmean_1stdevneg,hotmean_1stdevpos,'linestyle','-','linewidth',1.3,'color',tcolor,'CapSize',0);

            set(gca,'ytick',1:5,'yticklabel','','ticklength',[0.04 0.01],'ylim',[1 5]);
            set(gca,'xlim',[15 45],'xtick',20:10:40,'xticklabelrotation',0);
            set(gca,'position',[lefts(subplotc)+0.14*w_big,bottoms(subplotc)+0.09*h_big,w_big/3,h_big/3]);
            set(gca,'fontweight','bold','fontname','arial','fontsize',8);
            if row==1 && col==1;title('T','fontname','arial','fontsize',9,'fontweight','bold');end

            %q
            axes('position',[0.1 0.1 0.1 0.1]);hold on;
            allmeantoplot=fliplr(qprofile_alldays_era5_mean(row,2:6,coltopullfrom));
                allmean_1stdevneg=-1.*fliplr(qprofile_alldays_era5_stdev(row,2:6,coltopullfrom));
                allmean_1stdevpos=1.*fliplr(qprofile_alldays_era5_stdev(row,2:6,coltopullfrom));
            hotmeantoplot=fliplr(qprofile_twhotdays_era5_mean(row,2:6,coltopullfrom));
                hotmean_1stdevneg=-1.*fliplr(qprofile_twhotdays_era5_stdev(row,2:6,coltopullfrom));
                hotmean_1stdevpos=1.*fliplr(qprofile_twhotdays_era5_stdev(row,2:6,coltopullfrom));

            plot(allmeantoplot,1:5,'color',qcolor_pale,'linestyle',':','linewidth',1.7);
            p=errorbar(allmeantoplot,1:5,NaN,NaN,allmean_1stdevneg,allmean_1stdevpos,'linestyle','--','linewidth',1.3,'color',qcolor_pale,'CapSize',0);
                for i=1:5;p(i).Bar.LineStyle='dotted';end
            plot(hotmeantoplot,1:5,'color',qcolor,'linestyle','-','linewidth',1.7); 
            p=errorbar(hotmeantoplot,1:5,NaN,NaN,hotmean_1stdevneg,hotmean_1stdevpos,'linestyle','-','linewidth',1.3,'color',qcolor,'CapSize',0);

            set(gca,'ytick',1:5,'yticklabel','','ticklength',[0.04 0.01],'ylim',[1 5]);
            set(gca,'xlim',[3 26],'xtick',5:10:25,'xticklabelrotation',0);
            set(gca,'position',[lefts(subplotc)+0.66*w_big,bottoms(subplotc)+0.66*h_big,w_big/3,h_big/3]);
            set(gca,'fontweight','bold','fontname','arial','fontsize',8);
            if row==1 && col==1;title('q','fontname','arial','fontsize',9,'fontweight','bold');end

        end
    end
    set(gcf,'color','w');
    figname='fig5_rev';
    curpart=1;highqualityfiguresetup;curpart=2;highqualityfiguresetup;
end





if readutciandwbgtdata==1
    utci_era5=NaN.*ones(size(tw_era5));wbgt_era5=NaN.*ones(size(tw_era5));
    bigc=0;
    for y=2000:2020
        for m=6:8
            utci_tmp=ncread(strcat(era5loc,'utciwbgtoutput_abudhabi_',num2str(y),'0',num2str(m),'.nc'),'utci2m');invalid=abs(utci_tmp)>100;utci_tmp(invalid)=NaN;
            wbgt_tmp=ncread(strcat(era5loc,'utciwbgtoutput_abudhabi_',num2str(y),'0',num2str(m),'.nc'),'wbgt2m');invalid=abs(wbgt_tmp)>100;wbgt_tmp(invalid)=NaN;
            %lat_tmp=ncread(strcat(era5loc,'utcioutput_abudhabi_',num2str(y),'0',num2str(m),'.nc'),'lat');
            %lon_tmp=ncread(strcat(era5loc,'utcioutput_abudhabi_',num2str(y),'0',num2str(m),'.nc'),'lon');

            utci_tmp=permute(utci_tmp,[2 1 3]);
            wbgt_tmp=permute(wbgt_tmp,[2 1 3]);

            %Save into big array
            hoursthismonth=size(utci_tmp,3);
            utci_era5(:,:,bigc+1:bigc+hoursthismonth)=utci_tmp;
            wbgt_era5(:,:,bigc+1:bigc+hoursthismonth)=wbgt_tmp;
            
            bigc=bigc+hoursthismonth;
        end
    end
    save(strcat(saveloc,'utciera5.mat'),'utci_era5');
    save(strcat(saveloc,'wbgtera5.mat'),'wbgt_era5');
end


if readlurompsheatindex==1
    lurompshi_era5=NaN.*ones(size(tw_era5));
    bigc=0;
    for y=2000:2020
        for m=6:8
            lurompshi_tmp=ncread(strcat('/Volumes/ExternalDriveD/Various_LCLUC_arrays/lurompshioutput_abudhabi_',num2str(y),'0',num2str(m),'.nc'),'lurompshi2m');
            invalid=abs(lurompshi_tmp)>100;lurompshi_tmp(invalid)=NaN;

            lurompshi_tmp=permute(lurompshi_tmp,[2 1 3]);

            %Save into big array
            hoursthismonth=size(lurompshi_tmp,3);
            lurompshi_era5(:,:,bigc+1:bigc+hoursthismonth)=lurompshi_tmp;
            
            bigc=bigc+hoursthismonth;
        end
    end
    save(strcat(saveloc,'lurompshiera5.mat'),'lurompshi_era5');
end


if calcheatindex==1
    needtorecalcrh=0;
    if needtorecalcrh==1

        needtorecalctd=0;
        if needtorecalctd==1
            exist t_era5;
            if ans==0
                tmp=load(strcat(saveloc,'tera5.mat'));t_era5=tmp.t_era5;clear tmp; %C
                tmp=load(strcat(saveloc,'qera5.mat'));q_era5=tmp.q_era5;clear tmp; %g/kg
                tmp=load(strcat(saveloc,'psfcera5.mat'));psfc_era5=tmp.psfc_era5;clear tmp; %Pa
            end
            td_era5=calcTdfromq_dynamicP(q_era5,psfc_era5./100);clear q_era5;clear psfc_era5;
            save(strcat(saveloc,'tdera5.mat'),'td_era5');
        end
    
        exist td_era5;
        if ans==0;tmp=load(strcat(saveloc,'tdera5.mat'));td_era5=tmp.td_era5;clear tmp;end %C
    
        rh_era5=calcrhfromTandTd(t_era5,td_era5);clear td_era5; %percent
        save(strcat(saveloc,'rhera5.mat'),'rh_era5');
    end

    exist rh_era5;
    if ans==0;tmp=load(strcat(saveloc,'rhera5.mat'));rh_era5=tmp.rh_era5;clear tmp;end %percent

    %4 min to actually calculate HI, another 5 to save it
    hi_era5=(calcnwsheatindex((t_era5.*1.8)+32,rh_era5)-32).*5/9; %C
    
    clear t_era5;clear rh_era5;
    save(strcat(saveloc,'hiera5.mat'),'hi_era5');
end


if getaddlindexstats==1
    exist hi_era5;if ans==0;tmp=load(strcat(saveloc,'hiera5.mat'));hi_era5=tmp.hi_era5;clear tmp;end
    hi_era5_alltimemax=NaN.*ones(era5latsz,era5lonsz);hi_era5_p95=NaN.*ones(era5latsz,era5lonsz);
    for i=1:era5lonsz
        for j=1:era5latsz
            hi_era5_alltimemax(i,j)=max(hi_era5(i,j,:));
            hi_era5_p95(i,j)=quantile(hi_era5(i,j,:),0.95);
        end
    end

    exist utci_era5;if ans==0;tmp=load(strcat(saveloc,'utciera5.mat'));utci_era5=tmp.utci_era5;clear tmp;end
    utci_era5_alltimemax=NaN.*ones(era5latsz,era5lonsz);utci_era5_p95=NaN.*ones(era5latsz,era5lonsz);
    for i=1:era5lonsz
        for j=1:era5latsz
            utci_era5_alltimemax(i,j)=max(utci_era5(i,j,:));
            utci_era5_p95(i,j)=quantile(utci_era5(i,j,:),0.95);
        end
    end

    exist wbgt_era5;if ans==0;tmp=load(strcat(saveloc,'wbgtera5.mat'));wbgt_era5=tmp.wbgt_era5;clear tmp;end
    wbgt_era5_alltimemax=NaN.*ones(era5latsz,era5lonsz);wbgt_era5_p95=NaN.*ones(era5latsz,era5lonsz);
    for i=1:era5lonsz
        for j=1:era5latsz
            wbgt_era5_alltimemax(i,j)=max(wbgt_era5(i,j,:));
            wbgt_era5_p95(i,j)=quantile(wbgt_era5(i,j,:),0.95);
        end
    end

    exist lurompshi_era5;if ans==0;tmp=load(strcat(saveloc,'lurompshiera5.mat'));lurompshi_era5=tmp.lurompshi_era5;clear tmp;end
    lurompshi_era5_alltimemax=NaN.*ones(era5latsz,era5lonsz);lurompshi_era5_p95=NaN.*ones(era5latsz,era5lonsz);
    for i=1:era5lonsz
        for j=1:era5latsz
            lurompshi_era5_alltimemax(i,j)=max(lurompshi_era5(i,j,:));
            lurompshi_era5_p95(i,j)=quantile(lurompshi_era5(i,j,:),0.95);
        end
    end
    save(strcat(saveloc,'addlstats.mat'),'hi_era5_alltimemax','hi_era5_p95','lurompshi_era5_alltimemax','lurompshi_era5_p95',...
        'utci_era5_alltimemax','utci_era5_p95','wbgt_era5_alltimemax','wbgt_era5_p95');
end

if calcspreads_addlindices==1
    centralhours_hi_era5=NaN.*ones(era5lonsz,era5latsz);hrofmaxp95prob_hi_era5=NaN.*ones(era5lonsz,era5latsz);
    centralhours_utci_era5=NaN.*ones(era5lonsz,era5latsz);hrofmaxp95prob_utci_era5=NaN.*ones(era5lonsz,era5latsz);
    centralhours_wbgt_era5=NaN.*ones(era5lonsz,era5latsz);hrofmaxp95prob_wbgt_era5=NaN.*ones(era5lonsz,era5latsz);
    centralhours_lurompshi_era5=NaN.*ones(era5lonsz,era5latsz);hrofmaxp95prob_lurompshi_era5=NaN.*ones(era5lonsz,era5latsz);

    hi_days=reshape(hi_era5,[nlat nlon 24 climodaylen*nyr]);
    utci_days=reshape(utci_era5,[nlat nlon 24 climodaylen*nyr]);
    wbgt_days=reshape(wbgt_era5,[nlat nlon 24 climodaylen*nyr]);
    lurompshi_days=reshape(lurompshi_era5,[nlat nlon 24 climodaylen*nyr]);

    for i=1:era5lonsz
        for j=1:era5latsz
            %95th pctile across all hours at this gridcell
            clear probp95_hi_era5;
            for hr=1:24;datathishr=squeeze(hi_days(i,j,hr,:));probp95_hi_era5(hr)=sum(datathishr>quantile(hi_era5(i,j,:),0.95))/size(hi_days,4);end
            [~,hrofmaxp95prob_hi_era5(i,j)]=max(probp95_hi_era5);
            %Do at least 50% of p95 values occur within 3 hours of the peak?
            hrmax=hrofmaxp95prob_hi_era5(i,j);
            if hrmax>=22
                tosum=[probp95_hi_era5(hrmax-3:24) probp95_hi_era5(1:hrmax-21)];
            elseif hrmax<=3
                tosum=[probp95_hi_era5(hrmax+21) probp95_hi_era5(1:hrmax+3)];
            else
                tosum=[probp95_hi_era5(hrmax-3:hrmax+3)];
            end
            sum_within3=sum(tosum);
            sum_all=sum(probp95_hi_era5);
            frac_within3=sum_within3/sum_all;
            if frac_within3>=0.5
                centralhours_hi_era5(i,j)=hrmax;
            end
        end
    end
    a=centralhours_hi_era5+tzadj;a(a>=24)=a(a>=24)-24;centralhours_hi_local_era5=a;
    save(strcat(saveloc,'centralvarsera5.mat'),'centralhours_hi_local_era5','-append');


    for i=1:era5lonsz
        for j=1:era5latsz
            %95th pctile across all hours at this gridcell
            clear probp95_utci_era5;
            for hr=1:24;datathishr=squeeze(utci_days(i,j,hr,:));probp95_utci_era5(hr)=sum(datathishr>quantile(utci_era5(i,j,:),0.95))/size(utci_days,4);end
            [~,hrofmaxp95prob_utci_era5(i,j)]=max(probp95_utci_era5);
            %Do at least 50% of p95 values occur within 3 hours of the peak?
            hrmax=hrofmaxp95prob_utci_era5(i,j);
            if hrmax>=22
                tosum=[probp95_utci_era5(hrmax-3:24) probp95_utci_era5(1:hrmax-21)];
            elseif hrmax<=3
                tosum=[probp95_utci_era5(hrmax+21) probp95_utci_era5(1:hrmax+3)];
            else
                tosum=[probp95_utci_era5(hrmax-3:hrmax+3)];
            end
            sum_within3=sum(tosum);
            sum_all=sum(probp95_utci_era5);
            frac_within3=sum_within3/sum_all;
            if frac_within3>=0.5
                centralhours_utci_era5(i,j)=hrmax;
            end
        end
    end
    a=centralhours_utci_era5+tzadj;a(a>=24)=a(a>=24)-24;centralhours_utci_local_era5=a;
    save(strcat(saveloc,'centralvarsera5.mat'),'centralhours_utci_local_era5','-append');

    for i=1:era5lonsz
        for j=1:era5latsz
            %95th pctile across all hours at this gridcell
            clear probp95_wbgt_era5;
            for hr=1:24;datathishr=squeeze(wbgt_days(i,j,hr,:));probp95_wbgt_era5(hr)=sum(datathishr>quantile(wbgt_era5(i,j,:),0.95))/size(wbgt_days,4);end
            [~,hrofmaxp95prob_wbgt_era5(i,j)]=max(probp95_wbgt_era5);
            %Do at least 50% of p95 values occur within 3 hours of the peak?
            hrmax=hrofmaxp95prob_wbgt_era5(i,j);
            if hrmax>=22
                tosum=[probp95_wbgt_era5(hrmax-3:24) probp95_wbgt_era5(1:hrmax-21)];
            elseif hrmax<=3
                tosum=[probp95_wbgt_era5(hrmax+21) probp95_wbgt_era5(1:hrmax+3)];
            else
                tosum=[probp95_wbgt_era5(hrmax-3:hrmax+3)];
            end
            sum_within3=sum(tosum);
            sum_all=sum(probp95_wbgt_era5);
            frac_within3=sum_within3/sum_all;
            if frac_within3>=0.5
                centralhours_wbgt_era5(i,j)=hrmax;
            end
        end
    end
    a=centralhours_wbgt_era5+tzadj;a(a>=24)=a(a>=24)-24;centralhours_wbgt_local_era5=a;
    save(strcat(saveloc,'centralvarsera5.mat'),'centralhours_wbgt_local_era5','-append');

    for i=1:era5lonsz
        for j=1:era5latsz
            %95th pctile across all hours at this gridcell
            clear probp95_lurompshi_era5;
            for hr=1:24;datathishr=squeeze(lurompshi_days(i,j,hr,:));...
                    probp95_lurompshi_era5(hr)=sum(datathishr>quantile(lurompshi_era5(i,j,:),0.95))/size(lurompshi_days,4);end
            [~,hrofmaxp95prob_lurompshi_era5(i,j)]=max(probp95_lurompshi_era5);
            %Do at least 50% of p95 values occur within 3 hours of the peak?
            hrmax=hrofmaxp95prob_lurompshi_era5(i,j);
            if hrmax>=22
                tosum=[probp95_lurompshi_era5(hrmax-3:24) probp95_lurompshi_era5(1:hrmax-21)];
            elseif hrmax<=3
                tosum=[probp95_lurompshi_era5(hrmax+21) probp95_lurompshi_era5(1:hrmax+3)];
            else
                tosum=[probp95_lurompshi_era5(hrmax-3:hrmax+3)];
            end
            sum_within3=sum(tosum);
            sum_all=sum(probp95_lurompshi_era5);
            frac_within3=sum_within3/sum_all;
            if frac_within3>=0.5
                centralhours_lurompshi_era5(i,j)=hrmax;
            end
        end
    end
    a=centralhours_lurompshi_era5+tzadj;a(a>=24)=a(a>=24)-24;centralhours_lurompshi_local_era5=a;
    save(strcat(saveloc,'centralvarsera5.mat'),'centralhours_lurompshi_local_era5','-append');
end

if heatindexcomparison==1
    tarr=0:0.05:50;
    c=1;clear nws_rh75;clear nws_rh40;clear nws_rh15;
    for t=0:0.05:50
        nws_rh75(c)=(calcnwsheatindex((t*1.8)+32,75)-32)*5/9;
        nws_rh40(c)=(calcnwsheatindex((t*1.8)+32,40)-32)*5/9;
        nws_rh15(c)=(calcnwsheatindex((t*1.8)+32,15)-32)*5/9;
        c=c+1;
    end

    luromps_rh75=csvread(strcat(saveloc,'lurompsvalidation_rh75.csv'));
    luromps_rh40=csvread(strcat(saveloc,'lurompsvalidation_rh40.csv'));
    luromps_rh15=csvread(strcat(saveloc,'lurompsvalidation_rh15.csv'));

    figure(1);clf;
    subplot(10,10,100);hold on;plot(tarr,nws_rh75,'k','linewidth',2);plot(tarr(2:end-1),luromps_rh75,'color',colors('orange'),'linewidth',2);xlim([0 50]);
        title('RH=75%','fontsize',14,'fontweight','bold','fontname','arial');set(gca,'fontsize',12,'fontweight','bold','fontname','arial');
        set(gca,'position',[0.05 0.7 0.5 0.26]);
    subplot(3,1,2);hold on;plot(tarr,nws_rh40,'k','linewidth',2);plot(tarr(2:end-1),luromps_rh40,'color',colors('orange'),'linewidth',2);xlim([0 50]);
        title('RH=40%','fontsize',14,'fontweight','bold','fontname','arial');set(gca,'fontsize',12,'fontweight','bold','fontname','arial');
        set(gca,'position',[0.05 0.37 0.5 0.26]);
    subplot(3,1,3);hold on;plot(tarr,nws_rh15,'k','linewidth',2);plot(tarr(2:end-1),luromps_rh15,'color',colors('orange'),'linewidth',2);xlim([0 50]);
        title('RH=15%','fontsize',14,'fontweight','bold','fontname','arial');set(gca,'fontsize',12,'fontweight','bold','fontname','arial');
        set(gca,'position',[0.05 0.04 0.5 0.26]);

    figname='heatindexcomparison';curpart=1;highqualityfiguresetup;curpart=2;highqualityfiguresetup;
    
    set(gcf,'color','w');
end

if makefig1_suppversion_addlindices==1
    figure(900);clf;

    %HI in top row, UTCI in middle row, WBGT in bottom row
    clear datatoplot;
    datatoplot{1}{1}=hi_era5_p95;datatoplot{1}{2}=centralhours_hi_local_era5;
    datatoplot{2}{1}=utci_era5_p95;datatoplot{2}{2}=centralhours_utci_local_era5;
    datatoplot{3}{1}=wbgt_era5_p95;datatoplot{3}{2}=centralhours_wbgt_local_era5;
    lefts=[0.12;0.52;0.12;0.52;0.12;0.52];bottoms=[0.67;0.67;0.345;0.345;0.02;0.02];
    cmins=[42;1;40;1;28;1];cmaxs=[52;24;54;24;35;24];
    cmaps={colormaps('wbt','more','not');colormaps('mygbm','24','not');colormaps('wbt','more','not');colormaps('mygbm','24','not');...
        colormaps('wbt','more','not');colormaps('mygbm','24','not')};
    stepsizes=[0.5;1;0.7;1;0.5;1];

    c=1;
    for row=1:3
        for col=1:2
            subplot(10,10,100);
            data={lat_era5;lon_era5;datatoplot{row}{col}};
            cmin=cmins(c);cmax=cmaxs(c);cmap=cmaps{c};
            step=stepsizes(c);

            vararginnew={'mapproj';'mercator';'datatounderlay';data;'underlaycaxismin';cmin;'underlaycaxismax';cmax;'underlaystepsize';step;'underlaycolormap';cmap;
                'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
                'overlaynow';0;'conttoplot';'all';'nonewfig';1;'omitfirstsubplotcolorbar';0;...
                'colorbarfontsize';14};
            datatype='custom';region={wb;nb;eb;sb};plotModelData(data,region,vararginnew,datatype);

            set(gca,'position',[lefts(c) bottoms(c) 0.45 0.31]);

            c=c+1;
        end
    end


    figname=strcat('hiutciwbgtp95_withera5');
    set(gcf,'color','w');curpart=1;highqualityfiguresetup;curpart=2;highqualityfiguresetup;
end

clear t2m_days;clear q2m_days;


%Timing of wind shift on p95-Tw days in ERA5 vs stns, and across space (including different places along coast)
%panel a: hr of qmin in ERA5
%panel b: stns
%panel c: hr of onset of northerly winds in ERA5
%panel d: stns
%Takeaways: very sharp geographic difference is hard to get precisely right -- for
%example, comparing panels a and b of this figure
    %Also, none of the official stns in Fig 1c is far off, actually, though at first glance they may appear so
    
if suppfigwinddirs==1
    if moresuppcalcs==1
        %Some exploratory stuff
        subdailywinddir_stns=subdailywinddir_sPG;
        subdailytw_stns=subdailytw_sPG;
        subdailyt_stns=subdailyt_sPG;
        subdailytd_stns=subdailytd_sPG;
        for s=1:4
            winddircomposite=zeros(1932,5);
            for d=1:1932
                [peakval,peakhr]=max(subdailytw_stns(s,d,:));
                if peakhr>=3 && peakhr<=22
                    winddircomposite(d,:)=round2(subdailywinddir_stns(s,d,peakhr-2:peakhr+2),30);
                end
            end
            modewinddircomp_stns(s,:)=mode(winddircomposite);
        end
    
        test=squeeze(mode(subdailywinddir_stns,2));imsq(test);
    
    
        %Organize relevant station data
        clear winddirvec;clear windspdvec;clear twvec;clear tvec;clear qvec;clear allstnlats;clear allstnlons;stnc=0;
        for reg=1:3
            r=regnames{reg};
            twp95vec=eval(['tw_' r 'stns_p95;']);subdailytw=eval(['subdailytw_' r ';']);
            subdailyt=eval(['subdailyt_' r ';']);subdailytd=eval(['subdailytd_' r ';']);
            subdailywinddir=eval(['subdailywinddir_' r ';']);subdailywindspd=eval(['subdailywindspd_' r ';']);
            stninfo=eval(['stninfo_' r ';']);
            for stn=1:length(twp95vec)
                stnc=stnc+1;
                allstnlats(stnc)=stninfo(stn,1);allstnlons(stnc)=stninfo(stn,2);
                thisp95=twp95vec(stn);dayc=0;
                for hr=1:24
                    for day=1:size(subdailytw,2)
                        if subdailytw(stn,day,hr)>thisp95
                            dayc=dayc+1;
                            winddirvec{stnc}(dayc,:)=subdailywinddir(stn,day,:);
                            windspdvec{stnc}(dayc,:)=subdailywindspd(stn,day,:);
                            twvec{stnc}(dayc,:)=subdailytw(stn,day,:);
                            tvec{stnc}(dayc,:)=subdailyt(stn,day,:);
                            qvec{stnc}(dayc,:)=calcqfromTd(subdailytd(stn,day,:));
                        end
                    end
                end
            end
        end
    
        %Compute hr of min q/Td on p95-Tw days from stations
        clear hrofqmin_stns;
        for stn=1:length(qvec)
            [~,hrofqmin_stns(stn)]=min(mean(qvec{stn},'omitnan'));
        end
    
        %Compute median hr of onset of wind dir near-northerly (300-360-60) on p95-Tw days from stations
        clear hrofwindshift_stns;
        for stn=1:length(winddirvec)
            winddirs=winddirvec{stn};
            hrofwindshift_tmp=NaN.*ones(size(winddirs,1),1);
            for day=1:size(winddirs,1)
                keepgoing=1;
                [~,thisqminhr]=min(qvec{stn}(day,:));
                if thisqminhr>=3
                    for hr=thisqminhr-2:24 %expect answer to be around or just after the day's hrofqmin
                        if (winddirs(day,hr)>=300 || winddirs(day,hr)<=60) && keepgoing==1
                            hrofwindshift_tmp(day)=hr;keepgoing=0;
                        end
                    end
                end
            end
            hrofwindshift_stns(stn)=median(hrofwindshift_tmp,'omitnan');
        end
        
        %Repeat for ERA5 gridcells (15 sec)
        twholder=cell(3,1);tholder=cell(3,1);qholder=cell(3,1);winddirholder=cell(3,1);windspdholder=cell(3,1);c=zeros(3,1);
        peaktwhrbypt=NaN.*ones(nlon,nlat);hrofqmin=NaN.*ones(nlon,nlat);hrofwindshift=NaN.*ones(nlon,nlat);
        for i=1:nlon
            for j=1:nlat
                r=NaN;
                if myregs_era5(i,j)==1
                    r=1;
                elseif myregs_era5(i,j)==2
                    r=2;
                elseif myregs_era5(i,j)==3
                    r=3;
                end
                if ~isnan(r)
                    c(r)=c(r)+1;
                    twholder{r}(c(r),:)=tw2mp95_diurnalmean(i,j,:);
                    tholder{r}(c(r),:)=t2mp95_diurnalmean(i,j,:);
                    qholder{r}(c(r),:)=q2mp95_diurnalmean(i,j,:);
                    winddirholder{r}(c(r),:)=winddir10p95_diurnalmode(i,j,:);
                    windspdholder{r}(c(r),:)=windspd10p95_diurnalmean(i,j,:);
                end
                [maxtwval,peaktwhrbypt(i,j)]=max(tw2mp95_diurnalmean(i,j,:));
                [~,hrofqmin(i,j)]=min(q2mp95_diurnalmean(i,j,:)); %hr of min q on p95-Tw days
    
                %Compute mean hr of onset of wind dir near-northerly (300-360-60) on p95-Tw days
                winddirs=winddir10p95arr{i,j}; 
                hrofwindshift_tmp=NaN.*ones(size(winddirs,1),1);
                for day=1:size(winddirs,1)
                    keepgoing=1;
                    [~,thisqminhr]=min(q2mp95arr{i,j}(day,:));
                    if thisqminhr>=3
                        for hr=thisqminhr-2:24 %expect answer to be around or just after the day's hrofqmin
                            if (winddirs(day,hr)>=300 || winddirs(day,hr)<=60) && keepgoing==1
                                hrofwindshift_tmp(day)=hr;keepgoing=0;
                            end
                        end
                    end
                end
                hrofwindshift(i,j)=median(hrofwindshift_tmp,'omitnan');
            end
        end
        %Spatial means for each regions
        for r=1:3
            meantw_era5(r,:)=mean(twholder{r});
            meant_era5(r,:)=mean(tholder{r});
            meanq_era5(r,:)=mean(qholder{r});
            meanwinddir_era5(r,:)=mean(winddirholder{r});
            meanwindspd_era5(r,:)=mean(windspdholder{r});
        end
    
        %Convert arrays to be plotted to LST
        hrofqmin_stns_lst=NaN.*ones(length(allstnlats),1);
        hrofwindshift_stns_lst=NaN.*ones(length(allstnlats),1);
        for stn=1:length(allstnlats)
            stnhour=hrofqmin_stns(stn);
            stnhour_lst=stnhour+tzadj;if stnhour_lst>=24;stnhour_lst=stnhour_lst-24;end
            hrofqmin_stns_lst(stn)=stnhour_lst;
    
            stnhour=hrofwindshift_stns(stn);
            stnhour_lst=stnhour+tzadj;if stnhour_lst>=24;stnhour_lst=stnhour_lst-24;end
            hrofwindshift_stns_lst(stn)=stnhour_lst;
        end
    
        hrofqmin_lst=hrofqmin+tzadj;toadjust=hrofqmin_lst>=24;hrofqmin_lst(toadjust)=hrofqmin_lst(toadjust)-24;
        hrofwindshift_lst=hrofwindshift+tzadj;toadjust=hrofwindshift_lst>=24;hrofwindshift_lst(toadjust)=hrofwindshift_lst(toadjust)-24;
    end

    %Plot
    if plotwindshiftmap==1
        figure(345);clf;

        cmin=9;cmax=16; %focus on middle of day
        if cmax-cmin>=20
            mintick=4;tickstep=4;maxtick=24;
        elseif cmin==9 && cmax==16
            mintick=10;tickstep=2;maxtick=16;
        end

        %ERA5: hr of q min of p95 days
        ax=subplot(10,10,100);gca=ax;
        data={lat_era5;lon_era5;hrofqmin_lst};
        cmap=colormaps('mygbm','more','not');
        vararginnew={'mapproj';'mercator';'datatounderlay';data;'underlaycaxismin';cmin;'underlaycaxismax';cmax;'underlaystepsize';1;'underlaycolormap';cmap;
            'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
            'overlaynow';0;'conttoplot';'all';'nonewfig';1;'omitfirstsubplotcolorbar';0;...
            'colorbarfontsize';12;'colorbarticks';[mintick:tickstep:maxtick]};
        datatype='custom';region={wb;nb-0.95;eb-2.5;sb+1.45};plotModelData(data,region,vararginnew,datatype);
        set(gca,'position',[0.16 0.52 0.47 0.45]);
        t=text(-0.05,0.51,splabels{1},'units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',12);


        %Stns: hr of q min of p95 days
        ax=subplot(10,10,100);gca=ax;
        data={lat_era5;lon_era5;NaN.*ones(nlon,nlat)};
        cmap=colormaps('mygbm','more','not');
        vararginnew={'mapproj';'mercator';'datatounderlay';data;'underlaycaxismin';cmin;'underlaycaxismax';cmax;'underlaystepsize';1;'underlaycolormap';cmap;
            'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
            'overlaynow';0;'conttoplot';'all';'nonewfig';1;'omitfirstsubplotcolorbar';0;...
            'colorbarfontsize';12;'colorbarticks';[mintick:tickstep:maxtick]};
        datatype='custom';region={wb;nb-1;eb-2.5;sb+1.5};plotModelData(data,region,vararginnew,datatype);hold on;

        cmap=colormaps('mygbm','24','not');
        for stn=1:length(allstnlats)
            stnhour=hrofqmin_stns_lst(stn);
            if ~isnan(stnhour)
                if stnhour<cmin
                    valcolor=cmap(1,:);
                elseif stnhour>cmax
                    valcolor=cmap(end,:);
                else
                    valcolor=cmap(round(size(cmap,1)*(stnhour-cmin)/(cmax-cmin)),:);
                end
            end
            geoshow(allstnlats(stn),allstnlons(stn),'DisplayType','Point','Marker','o','MarkerFaceColor',valcolor,'MarkerEdgeColor','k','MarkerSize',8,'linewidth',1.5);
        end
        set(gca,'position',[0.52 0.52 0.47 0.45]);
        t=text(-0.05,0.51,splabels{2},'units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',12);

        

        %ERA5: hr of shift to northerly winds on p95 days
        ax=subplot(10,10,100);gca=ax;
        data={lat_era5;lon_era5;hrofwindshift_lst};
        cmap=colormaps('mygbm','more','not');
        vararginnew={'mapproj';'mercator';'datatounderlay';data;'underlaycaxismin';cmin;'underlaycaxismax';cmax;'underlaystepsize';1;'underlaycolormap';cmap;
            'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
            'overlaynow';0;'conttoplot';'all';'nonewfig';1;'omitfirstsubplotcolorbar';0;...
            'colorbarfontsize';12;'colorbarticks';[mintick:tickstep:maxtick]};
        datatype='custom';region={wb;nb-0.95;eb-2.5;sb+1.45};plotModelData(data,region,vararginnew,datatype);
        set(gca,'position',[0.16 0.02 0.47 0.45]);
        t=text(-0.05,0.51,splabels{3},'units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',12);


        %Stns: hr of shift to northerly winds on p95 days
        ax=subplot(10,10,100);gca=ax;
        data={lat_era5;lon_era5;NaN.*ones(nlon,nlat)};
        cmap=colormaps('mygbm','more','not');
        vararginnew={'mapproj';'mercator';'datatounderlay';data;'underlaycaxismin';cmin;'underlaycaxismax';cmax;'underlaystepsize';1;'underlaycolormap';cmap;
            'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
            'overlaynow';0;'conttoplot';'all';'nonewfig';1;'omitfirstsubplotcolorbar';0;...
            'colorbarfontsize';12;'colorbarticks';[mintick:tickstep:maxtick]};
        datatype='custom';region={wb;nb-1;eb-2.5;sb+1.5};plotModelData(data,region,vararginnew,datatype);hold on;

        
        cmap=colormaps('mygbm','24','not');
        for stn=1:length(allstnlats)
            stnhour=hrofwindshift_stns_lst(stn);
            if ~isnan(stnhour)
                if stnhour<cmin
                    valcolor=cmap(1,:);
                elseif stnhour>cmax
                    valcolor=cmap(end,:);
                else
                    valcolor=cmap(round(size(cmap,1)*(stnhour-cmin)/(cmax-cmin)),:);
                end
            end
            geoshow(allstnlats(stn),allstnlons(stn),'DisplayType','Point','Marker','o','MarkerFaceColor',valcolor,'MarkerEdgeColor','k','MarkerSize',8,'linewidth',1.5);
        end
        set(gca,'position',[0.52 0.02 0.47 0.45]);
        t=text(-0.05,0.51,splabels{4},'units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',12);



        figname='windshiftmap_NEW';
        set(gcf,'color','w');curpart=1;highqualityfiguresetup;curpart=2;highqualityfiguresetup;
    end

    %More detailed Abu Dhabi-Dubai comparison
    if plotabudhabivsdubaidetail==1
        figure(77);clf;

        %Abu Dhabi: stn 1, Dubai: stn 2
        colo{1}=colors('moderate dark blue');colo{2}=colors('sand');lw=1.5;
        toplot={twvec;0;winddirvec;windspdvec;tvec;qvec};
        ylabs={'deg C';'';'deg';'m/s';'deg C';'g/kg';'percent';'% difference'};

        %Single closest ERA5 gridpt for each
        era5toplot={tw2mp95arr;0;winddir10p95arr;windspd10p95arr;t2mp95arr;q2mp95arr};


        for loop=1:6
            if loop==1;loopforsp=loop;else;loopforsp=loop-1;end
            if loop~=2
                subplot(4,2,loop);hold on;
                for i=1:2 %cities -- Abu Dhabi, Dubai
                    if i==3 %wind dir
                        a=mode(toplot{loop}{i});
                    else
                        a=mean(toplot{loop}{i},'omitnan');
                    end
                    b=[a(21:24) a(1:20)]; %LST conversion
                    plot(b,'linewidth',lw,'color',colo{i});
                    set(gca,'xlim',[1 24],'fontsize',10,'fontweight','bold','fontname','arial');
                    if loop==3;set(gca,'ylim',[0 360],'ytick',[120:120:360]);end

                    %Add ERA5 comparison
                    outarr=interpolate2dlatlonarray(allstnlats(i),allstnlons(i),lat_era5,lon_era5,'-180-180');
                    if i==3 %wind dir
                        era5data=mode(squeeze(era5toplot{loop}{outarr(1,1),outarr(1,2)}));
                    else
                        era5data=mean(squeeze(era5toplot{loop}{outarr(1,1),outarr(1,2)}),'omitnan');
                    end
                    era5datb=[era5data(21:24) era5data(1:20)]; %LST conversion
                    plot(era5datb,'linewidth',lw,'color',colo{i},'linestyle','--');
                end
                t=text(-0.07,0.51,splabels{loopforsp},'units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',12);
                ylabel(ylabs{loop},'fontname','arial','fontweight','bold','fontsize',11);box on;
            end
        end

        subplot(4,2,7);hold on;
        for i=1:2
            rh_city=calcrhfromTandTd(mean(tvec{i},'omitnan'),calcTdfromq(mean(qvec{i},'omitnan')));
            plot([rh_city(21:24) rh_city(1:20)],'linewidth',lw,'color',colo{i});

            outarr=interpolate2dlatlonarray(allstnlats(i),allstnlons(i),lat_era5,lon_era5,'-180-180');
            rh_city_era5=calcrhfromTandTd(mean(t2mp95arr{outarr(1,1),outarr(1,2)},'omitnan'),...
                calcTdfromq(mean(q2mp95arr{outarr(1,1),outarr(1,2)},'omitnan')));
            plot([rh_city_era5(21:24) rh_city_era5(1:20)],'linewidth',lw,'color',colo{i},'linestyle','--');
        end
        set(gca,'xlim',[1 24],'fontsize',10,'fontweight','bold','fontname','arial');
        t=text(-0.07,0.51,splabels{6},'units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',12);box on;
        ylabel(ylabs{7},'fontname','arial','fontweight','bold','fontsize',11);
        xlabel('Hour (LST)','fontname','arial','fontweight','bold','fontsize',11);

        %Accumulated northerly wind 'marine-air transport'
        %Starts from 05 UTC, which is before the shift happens in these cities
        marineaircdf=NaN.*ones(2,24);marineaircdf_era5=NaN.*ones(2,24);
        for i=1:2
            accum=0;accum_era5=0;
            winddirarr=mode(winddirvec{i});windspdarr=mean(windspdvec{i},'omitnan');
            for hr=6:22 %UTC time
                if winddirarr(hr)>=300 || winddirarr(hr)<=60 %northerly wind
                    accum=accum+windspdarr(hr);
                    marineaircdf(i,hr)=accum;
                end
            end

            outarr=interpolate2dlatlonarray(allstnlats(i),allstnlons(i),lat_era5,lon_era5,'-180-180');
            winddirarr_era5=mode(winddir10p95arr{outarr(1,1),outarr(1,2)});
            windspdarr_era5=mean(windspd10p95arr{outarr(1,1),outarr(1,2)},'omitnan');
            for hr=6:22 %UTC time
                if winddirarr_era5(hr)>=300 || winddirarr_era5(hr)<=60 %northerly wind
                    accum_era5=accum_era5+windspdarr_era5(hr);
                    marineaircdf_era5(i,hr)=accum_era5;
                end
            end
        end
        %Plot relative diff (Abu Dhabi vs Dubai)
        subplot(4,2,8);hold on;

        toosmalltoinclude=marineaircdf<5;marineaircdf(toosmalltoinclude)=NaN; %impose minimum of 5 m/s*hr to avoid misleadingly large pct diffs
        toosmalltoinclude=marineaircdf_era5<5;marineaircdf_era5(toosmalltoinclude)=NaN;

        reldiff=NaN.*ones(24,1);
        for hr=1:24;reldiff(hr)=100.*(marineaircdf(1,hr)-marineaircdf(2,hr))./marineaircdf(2,hr);end
        plot([reldiff(21:24);reldiff(1:20)],'color',colors('crimson'),'linewidth',lw);
        reldiff_era5=NaN.*ones(24,1);
        for hr=1:24;reldiff_era5(hr)=100.*(marineaircdf_era5(1,hr)-marineaircdf_era5(2,hr))./marineaircdf_era5(2,hr);end
        plot([reldiff_era5(21:24);reldiff_era5(1:20)],'color',colors('crimson'),'linewidth',lw,'linestyle','--');

        set(gca,'xlim',[1 24],'ylim',[-30 30],'fontsize',10,'fontweight','bold','fontname','arial');
        z=zeros(24,1);plot(z,'k--','linewidth',lw);
        t=text(-0.07,0.51,splabels{7},'units','normalized');set(t,'fontname','arial','fontweight','bold','fontsize',12);box on;
        ylabel(ylabs{8},'fontname','arial','fontweight','bold','fontsize',11);
        xlabel('Hour (LST)','fontname','arial','fontweight','bold','fontsize',11);
        
        set(gcf,'color','w');
        figname='abudhabivsdubai_detail';curpart=1;highqualityfiguresetup;curpart=2;highqualityfiguresetup;
    end

    %Bonus: map of diurnal q range
    if plotdiurnalqrange==1
        clear diurnalrange;
        for i=1:161;for j=1:161;diurnalrange(i,j)=max(q2mp95_diurnalmean(i,j,:))-min(q2mp95_diurnalmean(i,j,:));end;end
        figure(900);clf;
        data={lat_era5;lon_era5;diurnalrange};
        cmap=colormaps('whitelightpurpledarkpurple','more','not');cmin=0;cmax=12;
        vararginnew={'mapproj';'mercator';'datatounderlay';data;'underlaycaxismin';cmin;'underlaycaxismax';cmax;'underlaystepsize';1;'underlaycolormap';cmap;
            'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
            'overlaynow';0;'conttoplot';'all';'nonewfig';1;'omitfirstsubplotcolorbar';0;...
            'colorbarfontsize';12;'colorbarticks';[0:2:12]};
        datatype='custom';region={wb;nb-0.95;eb-2.5;sb+1.45};plotModelData(data,region,vararginnew,datatype);
        set(gcf,'color','w');
        figname='diurnalqrange';curpart=1;highqualityfiguresetup;curpart=2;highqualityfiguresetup;
    end
end
