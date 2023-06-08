function bds_DCM=ECI2ECEF_DCM_BDS(time_bds) %[YYYY,MM,DD,hh,mm,ss]
time_dt_bds=datetime(time_bds);
time_jd_bds=juliandate(time_dt_bds);
time_GMST_bds=siderealTime(time_jd_bds);

bds_UT=[2023, 2, 26, 0, 0, 0]; %Universal time on Apr 29, 2023
bds_UT_dt=datetime(bds_UT);
bds_UT_jd=juliandate(bds_UT_dt);
bds_UT_GMST=siderealTime(bds_UT_jd); 

w_earth=360/86160; %earth angular rate [deg/s], sidereal day in seconds

bds_t_t0=etime(time_bds,bds_UT); %[s]

bds_theta_g=bds_UT_GMST+w_earth*bds_t_t0;

bds_DCM=[cosd(bds_theta_g) sind(bds_theta_g) 0; -sind(bds_theta_g) cosd(bds_theta_g) 0; 0 0 1];
