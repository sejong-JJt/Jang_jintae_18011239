function bds_DCM=ECI2ECEF_DCM_BDS(time) %[YYYY,MM,DD,hh,mm,ss]
time_dt=datetime(time);
time_jd=juliandate(time_dt);
time_GMST=siderealTime(time_jd);

bds_UT=[2023, 5, 9, 2, 0, 0]; %Universal time on Apr 29, 2023
bds_UT_dt=datetime(bds_UT);
bds_UT_jd=juliandate(bds_UT_dt);
bds_UT_GMST=siderealTime(bds_UT_jd); 

w_earth=360/86160; %earth angular rate [deg/s], sidereal day in seconds

bds_t_t0=etime(time,bds_UT); %[s]

bds_theta_g=bds_UT_GMST+w_earth*bds_t_t0;

bds_DCM=[cosd(bds_theta_g) sind(bds_theta_g) 0; -sind(bds_theta_g) cosd(bds_theta_g) 0; 0 0 1];
