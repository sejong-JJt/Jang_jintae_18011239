function gps_DCM=ECI2ECEF_DCM_GPS(time_gps) %[YYYY,MM,DD,hh,mm,ss]
time_dt_gps=datetime(time_gps);
time_jd_gps=juliandate(time_dt_gps);
time_GMST_gps=siderealTime(time_jd_gps);

gps_UT=[2023, 4, 29, 0, 0, 0]; %Universal time on Apr 29, 2023
gps_UT_dt=datetime(gps_UT);
gps_UT_jd=juliandate(gps_UT_dt);
gps_UT_GMST=siderealTime(gps_UT_jd); 

w_earth=360/86160; %earth angular rate [deg/s], sidereal day in seconds

gps_t_t0=etime(time_gps,gps_UT); %[s]

gps_theta_g=gps_UT_GMST+w_earth*gps_t_t0;

gps_DCM=[cosd(gps_theta_g) sind(gps_theta_g) 0; -sind(gps_theta_g) cosd(gps_theta_g) 0; 0 0 1];
