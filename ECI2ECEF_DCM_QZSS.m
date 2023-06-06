function qzss_DCM=ECI2ECEF_DCM_QZSS(time) %[YYYY,MM,DD,hh,mm,ss]
time_dt=datetime(time);
time_jd=juliandate(time_dt);
time_GMST=siderealTime(time_jd);

qzss_UT=[2023, 5, 9, 2, 0, 0]; %Universal time on Apr 29, 2023
qzss_UT_dt=datetime(qzss_UT);
qzss_UT_jd=juliandate(qzss_UT_dt);
qzss_UT_GMST=siderealTime(qzss_UT_jd); 

w_earth=360/86160; %earth angular rate [deg/s], sidereal day in seconds

qzss_t_t0=etime(time,qzss_UT); %[s]

qzss_theta_g=qzss_UT_GMST+w_earth*qzss_t_t0;

qzss_DCM=[cosd(qzss_theta_g) sind(qzss_theta_g) 0; -sind(qzss_theta_g) cosd(qzss_theta_g) 0; 0 0 1];
