function qzss_DCM=ECI2ECEF_DCM_QZSS(time_qzss) %[YYYY,MM,DD,hh,mm,ss]
time_dt_qzss=datetime(time_qzss);
time_jd_qzss=juliandate(time_dt_qzss);
time_GMST_qzss=siderealTime(time_jd_qzss);

qzss_UT=[2023, 5, 9, 2, 0, 0]; %Universal time on Apr 29, 2023
qzss_UT_dt=datetime(qzss_UT);
qzss_UT_jd=juliandate(qzss_UT_dt);
qzss_UT_GMST=siderealTime(qzss_UT_jd); 

w_earth=360/86160; %earth angular rate [deg/s], sidereal day in seconds

qzss_t_t0=etime(time_qzss,qzss_UT); %[s]

qzss_theta_g=qzss_UT_GMST+w_earth*qzss_t_t0;

qzss_DCM=[cosd(qzss_theta_g) sind(qzss_theta_g) 0; -sind(qzss_theta_g) cosd(qzss_theta_g) 0; 0 0 1];
