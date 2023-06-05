function DCM=ECI2ECEF_DCM(time) %[YYYY,MM,DD,hh,mm,ss]
time_dt=datetime(time);
time_jd=juliandate(time_dt);
time_GMST=siderealTime(time_jd);

UT=[2023, 4, 29, 0, 0, 0]; %Universal time on Apr 29, 2023
UT_dt=datetime(UT);
UT_jd=juliandate(UT_dt);
UT_GMST=siderealTime(UT_jd); 

w_earth=360/86160; %earth angular rate [deg/s], sidereal day in seconds

t_t0=etime(time,UT); %[s]

theta_g=UT_GMST+w_earth*t_t0;

DCM=[cosd(theta_g) sind(theta_g) 0; -sind(theta_g) cosd(theta_g) 0; 0 0 1];
