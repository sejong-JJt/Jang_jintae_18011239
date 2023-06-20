close
clear
clc

load nav.mat
gps_M=zeros(1,1440);
gps_a=getfield(nav.GPS,'a');           %Semi-major axis [m]
gps_e=getfield(nav.GPS,'e');           %Eccentricity
gps_i=getfield(nav.GPS,'i');           %Inclination [rad]
gps_omega=getfield(nav.GPS,'omega');   %Argument of perigee,w [rad]
gps_M(1)=getfield(nav.GPS,'M0');       %Mean anomaly at toc [rad]
gps_toc=getfield(nav.GPS,'toc');       %Epoch,current time
gps_OMEGA=getfield(nav.GPS,'OMEGA');   %RAAN [rad]

startTm_gps=gps_toc;
startTime_gps=datetime(gps_toc);
w_earth=360/86160; %earth angular rate [deg/s]

gps_UT=[2023, 4, 29, 0, 0, 0]; %Universal time on Apr 29, 2023
gps_UT_dt=datetime(gps_UT);
gps_UT_jd=juliandate(gps_UT_dt);
gps_UT_GMST=siderealTime(gps_UT_jd); 

gps_ENU=[0 0 0];

u=3.986004418*10^14; %[m^3 s^-2]

gps_p=gps_a*(1-gps_e^2); %[m]

%ground station
gs_lat=37;      %[deg]
gs_lon=127;     %[deg]
gs_h=1000;      %[m]
el_mask=10;     %[deg]

i=1;
gps_E(1)=gps_M(1);
while (1)
    gps_E(i+1)=gps_M(1)+gps_e*sin(gps_E(i));
    if abs(gps_E(i+1)-gps_E(i)) < gps_e
        break
    end
    i=i+1;
end

gps_E(1)=gps_E(end);

for t=1:1:1440 %24hour minutes
    %new GMST
    gps_DCM=ECI2ECEF_DCM_GPS(gps_toc);
    
    %find true_anomaly at t
    gps_y(t)=(sqrt(1-gps_e^2)*sin(gps_E(t)))/(1-gps_e*cos(gps_E(t)));
    gps_x(t)=(cos(gps_E(t))-gps_e)/(1-gps_e*cos(gps_E(t)));

    gps_ta(t)=mod(atan2(gps_y(t),gps_x(t)),2*pi); %true_anomaly [rad]
    
    %PQW range & velocity with true_anomaly
    gps_r(t)=gps_p/(1+gps_e*cos(gps_ta(t))); %[m]
    gps_rangeInPQW = solve_rad_RangeInPerifocalFrame(gps_a, gps_e, gps_ta(t));
    gps_velocityInPQW = solve_rad_VelocityInPerifocalFrame(gps_a, gps_e, gps_ta(t));
    
    %change PQW -> ECI
    gps_C_pqw2eci = PQW2ECI(gps_omega, gps_i, gps_OMEGA);
    gps_rangeInECI=gps_C_pqw2eci*gps_rangeInPQW;
    gps_velocityInECI=gps_C_pqw2eci*gps_velocityInPQW;
    
    %ECI -> ECEF
    gps_rangeInECEF=gps_DCM*gps_rangeInECI; %[m]
    gps_velocityInECEF=gps_DCM*gps_velocityInECI; %[m/s]
    
    %ECEF-> geodetic lat,long
    earth=wgs84Ellipsoid('meter');
    [gps_lat(t), gps_lon(t), gps_h(t)]=ecef2geodetic(earth, gps_rangeInECEF(1), gps_rangeInECEF(2), gps_rangeInECEF(3)); %geodetic lat,long [deg,deg,m]
    
    
    %ENU -> Az & El
    [gps_east(t), gps_north(t), gps_up(t)]=ecef2enu(gps_rangeInECEF(1), gps_rangeInECEF(2), gps_rangeInECEF(3),gs_lat,gs_lon,gs_h,earth);
    gps_ENU=[gps_east(t), gps_north(t), gps_up(t)];
     
    gps_az(t)=azimuth(gps_ENU); %[deg]
    gps_el(t)=elevation(gps_ENU, el_mask); %[deg]
    

    %time pass "1 minute", at t=1-> 00:00:00/ t=2-> 00:01:00
    gps_toc(1,5)=gps_toc(1,5)+1;
    if gps_toc(1,5)>=60
        gps_toc(1,5)=gps_toc(1,5)-60;
        gps_toc(1,4)=gps_toc(1,4)+1;
    end

    if gps_toc(1,4)>=24
        gps_toc(1,4)=gps_toc(1,4)-24;
        gps_toc(1,3)=gps_toc(1,3)+1;
    end

    %find Mean & Eccentricity anomaly at t+1
    gps_t_t0=etime(gps_toc,startTm_gps); %t0=startTime
    gps_tau=sqrt(gps_a^3/u)*2*pi;
    gps_k=fix(gps_t_t0/gps_tau); 
    gps_teM(t+1)=gps_M(1)+sqrt(u/gps_a^3)*gps_t_t0-2*pi*gps_k; %sqrt(u/gps_a^3)=[m/s]
    if gps_teM(t+1)==0
        gps_teM(t+1)=gps_teM(t+1)+2*pi;
    end
    
    gps_q=fix(gps_teM(t+1)/(2*pi)); 
    gps_M(t+1)=gps_teM(t+1)-2*pi*gps_q;

    if gps_M(t+1)==0
        gps_M(t+1)=gps_M(t+1)+2*pi;
    end
    
    i=1;
    gps_tE(1)=gps_E(t);
    while (1)
        gps_tE(i+1)=gps_M(t+1)+gps_e*sin(gps_tE(i));
        if abs(gps_tE(i+1)-gps_tE(i)) < gps_e
            break
        end
        i=i+1;
    end
    gps_E(t+1)=gps_tE(end);

end


%plot groung track with 'geoplot'
figure(1)
geoplot(gps_lat, gps_lon,'b.')
geolimits([-70 70],[-180 180]) %earth angular rate must contain
legend('GPS')

%plot sky view with 'skyplot'
figure(2)
skyplot(gps_az, gps_el)
legend('GPS')

%plot satellite orbit
stopTime_gps=startTime_gps+days(1);
sampleTime=60; %seconds
sc_gps = satelliteScenario(startTime_gps,stopTime_gps,sampleTime);

gps_sc_i=rad2deg(gps_i);
gps_sc_OMEGA=rad2deg(gps_OMEGA);
gps_sc_omega=rad2deg(gps_omega);
gps_sc_ta=rad2deg(gps_ta(1));

sat_gps = satellite(sc_gps,gps_a,gps_e,gps_sc_i,gps_sc_OMEGA,gps_sc_omega,gps_sc_ta);
gs_gps=groundStation(sc_gps,gs_lat,gs_lon);
ac_gps=access(sat_gps,gs_gps);

leadTime = 24*3600; % seconds
trailTime = leadTime;
gt_gps = groundTrack(sat_gps,"LeadTime",leadTime,"TrailTime",trailTime);

satelliteScenarioViewer(sc_gps)

    