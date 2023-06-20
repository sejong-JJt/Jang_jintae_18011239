close
clear
clc

load nav.mat
qzss_M=zeros(1,1440);
qzss_a=getfield(nav.QZSS,'a');           %Semi-major axis [m]
qzss_e=getfield(nav.QZSS,'e');           %Eccentricity
qzss_i=getfield(nav.QZSS,'i');           %Inclination [rad]
qzss_omega=getfield(nav.QZSS,'omega');   %Argument of perigee,w [rad]
qzss_M(1)=getfield(nav.QZSS,'M0');       %Mean anomaly at toc [rad]
qzss_toc=getfield(nav.QZSS,'toc');       %Epoch,current time
qzss_OMEGA=getfield(nav.QZSS,'OMEGA');   %RAAN [rad]

startTm_qzss=qzss_toc;
startTime_qzss=datetime(qzss_toc);
w_earth=360/86160; %earth angular rate [deg/s]

qzss_UT=[2023, 5, 9, 2, 0, 0]; %Universal time on Apr 29, 2023
qzss_UT_dt=datetime(qzss_UT);
qzss_UT_jd=juliandate(qzss_UT_dt);
qzss_UT_GMST=siderealTime(qzss_UT_jd); 

qzss_ENU=[0 0 0];

u=3.986004418*10^14; %[m^3 s^-2]

qzss_p=qzss_a*(1-qzss_e^2); %[m]

%ground station
gs_lat=37;      %[deg]
gs_lon=127;     %[deg]
gs_h=1000;      %[m]
el_mask=10;     %[deg]

i=1;
qzss_E(1)=qzss_M(1);
while (1)
    qzss_E(i+1)=qzss_M(1)+qzss_e*sin(qzss_E(i));
    if abs(qzss_E(i+1)-qzss_E(i)) < qzss_e
        break
    end
    i=i+1;
end

qzss_E(1)=qzss_E(end);

for t=1:1:1440 %24hour minutes
    %new GMST
    qzss_DCM=ECI2ECEF_DCM_QZSS(qzss_toc);
    
    %find true_anomaly at t
    qzss_y(t)=(sqrt(1-qzss_e^2)*sin(qzss_E(t)))/(1-qzss_e*cos(qzss_E(t)));
    qzss_x(t)=(cos(qzss_E(t))-qzss_e)/(1-qzss_e*cos(qzss_E(t)));

    qzss_ta(t)=mod(atan2(qzss_y(t),qzss_x(t)),2*pi); %true_anomaly [rad]
    
    %PQW range & velocity with true_anomaly
    qzss_r(t)=qzss_p/(1+qzss_e*cos(qzss_ta(t))); %[m]
    qzss_rangeInPQW = solve_rad_RangeInPerifocalFrame(qzss_a, qzss_e, qzss_ta(t));
    qzss_velocityInPQW = solve_rad_VelocityInPerifocalFrame(qzss_a, qzss_e, qzss_ta(t));
    
    %change PQW -> ECI
    qzss_C_pqw2eci = PQW2ECI(qzss_omega, qzss_i, qzss_OMEGA);
    qzss_rangeInECI=qzss_C_pqw2eci*qzss_rangeInPQW;
    qzss_velocityInECI=qzss_C_pqw2eci*qzss_velocityInPQW;
    
    %ECI -> ECEF
    qzss_rangeInECEF=qzss_DCM*qzss_rangeInECI; %[m]
    qzss_velocityInECEF=qzss_DCM*qzss_velocityInECI; %[m/s]
    
    %ECEF-> geodetic lat,long
    earth=wgs84Ellipsoid('meter');
    [qzss_lat(t), qzss_lon(t), qzss_h(t)]=ecef2geodetic(earth, qzss_rangeInECEF(1), qzss_rangeInECEF(2), qzss_rangeInECEF(3)); %geodetic lat,long [deg,deg,m]
    
    
    %ENU -> Az & El
    [qzss_east(t), qzss_north(t), qzss_up(t)]=ecef2enu(qzss_rangeInECEF(1), qzss_rangeInECEF(2), qzss_rangeInECEF(3),gs_lat,gs_lon,gs_h,earth);
    qzss_ENU=[qzss_east(t), qzss_north(t), qzss_up(t)];
     
    qzss_az(t)=azimuth(qzss_ENU); %[deg]
    qzss_el(t)=elevation(qzss_ENU, el_mask); %[deg]
    

    %time pass "1 minute", at t=1-> 00:00:00/ t=2-> 00:01:00
    qzss_toc(1,5)=qzss_toc(1,5)+1;
    if qzss_toc(1,5)>=60
        qzss_toc(1,5)=qzss_toc(1,5)-60;
        qzss_toc(1,4)=qzss_toc(1,4)+1;
    end

    if qzss_toc(1,4)>=24
        qzss_toc(1,4)=qzss_toc(1,4)-24;
        qzss_toc(1,3)=qzss_toc(1,3)+1;
    end

    %find Mean & Eccentricity anomaly at t+1
    qzss_t_t0=etime(qzss_toc,startTm_qzss); %t0=startTime
    qzss_tau=sqrt(qzss_a^3/u)*2*pi;
    qzss_k=fix(qzss_t_t0/qzss_tau); 
    qzss_teM(t+1)=qzss_M(1)+sqrt(u/qzss_a^3)*qzss_t_t0-2*pi*qzss_k; %sqrt(u/gps_a^3)=[m/s]
    if qzss_teM(t+1)==0
        qzss_teM(t+1)=qzss_teM(t+1)+2*pi;
    end
    
    qzss_q=fix(qzss_teM(t+1)/(2*pi)); 
    qzss_M(t+1)=qzss_teM(t+1)-2*pi*qzss_q;

    if qzss_M(t+1)==0
        qzss_M(t+1)=qzss_M(t+1)+2*pi;
    end
    
    i=1;
    qzss_tE(1)=qzss_E(t);
    while (1)
        qzss_tE(i+1)=qzss_M(t+1)+qzss_e*sin(qzss_tE(i));
        if abs(qzss_tE(i+1)-qzss_tE(i)) < qzss_e*0.01 % higher accuracy
            break
        end
        i=i+1;
    end
    qzss_E(t+1)=qzss_tE(end);

end


%plot groung track with 'geoplot'
figure(1)
geoplot(qzss_lat, qzss_lon,'r.')
geolimits([-70 70],[-180 180]) %earth angular rate must contain
legend('QZSS')

%plot sky view with 'skyplot'
figure(2)
skyplot(qzss_az, qzss_el)
legend('QZSS')

%plot satellite orbit
stopTime_qzss=startTime_qzss+days(1);
sampleTime=60; %seconds
sc_qzss = satelliteScenario(startTime_qzss,stopTime_qzss,sampleTime);

qzss_sc_i=rad2deg(qzss_i);
qzss_sc_OMEGA=rad2deg(qzss_OMEGA);
qzss_sc_omega=rad2deg(qzss_omega);
qzss_sc_ta=rad2deg(qzss_ta(1));

sat_qzss = satellite(sc_qzss,qzss_a,qzss_e,qzss_sc_i,qzss_sc_OMEGA,qzss_sc_omega,qzss_sc_ta);
gs_qzss=groundStation(sc_qzss,gs_lat,gs_lon);
ac=access(sat_qzss,gs_qzss);

leadTime = 24*3600; % seconds
trailTime = leadTime;
gt_qzss = groundTrack(sat_qzss,"LeadTime",leadTime,"TrailTime",trailTime);

satelliteScenarioViewer(sc_qzss)
