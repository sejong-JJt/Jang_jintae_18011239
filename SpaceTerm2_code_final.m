clear
close
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

qzss_M=zeros(1,1440);
qzss_a=getfield(nav.QZSS,'a');           %Semi-major axis [m]
qzss_e=getfield(nav.QZSS,'e');           %Eccentricity
qzss_i=getfield(nav.QZSS,'i');           %Inclination [rad]
qzss_omega=getfield(nav.QZSS,'omega');   %Argument of perigee,w [rad]
qzss_M(1)=getfield(nav.QZSS,'M0');       %Mean anomaly at toc [rad]
qzss_toc=getfield(nav.QZSS,'toc');       %Epoch,current time
qzss_OMEGA=getfield(nav.QZSS,'OMEGA');   %RAAN [rad]

bds_M=zeros(1,1440);
bds_a=getfield(nav.BDS,'a');           %Semi-major axis [m]
bds_e=getfield(nav.BDS,'e');           %Eccentricity
bds_i=getfield(nav.BDS,'i');           %Inclination [rad]
bds_omega=getfield(nav.BDS,'omega');   %Argument of perigee,w [rad]
bds_M(1)=getfield(nav.BDS,'M0');       %Mean anomaly at toc [rad]
bds_toc=getfield(nav.BDS,'toc');       %Epoch,current time
bds_OMEGA=getfield(nav.BDS,'OMEGA');   %RAAN [rad]



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



startTm_qzss=qzss_toc;
startTime_qzss=datetime(qzss_toc);

qzss_UT=[2023, 5, 9, 2, 0, 0]; %Universal time on Mar 9, 2023
qzss_UT_dt=datetime(qzss_UT);
qzss_UT_jd=juliandate(qzss_UT_dt);
qzss_UT_GMST=siderealTime(qzss_UT_jd); 

qzss_ENU=[0 0 0];
qzss_p=qzss_a*(1-qzss_e^2); %[m]



startTm_bds=bds_toc;
startTime_bds=datetime(bds_toc);

bds_UT=[2023, 2, 26, 0, 0, 0]; %Universal time on Apr 29, 2023
bds_UT_dt=datetime(bds_UT);
bds_UT_jd=juliandate(bds_UT_dt);
bds_UT_GMST=siderealTime(bds_UT_jd); 

bds_ENU=[0 0 0];
bds_p=bds_a*(1-bds_e^2); %[m]

%ground station
gs_lat=37;      %[deg]
gs_lon=127;     %[deg]
gs_h=1000;      %[m]
el_mask=10;     %[deg]

%%
%GPS
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


%%
%QZSS
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



%%
%BDS
i=1;
bds_E(1)=bds_M(1);
while (1)
    bds_E(i+1)=bds_M(1)+bds_e*sin(bds_E(i));
    if abs(bds_E(i+1)-bds_E(i)) < bds_e
        break
    end
    i=i+1;
end

bds_E(1)=bds_E(end);

for t=1:1:1440 %24hour minutes
    %new GMST
    bds_DCM=ECI2ECEF_DCM_BDS(bds_toc);
    
    %find true_anomaly at t
    bds_y(t)=(sqrt(1-bds_e^2)*sin(bds_E(t)))/(1-bds_e*cos(bds_E(t)));
    bds_x(t)=(cos(bds_E(t))-bds_e)/(1-bds_e*cos(bds_E(t)));

    bds_ta(t)=mod(atan2(bds_y(t),bds_x(t)),2*pi); %true_anomaly [rad]
    
    %PQW range & velocity with true_anomaly
    bds_r(t)=bds_p/(1+bds_e*cos(bds_ta(t))); %[m]
    bds_rangeInPQW = solve_rad_RangeInPerifocalFrame(bds_a, bds_e, bds_ta(t));
    bds_velocityInPQW = solve_rad_VelocityInPerifocalFrame(bds_a, bds_e, bds_ta(t));
    
    %change PQW -> ECI
    bds_C_pqw2eci = PQW2ECI(bds_omega, bds_i, bds_OMEGA);
    bds_rangeInECI=bds_C_pqw2eci*bds_rangeInPQW;
    bds_velocityInECI=bds_C_pqw2eci*bds_velocityInPQW;
    
    %ECI -> ECEF
    bds_rangeInECEF=bds_DCM*bds_rangeInECI; %[m]
    bds_velocityInECEF=bds_DCM*bds_velocityInECI; %[m/s]
    
    %ECEF-> geodetic lat,long
    earth=wgs84Ellipsoid('meter');
    [bds_lat(t), bds_lon(t), bds_h(t)]=ecef2geodetic(earth, bds_rangeInECEF(1), bds_rangeInECEF(2), bds_rangeInECEF(3)); %geodetic lat,long [deg,deg,m]
    
    
    %ENU -> Az & El
    [bds_east(t), bds_north(t), bds_up(t)]=ecef2enu(bds_rangeInECEF(1), bds_rangeInECEF(2), bds_rangeInECEF(3),gs_lat,gs_lon,gs_h,earth);
    bds_ENU=[bds_east(t), bds_north(t), bds_up(t)];
     
    bds_az(t)=azimuth(bds_ENU); %[deg]
    bds_el(t)=elevation(bds_ENU, el_mask); %[deg]
    

    %time pass "1 minute", at t=1-> 00:00:00/ t=2-> 00:01:00
    bds_toc(1,5)=bds_toc(1,5)+1;
    if bds_toc(1,5)>=60
        bds_toc(1,5)=bds_toc(1,5)-60;
        bds_toc(1,4)=bds_toc(1,4)+1;
    end

    if bds_toc(1,4)>=24
        bds_toc(1,4)=bds_toc(1,4)-24;
        bds_toc(1,3)=bds_toc(1,3)+1;
    end

    %find Mean & Eccentricity anomaly at t+1
    bds_t_t0=etime(bds_toc,startTm_bds); %t0=startTime
    bds_tau=sqrt(bds_a^3/u)*2*pi;
    bds_k=fix(bds_t_t0/bds_tau); 
    bds_teM(t+1)=bds_M(1)+sqrt(u/bds_a^3)*bds_t_t0-2*pi*bds_k; %sqrt(u/gps_a^3)=[m/s]
    if bds_teM(t+1)==0
        bds_teM(t+1)=bds_teM(t+1)+2*pi;
    end
    
    bds_q=fix(bds_teM(t+1)/(2*pi)); 
    bds_M(t+1)=bds_teM(t+1)-2*pi*bds_q;

    if bds_M(t+1)==0
        bds_M(t+1)=bds_M(t+1)+2*pi;
    end
    
    i=1;
    bds_tE(1)=bds_E(t);
    while (1)
        bds_tE(i+1)=bds_M(t+1)+bds_e*sin(bds_tE(i));
        if abs(bds_tE(i+1)-bds_tE(i)) < bds_e*0.01  % higher accuracy
            break
        end
        i=i+1;
    end
    bds_E(t+1)=bds_tE(end);

end


%%plot groung track with 'geoplot'
figure(1)
geoplot(gps_lat, gps_lon,'b.')
hold on
geoplot(qzss_lat, qzss_lon,'r.')
hold on
geoplot(bds_lat, bds_lon,'y.')
legend('GPS','QZSS','BDS')
geolimits([-70 70],[-180 180])

%%plot sky view with 'skyplot'
figure(2)
az=[gps_az,qzss_az,bds_az];
el=[gps_el,qzss_el,bds_el];
g=[ones(size(gps_az)) ones(size(qzss_az))+1 ones(size(bds_az))+2];

skyplot(az,el,GroupData=categorical(g))
legend('GPS','QZSS','BDS')

%%plot satellite orbit
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



stopTime_qzss=startTime_qzss+days(1);
sc_qzss = satelliteScenario(startTime_qzss,stopTime_qzss,sampleTime);

qzss_sc_i=rad2deg(qzss_i);
qzss_sc_OMEGA=rad2deg(qzss_OMEGA);
qzss_sc_omega=rad2deg(qzss_omega);
qzss_sc_ta=rad2deg(qzss_ta(1));

sat_qzss = satellite(sc_qzss,qzss_a,qzss_e,qzss_sc_i,qzss_sc_OMEGA,qzss_sc_omega,qzss_sc_ta);
gs_qzss=groundStation(sc_qzss,gs_lat,gs_lon);
ac=access(sat_qzss,gs_qzss);

gt_qzss = groundTrack(sat_qzss,"LeadTime",leadTime,"TrailTime",trailTime);

satelliteScenarioViewer(sc_qzss)



stopTime_bds=startTime_bds+days(1);
sc_bds = satelliteScenario(startTime_bds,stopTime_bds,sampleTime);

bds_sc_i=rad2deg(bds_i);
bds_sc_OMEGA=rad2deg(bds_OMEGA);
bds_sc_omega=rad2deg(bds_omega);
bds_sc_ta=rad2deg(bds_ta(1));

sat_bds = satellite(sc_bds,bds_a,bds_e,bds_sc_i,bds_sc_OMEGA,bds_sc_omega,bds_sc_ta);
gs_bds=groundStation(sc_bds,gs_lat,gs_lon);
ac=access(sat_bds,gs_bds);

gt_bds = groundTrack(sat_bds,"LeadTime",leadTime,"TrailTime",trailTime);

satelliteScenarioViewer(sc_bds)

