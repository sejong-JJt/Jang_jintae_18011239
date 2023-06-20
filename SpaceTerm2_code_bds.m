close
clear
clc

load nav.mat
bds_M=zeros(1,1440);
bds_a=getfield(nav.BDS,'a');           %Semi-major axis [m]
bds_e=getfield(nav.BDS,'e');           %Eccentricity
bds_i=getfield(nav.BDS,'i');           %Inclination [rad]
bds_omega=getfield(nav.BDS,'omega');   %Argument of perigee,w [rad]
bds_M(1)=getfield(nav.BDS,'M0');       %Mean anomaly at toc [rad]
bds_toc=getfield(nav.BDS,'toc');       %Epoch,current time
bds_OMEGA=getfield(nav.BDS,'OMEGA');   %RAAN [rad]

startTm_bds=bds_toc;
startTime_bds=datetime(bds_toc);
w_earth=360/86160; %earth angular rate [deg/s]

bds_UT=[2023, 2, 26, 0, 0, 0]; %Universal time on Feb 26, 2023
bds_UT_dt=datetime(bds_UT);
bds_UT_jd=juliandate(bds_UT_dt);
bds_UT_GMST=siderealTime(bds_UT_jd); 

bds_ENU=[0 0 0];

u=3.986004418*10^14; %[m^3 s^-2]

bds_p=bds_a*(1-bds_e^2); %[m]

%ground station
gs_lat=37;      %[deg]
gs_lon=127;     %[deg]
gs_h=1000;      %[m]
el_mask=10;     %[deg]

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


%plot groung track with 'geoplot'
figure(1)
geoplot(bds_lat, bds_lon,'y.')
geolimits([-70 70],[-180 180]) %earth angular rate must contain
legend('BDS')

%plot sky view with 'skyplot'
figure(2)
skyplot(bds_az, bds_el)
legend('BDS')

%plot satellite orbit
stopTime_bds=startTime_bds+days(1);
sampleTime=60; %seconds
sc_bds = satelliteScenario(startTime_bds,stopTime_bds,sampleTime);

bds_sc_i=rad2deg(bds_i);
bds_sc_OMEGA=rad2deg(bds_OMEGA);
bds_sc_omega=rad2deg(bds_omega);
bds_sc_ta=rad2deg(bds_ta(1));

sat_bds = satellite(sc_bds,bds_a,bds_e,bds_sc_i,bds_sc_OMEGA,bds_sc_omega,bds_sc_ta);
gs_bds=groundStation(sc_bds,gs_lat,gs_lon);
ac=access(sat_bds,gs_bds);

leadTime = 24*3600; % seconds
trailTime = leadTime;
gt_bds = groundTrack(sat_bds,"LeadTime",leadTime,"TrailTime",trailTime);

satelliteScenarioViewer(sc_bds)
