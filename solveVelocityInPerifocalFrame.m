function velocityInPQW = solveVelocityInPerifocalFrame(semimajor_axis, eccentricity, true_anomaly) %input unit-> [rad] & [m]
p=semimajor_axis*(1-eccentricity^2); %[m]
u=3.986004418*10^5; %[km^3 s^-2]

velocityInPQW=sqrt(u/p)*[-sin(true_anomaly*(pi/180)); eccentricity+cos(true_anomaly*(pi/180)); 0]; %[km/s]