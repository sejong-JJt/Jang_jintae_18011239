function velocityInPQW = solve_rad_VelocityInPerifocalFrame(semimajor_axis, eccentricity, true_anomaly) %input unit-> [rad] & [m]
p=semimajor_axis*(1-eccentricity^2); %[m]
u=3.986004418*10^14; %[m^3 s^-2]

velocityInPQW=sqrt(u/p)*[-sin(true_anomaly); eccentricity+cos(true_anomaly); 0]; %[m/s]