function rangeInPQW = solve_rad_RangeInPerifocalFrame(semimajor_axis, eccentricity, true_anomaly) %input unit-> [rad] & [m]
p=semimajor_axis*(1-eccentricity^2);  %[m]
r=p/(1+eccentricity*cos(true_anomaly));  %[m]

rangeInPQW=[r*cos(true_anomaly); r*sin(true_anomaly); 0]; %[m]