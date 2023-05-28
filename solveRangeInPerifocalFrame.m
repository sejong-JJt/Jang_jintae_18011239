function rangeInPQW = solveRangeInPerifocalFrame(semimajor_axis, eccentricity, true_anomaly) %input unit-> [deg] & [km]
p=semimajor_axis*(1-eccentricity^2);  %[km]
r=p/(1+eccentricity*cos(true_anomaly*(pi/180)));  %[km]

rangeInPQW=[r*cos(true_anomaly*(pi/180)); r*sin(true_anomaly*(pi/180)); 0]; %[km]