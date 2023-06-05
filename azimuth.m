function az = azimuth(ENU) %ENU has time parameter(t) -> [East(t), North(t), Up(t)] [km]
az_rad=acos(ENU(2)/sqrt(ENU(1)^2+ENU(2)^2)); %[rad]
if ENU(1)>0
    az_rad=az_rad;
else
    az_rad=2*pi-az_rad;
end                      %To express az_rad over pi

az=rad2deg(az_rad); %[deg]
