function el = elevation(ENU, el_mask) %ENU has time parameter(t) -> [East(t), North(t), Up(t)] [km]&[deg]
el_rad=asin(ENU(3)/sqrt(ENU(1)^2+ENU(2)^2+ENU(3)^2)); %[rad]
el=rad2deg(el_rad); %[deg]

if el>el_mask
    el=el;
else
    el=NaN;
end
