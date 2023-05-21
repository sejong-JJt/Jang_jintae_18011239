function C_pqw2eci = PQW2ECI(arg_prg, inc_angle, RAAN)
Rw=[cos(arg_prg) sin(arg_prg) 0; -sin(arg_prg) cos(arg_prg) 0; 0 0 1];
Ri=[1 0 0; 0 cos(inc_angle) sin(inc_angle); 0 -sin(inc_angle) cos(inc_angle)];
Ro=[cos(RAAN) sin(RAAN) 0; -sin(RAAN) cos(RAAN) 0; 0 0 1];

C_pqw2eci_beforetranspose=Ro*Ri*Rw;
C_pqw2eci=C_pqw2eci_beforetranspose.';
