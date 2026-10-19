
function dz = motorODE(t,z,u,params)

Jeq = params.Jeq;
B = params.B;
Km = params.Km;
Lm = params.Lm;
Rm = params.Rm;

omega = z(3:4);
i = z(5:6);

domega = (Km*i - B*omega)/Jeq;
di = (u(t,z) - Rm*i - Km*omega)/Lm;

dz = [omega; domega; di];

end
