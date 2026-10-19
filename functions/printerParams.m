
function params = printerParams()

rm = 0.006;
Jm = 5.8e-7;
m = 0.250;
Lm = 264e-6;
Rm = 7.14;
Km = 0.0286;

mu_r = 1;
mu_lx = 2;
mu_ly = 2;

Jeq = Jm + m*rm^2/2;

Br = 4*mu_r*rm^2*eye(2);

Bl = rm^2/4 * ...
    [2*mu_lx+2*mu_ly, 2*mu_lx-2*mu_ly;
     2*mu_lx-2*mu_ly, 2*mu_lx+2*mu_ly];

B = Bl + Br;

params = struct('rm',rm,'Jeq',Jeq,'B',B,'Km',Km,'Lm',Lm,'Rm',Rm);

end
