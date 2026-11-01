
function [p,x0] = hybridParams()

params = printerParams();

Jeq = params.Jeq;
B11 = params.B(1,1);
B12 = params.B(1,2);
Km = params.Km;
Lm = params.Lm;
Rm = params.Rm;

Ln = 0.001;
Cn = 0.001;
Rn = 14.4;
Ct = 10;
Rt = 8;
Te = 22;

p = [Jeq; B11; B12; Km; Lm; Rm; ...
     Ln; Cn; Rn; Ct; Rt; Te];

xstart = 0.10;
ystart = 0.10;
rm = params.rm;

Tprint = 200;
Vhold = sqrt(Rn*(Tprint-Te)/Rt);

x0 = [(xstart+ystart)/rm;
      (xstart-ystart)/rm;
      zeros(4,1);
      Vhold/Rn;
      Vhold;
      Tprint]


end
