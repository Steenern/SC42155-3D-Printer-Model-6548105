
function [Uff,theta,omega,i_des] = feedforward(t,x,y,params)

rm = params.rm;
Jeq = params.Jeq;
B = params.B;
Km = params.Km;
Lm = params.Lm;
Rm = params.Rm;

dx = gradient(x,t);
dy = gradient(y,t);

ddx = gradient(dx,t);
ddy = gradient(dy,t);

dddx = gradient(ddx,t);
dddy = gradient(ddy,t);

theta = [(x+y)/rm; (x-y)/rm];
omega = [(dx+dy)/rm; (dx-dy)/rm];
alpha = [(ddx+ddy)/rm; (ddx-ddy)/rm];
jerk = [(dddx+dddy)/rm; (dddx-dddy)/rm];

i_des = (Jeq*alpha+B*omega)/Km;
di_des = (Jeq*jerk+B*alpha)/Km;

Uff = Lm*di_des+Rm*i_des+Km*omega;

end
