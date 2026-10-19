clear; clc; close all
scriptFolder = fileparts(mfilename('fullpath'));
functionsFolder = fullfile(scriptFolder, '..', 'functions');
addpath(functionsFolder)

% Desired elliptical path

a = 0.04;
b = 0.02;
xc = 0.10;
yc = 0.10;
bed = 0.20;
T = 4;

w = 2*pi/T;
t = linspace(0,T,1000);

x_des = xc - a*cos(w*t);
y_des = yc + b*sin(w*t);


% Approach trajectory

T_approach = 3;

t = linspace(0,T_approach+T,3501);

x = zeros(size(t));
y = zeros(size(t));

idx1 = t <= T_approach;
idx2 = t > T_approach;

s1 = t(idx1)/T_approach;
p1 = 10*s1.^3 - 15*s1.^4 + 6*s1.^5;

x(idx1) = (xc-a)*p1;
y(idx1) = yc*p1;

s2 = (t(idx2)-T_approach)/T;
p2 = 10*s2.^3 - 15*s2.^4 + 6*s2.^5;

phi = 2*pi*p2;

x(idx2) = xc-a*cos(phi);
y(idx2) = yc+b*sin(phi);


% Feedforward printer head path

params = printerParams();
rm = params.rm;

[Uff,theta,omega,i_des] = feedforward(t,x,y,params);

u_ff = @(tt,z) interp1(t,Uff.',tt,'pchip').';

z0 = [theta(:,1); omega(:,1); i_des(:,1)];

[t_sim,z] = ode45(@(tt,z) motorODE(tt,z,u_ff,params), ...
    [0 T_approach+T],z0);

x_sim = rm/2*(z(:,1)+z(:,2));
y_sim = rm/2*(z(:,1)-z(:,2));

% Path tracking error

x_ref = interp1(t,x,t_sim,'pchip');
y_ref = interp1(t,y,t_sim,'pchip');

path_error = hypot(x_sim-x_ref,y_sim-y_ref);

idx = t_sim >= T_approach;
te = t_sim(idx);
ee = path_error(idx);

rms_error = 1000*sqrt(trapz(te,ee.^2)/T);
max_error = 1000*max(ee);

fprintf('RMS ellipse error: %.6f mm\n',rms_error)
fprintf('Maximum ellipse error: %.6f mm\n',max_error)


% Plot the paths

figure
hold on

rectangle('Position',[0 0 bed bed])
plot(x(idx2),y(idx2),':','LineWidth',1.5)
plot(x_sim,y_sim,'LineWidth',1.5)
plot(0,0,'ko','MarkerFaceColor','k')

xlabel('x [m]')
ylabel('y [m]')
legend('Desired','Feedforward','Origin','Location','best')
axis equal
xlim([0 bed])
ylim([0 bed])
grid on
