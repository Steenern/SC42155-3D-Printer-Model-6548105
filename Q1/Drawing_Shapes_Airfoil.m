clear; clc; close all
scriptFolder = fileparts(mfilename('fullpath'));
functionsFolder = fullfile(scriptFolder, '..', 'functions');
addpath(functionsFolder)

% Desired airfoil path

c = 0.16;
tau = 0.16;
xc = 0.02;
yc = 0.10;
bed = 0.20;
T = 4;

n = 400;
x_points = linspace(0,c,n/2);
xr = x_points/c;

yt = 5*tau*c*(0.2969*sqrt(xr) - 0.1260*xr ...
    - 0.3516*xr.^2 + 0.2843*xr.^3 - 0.1036*xr.^4);

x_des = xc + [x_points, fliplr(x_points)];
y_des = yc + [yt, -fliplr(yt)];

% Approach trajectory

T_approach = 2;

t = linspace(0,T_approach+T,3000);

x = zeros(size(t));
y = zeros(size(t));

idx1 = t <= T_approach;
idx2 = t > T_approach;

s1 = t(idx1)/T_approach;
p1 = 10*s1.^3 - 15*s1.^4 + 6*s1.^5;

x(idx1) = xc*p1;
y(idx1) = yc*p1;

distance = [0 cumsum(hypot(diff(x_des),diff(y_des)))];
distance = distance/distance(end);

[distance_unique,idx_unique] = unique(distance);

s2 = (t(idx2)-T_approach)/T;
p2 = 10*s2.^3 - 15*s2.^4 + 6*s2.^5;

x(idx2) = interp1(distance_unique, ...
    x_des(idx_unique),p2,'pchip');

y(idx2) = interp1(distance_unique, ...
    y_des(idx_unique),p2,'pchip');

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

fprintf('RMS airfoil error: %.6f mm\n',rms_error)
fprintf('Maximum airfoil error: %.6f mm\n',max_error)


% Plots

figure
hold on

rectangle('Position',[0 0 bed bed])
plot(x_des,y_des,':','LineWidth',2.5)
plot(x_sim,y_sim,'LineWidth',1.5)

[~,k] = max(ee);
t_max = te(k);

x_max = interp1(t_sim,x_sim,t_max);
y_max = interp1(t_sim,y_sim,t_max);

plot(x_max,y_max,'rx','MarkerSize',9,'LineWidth',1.5)
plot(0,0,'ko','MarkerFaceColor','k')

xlabel('x [m]')
ylabel('y [m]')
legend('Desired','Printer Path','Max Error','Origin', ...
    'Location','northeast','FontSize',12)
axis equal
xlim([0 bed])
ylim([0 bed])
grid on


figure
hold on

plot(te-T_approach,1000*ee,'LineWidth',1.5)
plot(t_max-T_approach,1000*ee(k),'rx', ...
    'MarkerSize',9,'LineWidth',1.5)

xlabel('Time [s]')
ylabel('Tracking error [mm]')
legend('Tracking Error','Max Error', ...
    'Location','northeast','FontSize',12)
set(gca,'YScale','log')
xlim([0 T])
grid on