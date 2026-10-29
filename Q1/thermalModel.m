
clear; clc; close all

%% Parameters

rho = 1240;
cp = 1800;
k = 0.13;
w = 0.0004;
a = 0.08;
b = 0.04;

h = 10;
hb = 100;
Te = 22;
Tb = 60;
T0 = 200;

Tprint = 4;
dt = 0.01;
tmax = 180;
tol = 1;

%% Filament geometry

L = pi*(3*(a+b) - sqrt((3*a+b)*(a+3*b)));
N = 200;
ds = L/N;
s = (0:N-1)'*ds;
T = nan(N,1);
nPrinted = 0;

%% Heat equation

alpha = k/(rho*cp);
gamma1 = alpha/ds^2;
beta = (3*h + hb)/(rho*cp*w);
dt_max = 1/(2*gamma1 + beta);
assert(dt < dt_max, 'Time step too large')

%% Simulation 
% Finite volume method with forward Euler

time = [];
T_history = [];
temperatureDifference = [];

for t = 0:dt:tmax

    nActive = min(N, floor(t/Tprint*(N-1)) + 1);
    T(nPrinted+1:nActive) = T0;
    nPrinted = nActive;
    temp = T(1:nActive);
    Phi = zeros(nActive+1,1);

    if nActive > 1
        Phi(2:end-1) = ...
            (temp(1:end-1)-temp(2:end))*gamma1;
    end

    if nActive == N
        Phi(1) = (temp(end)-temp(1))*gamma1;
        Phi(end) = Phi(1);
    end

    f_T = (3*h*(Te-temp) + hb*(Tb-temp))/(rho*cp*w);
    
    T(1:nActive) = temp + dt*(...
        Phi(1:end-1) - Phi(2:end) + f_T);

    time(end+1) = t+dt;
    T_history(:,end+1) = T;

    temperatureDifference(end+1) = ...
        max(T(1:nActive)) - min(T(1:nActive));

    if t >= Tprint && temperatureDifference(end) < tol
        break
    end
end

%% Results

fprintf('Filament length: %.4f m\n', L)
fprintf('Time until uniform: %.2f s\n', time(end))
fprintf('Cooling after printing: %.2f s\n', time(end)-Tprint)
fprintf('Final temperature difference: %.3f C\n', ...
    temperatureDifference(end))

%% Temperature heatmap

figure

hImg = imagesc(s, time, T_history');
set(hImg, 'AlphaData', double(~isnan(T_history')))
set(gca, 'YDir', 'normal', 'Color', [0.8 0.8 0.8])
colormap(hot)
clim([Te T0])
colorbar
xlabel('Distance along filament (m)')
ylabel('Time (s)')
hold on
yline(Tprint, 'w--', 'Printing complete', 'LineWidth', 1.5)

%% Temperature uniformity

figure
plot(time, temperatureDifference, 'LineWidth', 1.5)
hold on
yline(tol,'--','Uniformity criterion')
xline(Tprint,'--','Printing complete')
xlabel('Time (s)')
ylabel('Maximum temperature difference (C)')
grid on
