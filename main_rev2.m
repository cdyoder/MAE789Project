% yoder code to run TRREx simulations
clear
close all
clc
fclose('all');

%{
To do;
1. Need to 
%}

% derive eoms
statesname = 'TRREx_SimFile_rev2';           % file to make 
% need to change this down at the ode45 line
% keyboard

% prevent overwrite
if exist([statesname, '.m'], "file") ~= 2
    tic
    statesfile = DeriveEOMS(statesname);
    t5 = toc
end

% variable schedule
% ----------------------
% x(1) = thB (rad)      % x(5) = gam2  (rad)     % x(9) = gam4   (rad)
% x(2) = dthB (rad/s)   % x(6) = dgam2 (rad/s)   % x(10) = dgam4 (rad/s)
% x(3) = gam1 (rad)     % x(7) = gam3  (rad)     
% x(4) = dgam1 (rad/s)  % x(8) = dgam3 (rad/s)

% ics
xics = zeros(10, 1);

% time
tv = 0:0.05:20;

% parameters
Crr_nom = 0.07;
th_trig = 1e-3;

% gamma schedule
% make interpolation objects to get gamma as a function of time
Gam1dd = griddedInterpolant([0, 10000], [0, 0]);
Gam2dd = Gam1dd;
[~, ~, ~, ~, ~, ~, Gam4dd] = MakePolys(3, 3, 0, 30*pi/180);
[~, ~, ~, ~, ~, ~, Gam3dd] = MakePolys(4, 3, 0, 30*pi/180);

% run sim
disp('Simulation started');
optz = odeset('Stats', 'on');
tic
% [~, outsB] = ode45(@(tt, xx)statesname(tt, xx, Crr_nom, th_trig, ...
%     Gam1dd, Gam2dd, Gam3dd, Gam4dd, 0), tv, xics);
[~, outsB] = ode45(@(tt, xx)TRREx_SimFile_rev2(tt, xx, Crr_nom, th_trig, ...
    Gam1dd, Gam2dd, Gam3dd, Gam4dd,0), tv, xics, optz);
t6 = toc
disp('Simulation done');

% get quantities afterwards
Ff = NaN(length(tv), 1);
Frr = Ff;
Fn = Ff;
for i1 = 1:length(tv)
    xics = outsB(i1, :);
    xd = TRREx_SimFile_rev2(tv(i1), xics, Crr_nom, th_trig, ...
        Gam1dd, Gam2dd, Gam3dd, Gam4dd, 1);
    Ff(i1) = xd.Ff;
    Fn(i1) = xd.Fn;
    Frr(i1) = xd.Frr;
end

% plot
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
figure('color', 'w');
hold on
yyaxis left
plot(tv, outsB(:, 1)*180/pi);
ylabel('Angle [deg]', 'interpreter', 'latex');
yyaxis right 
plot(tv, outsB(:, 2)*180/pi);
xlabel('Time [s]', 'interpreter', 'latex');
ylabel('Angle rate [deg/s]', 'interpreter', 'latex');
grid on

figure('color', 'w');
hold on
plot(tv, outsB(:, 3)*180/pi);
plot(tv, outsB(:, 5)*180/pi);
plot(tv, outsB(:, 7)*180/pi);
plot(tv, outsB(:, 9)*180/pi);
xlabel('Time [s]', 'interpreter', 'latex');
ylabel('Angle [deg]', 'interpreter', 'latex');
grid on
legend('$\gamma_1$', '$\gamma_2$', '$\gamma_3$', '$\gamma_4$', 'interpreter', 'latex');