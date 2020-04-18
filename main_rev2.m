% yoder code to run TRREx simulations
clear
close all
clc
fclose('all');

%{
To do;
1. Need to fix math error in the code
2. Need to change the differene in gamma angles for trrex
%}

% derive eoms
statesname = 'TRREx_SimFile_rev2';           % file to make 
% need to change this down at the ode45 line
% keyboard

% prevent overwrite
if exist([statesname, '.m'], "file") ~= 2
    tic
    disp('States file not found, thus starting derivation');
    statesfile = DeriveEOMS(statesname);
    t5 = toc
    disp('Derivation complete');
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
tv = 0:0.05:1;

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
[ta, outsB] = ode45(@(tt, xx)TRREx_SimFile_rev2(tt, xx, Crr_nom, th_trig, ...
    Gam1dd, Gam2dd, Gam3dd, Gam4dd,0), tv, xics, optz);
t6 = toc
disp('Simulation done');

% get quantities afterwards
Ff = NaN(length(ta), 1);
Frr = Ff;
Fn = Ff;
for i1 = 1:length(ta)
    xics = outsB(i1, :);
    xd = TRREx_SimFile_rev2(ta(i1), xics, Crr_nom, th_trig, ...
        Gam1dd, Gam2dd, Gam3dd, Gam4dd, 1);
    Ff(i1) = xd.Ffr;
    Fn(i1) = xd.Fn;
    Frr(i1) = xd.Frr;
end




% plot setup
figdir = 'bin';
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultFigureUnits', 'inches');
pp = [0, 0, 3, 2.5];
fs = 8;


% plot of angle, angle dot
figure('color', 'w');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', pp(3:4));
set(gcf, 'PaperPosition', pp);
set(gcf, 'Position', [3, 3, pp(3), pp(4)]);
hold on
yyaxis left
plot(ta, outsB(:, 1)*180/pi);
ylabel('$\theta$ [deg]', 'interpreter', 'latex');
yyaxis right 
plot(ta, outsB(:, 2)*180/pi);
xlabel('Time [s]', 'interpreter', 'latex');
ylabel('$\dot{\theta}$ [deg/s]', 'interpreter', 'latex');
grid on
set(gca, 'FontSize', fs);
figname = 'angles';
savefig(gcf, fullfile(figdir, [figname, '.fig']));
print(fullfile(figdir, figname), '-dpdf');
print(fullfile(figdir, figname), '-dpng');


% plot of x, xdot
rCH = 0.3937;
figure('color', 'w');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', pp(3:4));
set(gcf, 'PaperPosition', pp);
set(gcf, 'Position', [3, 3, pp(3), pp(4)]);

hold on
yyaxis left
plot(ta, outsB(:, 1)*rCH);
ylabel('$x$ [m]', 'interpreter', 'latex');
yyaxis right 
plot(ta, outsB(:, 2)*rCH);
xlabel('Time [s]', 'interpreter', 'latex');
ylabel('$\dot{x}$ [m/s]', 'interpreter', 'latex');
grid on
set(gca, 'FontSize', fs);
figname = 'posn';
savefig(gcf, fullfile(figdir, [figname, '.fig']));
print(fullfile(figdir, figname), '-dpdf');
print(fullfile(figdir, figname), '-dpng');


% plot of arm angles
figure('color', 'w');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', pp(3:4));
set(gcf, 'PaperPosition', pp);
set(gcf, 'Position', [3, 3, pp(3), pp(4)]);

hold on
plot(ta, outsB(:, 3)*180/pi);
plot(ta, outsB(:, 5)*180/pi);
plot(ta, outsB(:, 7)*180/pi);
plot(ta, outsB(:, 9)*180/pi);
xlabel('Time [s]', 'interpreter', 'latex');
ylabel('Angle [deg]', 'interpreter', 'latex');
grid on
legend('$\gamma_1$', '$\gamma_2$', '$\gamma_3$', '$\gamma_4$', 'interpreter', 'latex');
figname = 'arms';
savefig(gcf, fullfile(figdir, [figname, '.fig']));
print(fullfile(figdir, figname), '-dpdf');
print(fullfile(figdir, figname), '-dpng');


% plot of forces
figure('color', 'w');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', pp(3:4));
set(gcf, 'PaperPosition', pp);
set(gcf, 'Position', [3, 3, pp(3), pp(4)]);

hold on
plot(ta, Ff);
plot(ta, Fn);
plot(ta, Frr);
xlabel('Time [s]', 'interpreter', 'latex');
ylabel('Force [N]', 'interpreter', 'latex');
grid on
legend('$F_f$', '$F_N$', '$F_{rr}$', 'interpreter', 'latex');
figname = 'forces';
savefig(gcf, fullfile(figdir, [figname, '.fig']));
print(fullfile(figdir, figname), '-dpdf');
print(fullfile(figdir, figname), '-dpng');