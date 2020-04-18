function [tv, gam, gamd, gamdd, Fg, Fgd, Fgdd] = MakePolys(ts, dt, ang1, ang2)
% Make a polynomial as function of time for two angle path

% make path from ang1 to ang2
syms t a1 a2 a3 a4 a5 a6 a7
x = a1 + a2*t + a3*(t^2) + a4*(t^3) + a5*(t^4) + a6*(t^5) + a7*(t^6);
xd = diff(x, t);
xdd = diff(xd, t);

% assign boundary conditions
eq1 = subs(x, t, 0) - ang1 == 0;
eq2 = subs(x, t, dt) - ang2 == 0;
eq3 = subs(xd, t, 0) - 0 == 0;
eq4 = subs(xd, t, dt) - 0 == 0;
eq5 = subs(xdd, t, 0) - 0 == 0;
eq6 = subs(xdd, t, dt) - 0 == 0;
eq7 = subs(x, t, dt*0.5) - (ang2 - ang1)*0.5 + ang1 == 0;
% soln = vpasolve([eq1, eq2, eq3, eq4, eq5, eq6, eq7], ...
%     [a1 a2 a3 a4 a5 a6 a7]);
[b1, b2, b3, b4, b5, b6, b7] = vpasolve([eq1, eq2, eq3, eq4, eq5, eq6, eq7], ...
    [a1 a2 a3 a4 a5 a6 a7]);

% make first interpolation object
x1 = linspace(0, dt, 17);
y1 = double(b1 + b2*x1 + b3*(x1.^2) + b4*(x1.^3) + b5*(x1.^4) + b6*(x1.^5) + b7*(x1.^6));
y1d = double(6*b7*(x1.^5) + 5*b6*(x1.^4) + 4*b5*(x1.^3) + 3*b4*(x1.^2) + 2*b3*x1 + b2);
y1dd = double(30*b7*(x1.^4) + 20*b6*(x1.^3) + 12*b5*(x1.^2) + 6*b4*x1 + 2*b3);
% figure('color', 'w');
% hold on
% plot(x1, y1*180/pi, '-o')
% xlabel('Time [s]');
% ylabel('Angle [deg]');
% grid on
% keyboardn

% second boundary conditions
y2 = fliplr(y1);
y2d = fliplr(y1d);
y2dd = fliplr(y1dd);
x2 = x1(end) + x1;
tv = [x1, x2(2:end)] + ts;
gam = [y1, y2(2:end)];
gamd = [y1d, y2d(2:end)];
gamdd = [y1dd, y2dd(2:end)];
% figure('color', 'w');
% hold on
% plot(tv, gam*180/pi, '-o')
% plot(tv, gamd*180/pi, '-s')
% plot(tv, gamdd*180/pi, '-^')
% xlabel('Time [s]');
% ylabel('function');
% legend('Angle [deg]', 'Velocity [deg/s]', 'Acceleration [deg/s2]', ...
%     'location', 'southeast');
% grid on
% keyboard

% make interp options
tv = [0, tv, 10000];
gam = [ang1, gam, ang1];
gamd = [0, gamd, 0];
gamdd = [0, gamdd, 0];
Fg = griddedInterpolant(tv, gam);
Fgd = griddedInterpolant(tv, gamd);
Fgdd = griddedInterpolant(tv, gamdd);
% figure('color', 'w');
% tb = 0:0.1:10;
% hold on
% plot(tb, Fg(tb)*180/pi, '-o')
% plot(tb, Fgd(tb)*180/pi, '-s')
% plot(tb, Fgdd(tb)*180/pi, '-^')
% xlabel('Time [s]');
% ylabel('function');
% legend('Angle [deg]', 'Velocity [deg/s]', 'Acceleration [deg/s2]', ...
%     'location', 'southeast');
% grid on
% keyboard

end