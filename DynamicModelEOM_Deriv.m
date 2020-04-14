%% Clear the workspace
clear all;

%% Some basic Inputs
mFileName = 'TRREx_SimFile.m';

%% symbolic variables go here
syms thetaB(t) betaG gamma1(t) gamma2(t) gamma3(t) gamma4(t) ...
    ...
    M mL ...    % masses
    ...
    g ...       % gravitational constant
    ...
    rCH...      % radius of chassis
    ...
    rBO_x_B(t) rBO_y_B(t) rBO_z_B(t) ... % Vector from O to B
    rh1B_x_B rh1B_y_B rh1B_z_B ...
    rh2B_x_B rh2B_y_B rh2B_z_B ...
    rh3B_x_B rh3B_y_B rh3B_z_B ...
    rh4B_x_B rh4B_y_B rh4B_z_B...     % vectors from B to each hinge point
    rC1h1_x_C1 rC1h1_y_C1 rC1h1_z_C1 ...
    rC2h2_x_C2 rC2h2_y_C2 rC2h2_z_C2 ...
    rC3h3_x_C3 rC3h3_y_C3 rC3h3_z_C3 ...
    rC4h4_x_C4 rC4h4_y_C4 rC4h4_z_C4 ... % vectors from each hinge to each leg CM
    ...
    I_Ch_xx I_Ch_xy I_Ch_xz ...
    I_Ch_yx I_Ch_yy I_Ch_yz ...
    I_Ch_zx I_Ch_zy I_Ch_zz ...
    ...
    I_L_xx I_L_xy I_L_xz ...
    I_L_yx I_L_yy I_L_yz ...
    I_L_zx I_L_zy I_L_zz ...
    ...
    Fn Ffr Frr Fd...        % force terms which can be calculated in a different
                            % function

%% Rotation Matrices

OcB = [
    cos(thetaB(t)) -sin(thetaB(t)) 0
    sin(thetaB(t)) cos(thetaB(t)) 0
    0 0 1
    ];

OcG = [
    cos(-betaG) -sin(-betaG) 0
    sin(-betaG) cos(-betaG) 0
    0 0 1
    ];

BcC1 = [
    cos(gamma1(t)) -sin(gamma1(t)) 0
    sin(gamma1(t)) cos(gamma1(t)) 0
    0 0 1
    ];

BcC2 = [
    cos(gamma2(t)) -sin(gamma2(t)) 0
    sin(gamma2(t)) cos(gamma2(t)) 0
    0 0 1
    ];

BcC3 = [
    cos(gamma3(t)) -sin(gamma3(t)) 0
    sin(gamma3(t)) cos(gamma3(t)) 0
    0 0 1
    ];

BcC4 = [
    cos(gamma4(t)) -sin(gamma4(t)) 0
    sin(gamma4(t)) cos(gamma4(t)) 0
    0 0 1
    ];

OcC1 = simplify( OcB * BcC1 );
OcC2 = simplify( OcB * BcC2 );
OcC3 = simplify( OcB * BcC3 );
OcC4 = simplify( OcB * BcC4 );

BcO = transpose(OcB);
GcO = transpose(OcG);
C1cB = transpose(BcC1);
C2cB = transpose(BcC2);
C3cB = transpose(BcC3);
C4cB = transpose(BcC4);
C1cO = transpose(OcC1);
C2cO = transpose(OcC2);
C3cO = transpose(OcC3);
C4cO = transpose(OcC4);

BcG = simplify(BcO * OcG);
GcB = transpose(BcG);

%% Angular Rates
OcB_dot = simplify(diff(OcB,t));
BcC1_dot = simplify(diff(BcC1,t));
BcC2_dot = simplify(diff(BcC2,t));
BcC3_dot = simplify(diff(BcC3,t));
BcC4_dot = simplify(diff(BcC4,t));
OcC1_dot = simplify(diff(OcC1,t));
OcC2_dot = simplify(diff(OcC2,t));
OcC3_dot = simplify(diff(OcC3,t));
OcC4_dot = simplify(diff(OcC4,t));

O_omega_B_x_B = [0 0 1]*BcO*OcB_dot*[0;1;0];
O_omega_B_y_B = [1 0 0]*BcO*OcB_dot*[0;0;1];
O_omega_B_z_B = [0 1 0]*BcO*OcB_dot*[1;0;0];

B_omega_C1_x_C1 = [0 0 1]*BcC1*BcC1_dot*[0;1;0];
B_omega_C1_y_C1 = [1 0 0]*BcC1*BcC1_dot*[0;0;1];
B_omega_C1_z_C1 = [0 1 0]*BcC1*BcC1_dot*[1;0;0];

B_omega_C2_x_C2 = [0 0 1]*BcC2*BcC2_dot*[0;1;0];
B_omega_C2_y_C2 = [1 0 0]*BcC2*BcC2_dot*[0;0;1];
B_omega_C2_z_C2 = [0 1 0]*BcC2*BcC2_dot*[1;0;0];

B_omega_C3_x_C3 = [0 0 1]*BcC3*BcC3_dot*[0;1;0];
B_omega_C3_y_C3 = [1 0 0]*BcC3*BcC3_dot*[0;0;1];
B_omega_C3_z_C3 = [0 1 0]*BcC3*BcC3_dot*[1;0;0];

B_omega_C4_x_C4 = [0 0 1]*BcC4*BcC4_dot*[0;1;0];
B_omega_C4_y_C4 = [1 0 0]*BcC4*BcC4_dot*[0;0;1];
B_omega_C4_z_C4 = [0 1 0]*BcC4*BcC4_dot*[1;0;0];

O_omega_C1_x_C1 = [0 0 1]*OcC1*OcC1_dot*[0;1;0];
O_omega_C1_y_C1 = [1 0 0]*OcC1*OcC1_dot*[0;0;1];
O_omega_C1_z_C1 = [0 1 0]*OcC1*OcC1_dot*[1;0;0];

O_omega_C2_x_C2 = [0 0 1]*OcC2*OcC2_dot*[0;1;0];
O_omega_C2_y_C2 = [1 0 0]*OcC2*OcC2_dot*[0;0;1];
O_omega_C2_z_C2 = [0 1 0]*OcC2*OcC2_dot*[1;0;0];

O_omega_C3_x_C3 = [0 0 1]*OcC3*OcC3_dot*[0;1;0];
O_omega_C3_y_C3 = [1 0 0]*OcC3*OcC3_dot*[0;0;1];
O_omega_C3_z_C3 = [0 1 0]*OcC3*OcC3_dot*[1;0;0];

O_omega_C4_x_C4 = [0 0 1]*OcC4*OcC4_dot*[0;1;0];
O_omega_C4_y_C4 = [1 0 0]*OcC4*OcC4_dot*[0;0;1];
O_omega_C4_z_C4 = [0 1 0]*OcC4*OcC4_dot*[1;0;0];

O_omega_B_B = simplify(transpose([O_omega_B_x_B O_omega_B_y_B O_omega_B_z_B]));
O_alpha_B = simplify(diff(O_omega_B_B));

B_omega_C1_C1 = simplify(transpose([B_omega_C1_x_C1 B_omega_C1_y_C1 B_omega_C1_z_C1]));
B_omega_C2_C2 = simplify(transpose([B_omega_C2_x_C2 B_omega_C2_y_C2 B_omega_C2_z_C2]));
B_omega_C3_C3 = simplify(transpose([B_omega_C3_x_C3 B_omega_C3_y_C3 B_omega_C3_z_C3]));
B_omega_C4_C4 = simplify(transpose([B_omega_C4_x_C4 B_omega_C4_y_C4 B_omega_C4_z_C4]));
B_alpha_C1 = simplify(diff(B_omega_C1_C1));
B_alpha_C2 = simplify(diff(B_omega_C2_C2));
B_alpha_C3 = simplify(diff(B_omega_C3_C3));
B_alpha_C4 = simplify(diff(B_omega_C4_C4));

O_omega_C1_C1 = simplify(transpose([O_omega_C1_x_C1 O_omega_C1_y_C1 O_omega_C1_z_C1]));
O_omega_C2_C2 = simplify(transpose([O_omega_C2_x_C2 O_omega_C2_y_C2 O_omega_C2_z_C2]));
O_omega_C3_C3 = simplify(transpose([O_omega_C3_x_C3 O_omega_C3_y_C3 O_omega_C3_z_C3]));
O_omega_C4_C4 = simplify(transpose([O_omega_C4_x_C4 O_omega_C4_y_C4 O_omega_C4_z_C4]));
O_alpha_C1 = simplify(diff(O_omega_C1_C1));
O_alpha_C2 = simplify(diff(O_omega_C2_C2));
O_alpha_C3 = simplify(diff(O_omega_C3_C3));
O_alpha_C4 = simplify(diff(O_omega_C4_C4));


%% Position Vectors

% rBO (B frame)
rBO_B = transpose([rBO_x_B(t) rBO_y_B(t) rBO_z_B(t)]);

% rBO (O frame)
rBO_O = VectRotation(OcB,rBO_B);

% rBO (CM frames)
rBO_C1 = VectRotation(C1cB,rBO_B);
rBO_C2 = VectRotation(C2cB,rBO_B);
rBO_C3 = VectRotation(C3cB,rBO_B);
rBO_C4 = VectRotation(C4cB,rBO_B);

% rPB (G frame)
rPB_G = [0; rCH; 0];

% rPB (O frame)
rPB_O = VectRotation(OcG,rPB_G);

% rPB (B frame)
rPB_B = VectRotation(BcO,rPB_O);

% leg hinge points from B (B frame)
rh1B_B = transpose([rh1B_x_B rh1B_y_B rh1B_z_B]);
rh2B_B = transpose([rh2B_x_B rh2B_y_B rh2B_z_B]);
rh3B_B = transpose([rh3B_x_B rh3B_y_B rh3B_z_B]);
rh4B_B = transpose([rh4B_x_B rh4B_y_B rh4B_z_B]);

% leg hinge points from B (CM frames)
rh1B_C1 = VectRotation(C1cB,rh1B_B);
rh2B_C2 = VectRotation(C2cB,rh2B_B);
rh3B_C3 = VectRotation(C3cB,rh3B_B);
rh4B_C4 = VectRotation(C4cB,rh4B_B);

% leg CM points from hinge points (CM frames)
rC1h1_C1 = transpose([rC1h1_x_C1 rC1h1_y_C1 rC1h1_z_C1]);
rC2h2_C2 = transpose([rC2h2_x_C2 rC2h2_y_C2 rC2h2_z_C2]);
rC3h3_C3 = transpose([rC3h3_x_C3 rC3h3_y_C3 rC3h3_z_C3]);
rC4h4_C4 = transpose([rC4h4_x_C4 rC4h4_y_C4 rC4h4_z_C4]);

% leg CM points from B (CM frames)
rC1B_C1 = simplify(rh1B_C1 + rC1h1_C1);
rC2B_C2 = simplify(rh2B_C2 + rC2h2_C2);
rC3B_C3 = simplify(rh3B_C3 + rC3h3_C3);
rC4B_C4 = simplify(rh4B_C4 + rC4h4_C4);

rC1B_B = VectRotation(BcC1,rC1B_C1);
rC2B_B = VectRotation(BcC2,rC2B_C2);
rC3B_B = VectRotation(BcC3,rC3B_C3);
rC4B_B = VectRotation(BcC4,rC4B_C4);

% leg CM points from O (CM frames)
rC1O_C1 = simplify(rBO_C1 + rC1B_C1);
rC2O_C2 = simplify(rBO_C2 + rC2B_C2);
rC3O_C3 = simplify(rBO_C3 + rC3B_C3);
rC4O_C4 = simplify(rBO_C4 + rC4B_C4);

%% Velocity Vectors

O_v_C1O_C1 = firstTransport(rC1O_C1,O_omega_C1_C1);
O_v_C2O_C2 = firstTransport(rC2O_C2,O_omega_C2_C2);
O_v_C3O_C3 = firstTransport(rC3O_C3,O_omega_C3_C3);
O_v_C4O_C4 = firstTransport(rC4O_C4,O_omega_C4_C4);

O_v_C1B_B = firstTransport(rC1B_B,O_omega_B_B);
O_v_C2B_B = firstTransport(rC2B_B,O_omega_B_B);
O_v_C3B_B = firstTransport(rC3B_B,O_omega_B_B);
O_v_C4B_B = firstTransport(rC4B_B,O_omega_B_B);

O_v_BO_B = firstTransport(rBO_B,O_omega_B_B);

%% Acceleration Vectors

O_a_BO_B = secondTransport(rBO_B, O_omega_B_B, O_alpha_B);

O_a_C1B_C1 = secondTransport(rC1B_C1, O_omega_C1_C1, O_alpha_C1);
O_a_C2B_C2 = secondTransport(rC2B_C2, O_omega_C2_C2, O_alpha_C2);
O_a_C3B_C3 = secondTransport(rC3B_C3, O_omega_C3_C3, O_alpha_C3);
O_a_C4B_C4 = secondTransport(rC4B_C4, O_omega_C4_C4, O_alpha_C4);

O_a_C1B_B = VectRotation(BcC1,O_a_C1B_C1);
O_a_C2B_B = VectRotation(BcC2,O_a_C2B_C2);
O_a_C3B_B = VectRotation(BcC3,O_a_C3B_C3);
O_a_C4B_B = VectRotation(BcC4,O_a_C4B_C4);

%% Moments of Inertia

% Chassis MOI (B frame)
I_Ch_B = [
    I_Ch_xx I_Ch_xy I_Ch_xz
    I_Ch_yx I_Ch_yy I_Ch_yz
    I_Ch_zx I_Ch_zy I_Ch_zz
    ];

% Leg MOIs (CM frames)
I_L1_C1 = [
    I_L_xx I_L_xy I_L_xz
    I_L_yx I_L_yy I_L_yz
    I_L_zx I_L_zy I_L_zz
    ];
I_L2_C2 = I_L1_C1;
I_L3_C3 = I_L1_C1;
I_L4_C4 = I_L1_C1;

%% Angular Momentum Components and their derivatives

% Components of O_O_h_Bsys, in their respective Frames ( _C1 suffix means
% expressed in C1 Frame)
gamma1_C1 = simplify((I_L1_C1 * O_omega_C1_C1) + transpose((mL * cross(transpose(rC1B_C1), transpose(O_v_C1O_C1)))));
gamma2_C2 = simplify((I_L2_C2 * O_omega_C2_C2) + transpose((mL * cross(transpose(rC2B_C2), transpose(O_v_C2O_C2)))));
gamma3_C3 = simplify((I_L3_C3 * O_omega_C3_C3) + transpose((mL * cross(transpose(rC3B_C3), transpose(O_v_C3O_C3)))));
gamma4_C4 = simplify((I_L4_C4 * O_omega_C4_C4) + transpose((mL * cross(transpose(rC4B_C4), transpose(O_v_C4O_C4)))));
gammaB_B = simplify(I_Ch_B * O_omega_B_B);

% take O-derivative of each term
d_gamma1_C1 = firstTransport(gamma1_C1, O_omega_C1_C1);
d_gamma2_C2 = firstTransport(gamma2_C2, O_omega_C2_C2);
d_gamma3_C3 = firstTransport(gamma3_C3, O_omega_C3_C3);
d_gamma4_C4 = firstTransport(gamma4_C4, O_omega_C4_C4);
d_gammaB_B = firstTransport(gammaB_B, O_omega_B_B);

% rotate them all into a common frame (B-frame)
d_gamma1_B = VectRotation(BcC1,d_gamma1_C1);
d_gamma2_B = VectRotation(BcC2,d_gamma2_C2);
d_gamma3_B = VectRotation(BcC3,d_gamma3_C3);
d_gamma4_B = VectRotation(BcC4,d_gamma4_C4);

% sum them together to get ddt_O_O_h_Bsys (in B frame)
ddt_O_O_h_Bsys_B = simplify(d_gamma1_B + d_gamma2_B + d_gamma3_B + d_gamma4_B + d_gammaB_B);

% additional term of the RHS of sum of torques expression
derp1 = mL * O_v_C1B_B;
derp2 = mL * O_v_C2B_B;
derp3 = mL * O_v_C3B_B;
derp4 = mL * O_v_C4B_B;

Mtotal = M + (4*mL);

O_v_CMB_B = simplify((derp1 + derp2 + derp3 + derp4) / Mtotal);

addlTerm_B = simplify(transpose(Mtotal*cross(transpose(O_v_BO_B), transpose(O_v_CMB_B))));

%%  Final sums of ext forces and sums of ext moments expressions (RHS)
%***************************************************************
%***************************************************************
%***************************************************************
sum_tau_Bsys_B = ddt_O_O_h_Bsys_B + addlTerm_B;

sum_F_Bsys_B = (Mtotal*O_a_BO_B) + mL*(O_a_C1B_B + O_a_C2B_B + O_a_C3B_B + ...
    O_a_C4B_B);


%% Derive force terms
% Gravitational Forces in the O frame
Fg_B_O = [0; -M*g; 0]; % g is negative, and jO is pointing down
Fg_C1_O = [0; -mL*g; 0];
Fg_C2_O = Fg_C1_O;
Fg_C3_O = Fg_C1_O;
Fg_C4_O = Fg_C1_O;

% Gravitational Forces in the B frame
Fg_B_B = VectRotation(BcO,Fg_B_O);
Fg_C1_B = VectRotation(BcO,Fg_C1_O);
Fg_C2_B = VectRotation(BcO,Fg_C2_O);
Fg_C3_B = VectRotation(BcO,Fg_C3_O);
Fg_C4_B = VectRotation(BcO,Fg_C4_O);

% Normal Force in the O frame, and B frame
Fn_P_G = [0; Fn; 0];
Fn_P_B = VectRotation(BcG,Fn_P_G);
Fn_P_O = VectRotation(OcG,Fn_P_G);

% Friction and Force Rolling Resistance (O frame)
Ffr_P_G = [Ffr; 0; 0];
Frr_P_G = [Frr; 0; 0];

% Friction and Force Rolling Resistance (B frame)
Ffr_P_B = VectRotation(BcG,Ffr_P_G);
Frr_P_B = VectRotation(BcG,Frr_P_G);

% Disturbance force (O-Frame and B-Frame)
Fd_B_G = [Fd; 0; 0];
Fd_B_B = VectRotation(BcG,Fd_B_G);

%% Derive Torque Terms
% gravitational torques from each leg (B frame)
tau_C1_B_B = simplify(transpose(cross(transpose(rC1B_B),transpose(Fg_C1_B))));
tau_C2_B_B = simplify(transpose(cross(transpose(rC2B_B),transpose(Fg_C2_B))));
tau_C3_B_B = simplify(transpose(cross(transpose(rC3B_B),transpose(Fg_C3_B))));
tau_C4_B_B = simplify(transpose(cross(transpose(rC4B_B),transpose(Fg_C4_B))));

% Friction Torque Term (B-frame)
tau_fr_B_B = simplify(transpose(cross(transpose(rPB_B),transpose(Ffr_P_B))));


%% Build Sum of Forces and Sum of Moments LHS
%***********************************************************
%***********************************************************
%***********************************************************
SumF_B = simplify(Fg_B_B + Fg_C1_B + Fg_C2_B + Fg_C3_B + Fg_C4_B + Fn_P_B + Ffr_P_B ...
    + Frr_P_B + Fd_B_B);

SumTau_B = simplify(tau_C1_B_B + tau_C2_B_B + tau_C3_B_B + tau_C4_B_B + tau_fr_B_B);


%% EXTRA SUPPORT STUFF
syms Crr Cf dthB_tigger
% velocity of B wrt to O, in the G frame
O_v_BO_G = VectRotation(GcB,O_v_BO_B);
% Finding a numerically integratable Crr (cuts out to zero as angular
% velocity gets really small)
Crr_bar = Crr*tanh(str2sym('diff(thetaB(t), t)')/dthB_tigger);
%RHS of equation for solving for rolling resistance  (G frame)
Frr_EQ_G = [(-Crr_bar * abs(Fn) * (O_v_BO_G(1) / abs(O_v_BO_G(1))))
    0
    0];
%RHS of equation for solving for rolling resistance  (B frame)
Frr_EQ_B = VectRotation(BcG,Frr_EQ_G);



%% Substitute time derivatives for variable names for solving
% new syms to replace the time-dependent vars and their derivatives
syms thB dthB ddthB...
    gam1 dgam1 ddgam1...
    gam2 dgam2 ddgam2...
    gam3 dgam3 ddgam3...
    gam4 dgam4 ddgam4...
    ...
    rBOx drBOx ddrBOx...
    rBOy drBOy ddrBOy...
    rBOz drBOz ddrBOz...
        
% vector holding the original syms to be replaced
regVect = [...
    str2sym('thetaB(t)') str2sym('diff(thetaB(t), t)') str2sym('diff(thetaB(t), t, t)') ...
    str2sym('gamma1(t)') str2sym('diff(gamma1(t), t)') str2sym('diff(gamma1(t), t, t)') ...
    str2sym('gamma2(t)') str2sym('diff(gamma2(t), t)') str2sym('diff(gamma2(t), t, t)') ...
    str2sym('gamma3(t)') str2sym('diff(gamma3(t), t)') str2sym('diff(gamma3(t), t, t)') ...
    str2sym('gamma4(t)') str2sym('diff(gamma4(t), t)') str2sym('diff(gamma4(t), t, t)') ...
    ...
    str2sym('rBO_x_B(t)') str2sym('diff(rBO_x_B(t), t)') str2sym('diff(rBO_x_B(t), t, t)') ...
    str2sym('rBO_y_B(t)') str2sym('diff(rBO_y_B(t), t)') str2sym('diff(rBO_y_B(t), t, t)') ...
    str2sym('rBO_z_B(t)') str2sym('diff(rBO_z_B(t), t)') str2sym('diff(rBO_z_B(t), t, t)') ...
    ...
    str2sym('M') str2sym('mL') str2sym('g') ...    % masses
    ...
    str2sym('rCH')...      % radius of chassis
    ...% vectors from B to each hinge point
    str2sym('rh1B_x_B') str2sym('rh1B_y_B') str2sym('rh1B_z_B') ...
    str2sym('rh2B_x_B') str2sym('rh2B_y_B') str2sym('rh2B_z_B') ...
    str2sym('rh3B_x_B') str2sym('rh3B_y_B') str2sym('rh3B_z_B') ...
    str2sym('rh4B_x_B') str2sym('rh4B_y_B') str2sym('rh4B_z_B')...
    ...% vectors from each hinge to each leg CM
    str2sym('rC1h1_x_C1') str2sym('rC1h1_y_C1') str2sym('rC1h1_z_C1') ...
    str2sym('rC2h2_x_C2') str2sym('rC2h2_y_C2') str2sym('rC2h2_z_C2') ...
    str2sym('rC3h3_x_C3') str2sym('rC3h3_y_C3') str2sym('rC3h3_z_C3') ...
    str2sym('rC4h4_x_C4') str2sym('rC4h4_y_C4') str2sym('rC4h4_z_C4') ... 
    ...% MOI Chassis
    str2sym('I_Ch_xx') str2sym('I_Ch_xy') str2sym('I_Ch_xz') ...
    str2sym('I_Ch_yx') str2sym('I_Ch_yy') str2sym('I_Ch_yz') ...
    str2sym('I_Ch_zx') str2sym('I_Ch_zy') str2sym('I_Ch_zz') ...
    ...MOI Each Leg
    str2sym('I_L_xx') str2sym('I_L_xy') str2sym('I_L_xz') ...
    str2sym('I_L_yx') str2sym('I_L_yy') str2sym('I_L_yz') ...
    str2sym('I_L_zx') str2sym('I_L_zy') str2sym('I_L_zz') ...
    ];

% vector of vars or vals to be substituted into their corresponding var in 
% regVect
subVect = [...
    str2sym('theB') str2sym('dthB') str2sym('ddthB') ...
    str2sym('gam1') str2sym('dgam1') str2sym('ddgam1') ...
    str2sym('gam2') str2sym('dgam2') str2sym('ddgam2') ...
    str2sym('gam3') str2sym('dgam3') str2sym('ddgam3') ...
    str2sym('gam4') str2sym('dgam4') str2sym('dgam4') ...
    ...
    str2sym('rBOx') str2sym('drBOx') str2sym('ddrBOx') ...
    str2sym('rBOy') str2sym('drBOy') str2sym('ddrBOy') ...
    0 0 0 ... Get rid of all z linear terms
    ...
    24 0.8720 -9.807 ...   % M, mL, g
    ...
    0.3937...      % radius of chassis
    ...  % vectors from B to each hinge point
    0.3364 0.1259 0 ... h1
    -0.3364 0.1259 0 ... h2
    -0.3364 -0.1259 0 ... h3
    0.3364 -0.1259 0 ...  h4   
    ... % vectors from each hinge to each leg CM
    -0.1967 0.0364 0 ... C1
    -0.1967 0.0364 0 ... C2
    -0.1967 0.0364 0 ... C3
    -0.1967 0.0364 0 ... C4
    ... % MOI Chassis
    0 0 0 ... X-row
    0 0 0 ... Y-row
    0 0 1.2550 ... Z-row
    ... MOI Each Leg
    0 0 0 ... X-row
    0 0 0 ... Y-row
    0 0 0.0082 ... z-row
    ];

% now use subs to replace variables for any expression
    %RHS
sum_F_Bsys_B_sub = simplify(subs(sum_F_Bsys_B,regVect,subVect));
sum_tau_Bsys_B_sub = simplify(subs(sum_tau_Bsys_B,regVect,subVect));
%rBO_O_sub = simplify(subs(rBO_O,regVect,subVect));
    %LHS
SumF_B_sub = simplify(subs(SumF_B,regVect,subVect));
SumTau_B_sub = simplify(subs(SumTau_B,regVect,subVect));
%rBO_O_RHS_sub = simplify(subs(rBO_O_RHS_sub,regVect,subVect));

Frr_EQ_B_sub = simplify(subs(Frr_EQ_B,regVect,subVect));


%% Solve equations for the state derivatives and simplify them
LHS_expr = [SumF_B_sub; SumTau_B_sub ; Frr_P_B];
RHS_expr = [sum_F_Bsys_B_sub; sum_tau_Bsys_B_sub ; Frr_EQ_B_sub];
varMatrix = [ddthB ddrBOx ddrBOy Frr Ffr Fn];

stateDeriv_obj = solve(LHS_expr == RHS_expr,varMatrix);

StateDeriv_orig = [
    simplify(stateDeriv_obj.ddthB);
    simplify(stateDeriv_obj.ddrBOx);
    simplify(stateDeriv_obj.ddrBOy);
    ];

ForceDeriv_orig = [
    simplify(stateDeriv_obj.Frr);
    simplify(stateDeriv_obj.Ffr);
    simplify(stateDeriv_obj.Fn);
    ];


%% Substitute again and put into state space formulation
% x(1) = thB        % x(5) = rBOy       % x(9) = gam2       % x(13) = gam4
% x(2) = dthB       % x(6) = drBOy      % x(10)= dgam2      % x(14) = dgam4
% x(3) = rBOx       % x(7) = gam1       % x(11)= gam3       
% x(4) = drBOx      % x(8) = dgam1      % x(12)= dgam3      
%*******************************************************

regVect = [ ...
    thB dthB ...
    rBOx drBOx ...
    rBOy drBOy ...
    gam1 dgam1 ddgam1 ...
    gam2 dgam2 ddgam2 ...
    gam3 dgam3 ddgam3 ...
    gam4 dgam4 ddgam4 ...
    ];

subVect = [ ...
    str2sym('x(1)') str2sym('x(2)') ...
    str2sym('x(3)') str2sym('x(4)') ...
    str2sym('x(5)') str2sym('x(6)') ...
    str2sym('x(7)') str2sym('x(8)') str2sym('u(1)') ...
    str2sym('x(9)') str2sym('x(10)') str2sym('u(2)') ...
    str2sym('x(11)') str2sym('x(12)') str2sym('u(3)') ...
    str2sym('x(13)') str2sym('x(14)') str2sym('u(4)') ...
    ];

EqVect_LHS = transpose([ ...
    str2sym('xdot(1)') str2sym('xdot(2)') str2sym('xdot(3)') ...
    str2sym('xdot(4)') str2sym('xdot(5)') str2sym('xdot(6)') ...
    str2sym('xdot(7)') str2sym('xdot(8)') str2sym('xdot(9)') ...
    str2sym('xdot(10)') str2sym('xdot(11)') str2sym('xdot(12)') ...
    str2sym('xdot(13)') str2sym('xdot(14)')
    ]);

ForceEQVect_LHS = transpose([ ...
    str2sym('Frr') str2sym('Ffr') str2sym('Fn') ...
    ]);

StateDeriv_SS = subs(StateDeriv_orig,regVect,subVect);

% build the system of state derivative equations as they would go in an
% m-file
EQVect_RHS = [
    str2sym('x(2)')
    StateDeriv_SS(1)
    str2sym('x(4)')
    StateDeriv_SS(2)
    str2sym('x(6)')
    StateDeriv_SS(3)
    str2sym('x(8)')
    str2sym('u(1)')
    str2sym('x(10)')
    str2sym('u(2)')
    str2sym('x(12)')
    str2sym('u(3)')
    str2sym('x(14)')
    str2sym('u(4)')
    ];

ForceEQVect_RHS = [
    ForceDeriv_orig(1)
    ForceDeriv_orig(2)
    ForceDeriv_orig(3)
    ];

EQString = append(string(EqVect_LHS), " = ", string(EQVect_RHS));
ForceEQ_string = append(string(ForceEQVect_LHS), " = ", string(ForceEQVect_RHS));


%% Finally, write the basic M-file
    FID = fopen(mFileName, 'w');
    
    fprintf(FID,"function [xdot] = TRREx_SimFile(t,x,u)\n\n");
    
    

    fprintf(FID,"\tFd  = getDisturbance(t,x,u);\n\n");
    
    for i = 1:size(ForceEQ_string,1)
        fprintf(FID,"\t%s;\n",ForceEQ_string(i,:));
    end
    fprintf(FID,"\n\n");
    
    for i = 1:size(EQString,1)
        fprintf(FID,"\t%s;\n",EQString(i,:));
    end
    fprintf(FID,"\n\n");
    
    fprintf(FID,"\txdot=xdot';\n\n");

    fprintf(FID,"end\n\n\n");
    fclose(FID);
    
%% Support Functions

function [rotatedVect] = VectRotation(rotationMat,vect2Rotate)
    rotatedVect = simplify(rotationMat * vect2Rotate);
end

function [derivVect] = firstTransport(q,omega)
    syms t;
    dq = simplify(diff(q,t));
    crossProd = simplify(transpose(cross(transpose(omega),transpose(q))));
    derivVect = simplify(dq + crossProd);
end

function [derivVect] = secondTransport(q,omega,alpha)
    syms t;
    ddq = simplify(diff(q,t,2));
    dq = simplify(diff(q,t));
    crossProd_1 = simplify(transpose(2*cross(transpose(omega),transpose(dq))));
    crossProd_2 = simplify(transpose(cross(transpose(alpha),transpose(q))));
    crossProd_3_inner = simplify(transpose(cross(transpose(omega),transpose(q))));
    crossProd_3 = simplify(transpose(cross(transpose(omega),transpose(crossProd_3_inner))));
    
    derivVect = simplify(ddq + crossProd_1 + crossProd_2 + crossProd_3);
end