function [mFileName] = DeriveEOMS(basename)

    % define symbolics
    syms t thetaB(t) gamma1(t) gamma2(t) gamma3(t) gamma4(t)...
    ...
    M mL ...    % masses
    ...
    g ...       % gravitational constant
    ...
    rCH...      % radius of chassis
    ...
    rBO_x_B(t) ... % Vector from O to B
    rh1B_x_B rh1B_y_B ...
    rh2B_x_B rh2B_y_B ...
    rh3B_x_B rh3B_y_B ...
    rh4B_x_B rh4B_y_B...     % vectors from B to each hinge point
    rC1h1_x_C1 rC1h1_y_C1 ...
    rC2h2_x_C2 rC2h2_y_C2 ...
    rC3h3_x_C3 rC3h3_y_C3 ...
    rC4h4_x_C4 rC4h4_y_C4 ... % vectors from each hinge to each leg CM
    ...
    I_Ch_zz ...
    ...
    I_L_zz ...
    ...
    WZ(t) WZd WZdd mCH ...
    ...
    Fn Ffr Frr Fd Crr...        % force terms which can be calculated in a different
                            % function
                        
    % rotation matrices
    OcB = sym('OCB', [3, 3]);
    OcB(1, 1) = cos(thetaB);
    OcB(1, 2) = -sin(thetaB);
    OcB(1, 3) = 0;
    OcB(2, 1) = sin(thetaB);
    OcB(2, 2) = cos(thetaB);
    OcB(2, 3) = 0;
    OcB(3, 1) = 0;
    OcB(3, 2) = 0;
    OcB(3, 3) = 1;
    
    BcC1 = sym('BcC1', [3, 3]);
    BcC1(1, 1) = cos(gamma1);
    BcC1(1, 2) = -sin(gamma1);
    BcC1(1, 3) = 0;
    BcC1(2, 1) = sin(gamma1);
    BcC1(2, 2) = cos(gamma1);
    BcC1(2, 3) = 0;
    BcC1(3, 1) = 0;
    BcC1(3, 2) = 0;
    BcC1(3, 3) = 1;
    
    BcC2 = sym('BcC2', [3, 3]);
    BcC2(1, 1) = cos(gamma2);
    BcC2(1, 2) = -sin(gamma2);
    BcC2(1, 3) = 0;
    BcC2(2, 1) = sin(gamma2);
    BcC2(2, 2) = cos(gamma2);
    BcC2(2, 3) = 0;
    BcC2(3, 1) = 0;
    BcC2(3, 2) = 0;
    BcC2(3, 3) = 1;
    
    BcC3 = sym('BcC3', [3, 3]);
    BcC3(1, 1) = cos(gamma3);
    BcC3(1, 2) = -sin(gamma3);
    BcC3(1, 3) = 0;
    BcC3(2, 1) = sin(gamma3);
    BcC3(2, 2) = cos(gamma3);
    BcC3(2, 3) = 0;
    BcC3(3, 1) = 0;
    BcC3(3, 2) = 0;
    BcC3(3, 3) = 1;

    BcC4 = sym('BcC4', [3, 3]);
    BcC4(1, 1) = cos(gamma4);
    BcC4(1, 2) = -sin(gamma4);
    BcC4(1, 3) = 0;
    BcC4(2, 1) = sin(gamma4);
    BcC4(2, 2) = cos(gamma4);
    BcC4(2, 3) = 0;
    BcC4(3, 1) = 0;
    BcC4(3, 2) = 0;
    BcC4(3, 3) = 1;

    
    
    OcC1 = simplify( OcB * BcC1 );
    OcC2 = simplify( OcB * BcC2 );
    OcC3 = simplify( OcB * BcC3 );
    OcC4 = simplify( OcB * BcC4 );

    BcO = transpose(OcB);
    C1cB = transpose(BcC1);
    C2cB = transpose(BcC2);
    C3cB = transpose(BcC3);
    C4cB = transpose(BcC4);
    C1cO = transpose(OcC1);
    C2cO = transpose(OcC2);
    C3cO = transpose(OcC3);
    C4cO = transpose(OcC4);
    disp('Rotation matrices done');    
    
    
    
    
    
    
    
    OcB_dot = diff(OcB,t);
    BcC1_dot = diff(BcC1,t);
    BcC2_dot = diff(BcC2,t);
    BcC3_dot = diff(BcC3,t);
    BcC4_dot = diff(BcC4,t);
    OcC1_dot = diff(OcC1,t);
    OcC2_dot = diff(OcC2,t);
    OcC3_dot = diff(OcC3,t);
    OcC4_dot = diff(OcC4,t);

    % BcO is matrix from B to O
    % used with OwB velocities
    % from definition, OwB = []*BcO*OcB_d*[]
    mat0 = BcO*OcB_dot;
    O_omega_B_x_B = [0 0 1]*mat0*[0;1;0];
    O_omega_B_y_B = [1 0 0]*mat0*[0;0;1];
    O_omega_B_z_B = [0 1 0]*mat0*[1;0;0];
    O_omega_B_B = sym('O_omega_B_B', [3, 1]);
    O_omega_B_B(1) = 0;
    O_omega_B_B(2) = 0; 
    % O_omega_B_B(3) = WZ;
    O_omega_B_B(3) = O_omega_B_z_B;
    % eqn1 = simplify([O_omega_B_x_B; O_omega_B_y_B; O_omega_B_z_B]) - O_omega_B_B; % == 0
    O_alpha_B_B = simplify(diff(O_omega_B_B, t) + CrossMe(O_omega_B_B, O_omega_B_B));
    

    % BwC1 = []*C1cB*BcC1_dot*[]
    mat1 = C1cB*BcC1_dot;
    B_omega_C1_C1 = sym('B_omega_C1_C1', [3, 1]);
    B_omega_C1_C1(1) = [0 0 1]*mat1*[0;1;0];
    B_omega_C1_C1(2) = [1 0 0]*mat1*[0;0;1];
    B_omega_C1_C1(3) = [0 1 0]*mat1*[1;0;0];
    B_omega_C1_B = simplify(BcC1*B_omega_C1_C1);
    B_alpha_C1_B = simplify(diff(B_omega_C1_B, t) + CrossMe(B_omega_C1_B, B_omega_C1_B));

    % BwC1 = []*C1cB*BcC1_dot*[]
    mat2 = C2cB*BcC2_dot;
    B_omega_C2_C2 = sym('B_omega_C2_C2', [3, 1]);
    B_omega_C2_C2(1) = [0 0 1]*mat2*[0;1;0];
    B_omega_C2_C2(2) = [1 0 0]*mat2*[0;0;1];
    B_omega_C2_C2(3) = [0 1 0]*mat2*[1;0;0];
    B_omega_C2_B = simplify(BcC2*B_omega_C2_C2);
    B_alpha_C2_B = simplify(diff(B_omega_C2_B, t) + CrossMe(B_omega_C2_B, B_omega_C2_B));

    % BwC1 = []*C1cB*BcC1_dot*[]
    mat3 = C3cB*BcC3_dot;
    B_omega_C3_C3 = sym('B_omega_C3_C3', [3, 1]);
    B_omega_C3_C3(1) = [0 0 1]*mat3*[0;1;0];
    B_omega_C3_C3(2) = [1 0 0]*mat3*[0;0;1];
    B_omega_C3_C3(3) = [0 1 0]*mat3*[1;0;0];
    B_omega_C3_B = simplify(BcC3*B_omega_C3_C3);
    B_alpha_C3_B = simplify(diff(B_omega_C3_B, t) + CrossMe(B_omega_C3_B, B_omega_C3_B));

    % BwC1 = []*C1cB*BcC1_dot*[]
    mat4 = C4cB*BcC4_dot;
    B_omega_C4_C4 = sym('B_omega_C4_C4', [3, 1]);
    B_omega_C4_C4(1) = [0 0 1]*mat4*[0;1;0];
    B_omega_C4_C4(1) = [1 0 0]*mat4*[0;0;1];
    B_omega_C4_C4(1) = [0 1 0]*mat4*[1;0;0];
    B_omega_C4_B = simplify(BcC4*B_omega_C4_C4);
    B_alpha_C4_B = simplify(diff(B_omega_C4_B, t) + CrossMe(B_omega_C4_B, B_omega_C4_B));
    
    % additions
    O_omega_C1_B = O_omega_B_B + B_omega_C1_B;
    O_omega_C2_B = O_omega_B_B + B_omega_C2_B;
    O_omega_C3_B = O_omega_B_B + B_omega_C3_B;
    O_omega_C4_B = O_omega_B_B + B_omega_C4_B;   
    
    disp('Omega and alpha done');    

    
    
    
    
    
    % rBO (B frame
    rBO_B = sym('rBO_B', [3, 1]);
%     rBO_B(1) = rBO_x_B; % == R*theta
%     rBO_B(2) = rBO_y_B;
%     rBO_B(3) = rBO_z_B;
    rBO_B(1) = rCH*thetaB; % == R*theta
    rBO_B(2) = 0;
    rBO_B(3) = 0; 

    % % rBO (O frame)
    % rBO_O = VectRotation(OcB, rBO_B);

    % % rBO (CM frames)
    % rBO_C1 = VectRotation(C1cB, rBO_B);
    % rBO_C2 = VectRotation(C2cB, rBO_B);
    % rBO_C3 = VectRotation(C3cB, rBO_B);
    % rBO_C4 = VectRotation(C4cB, rBO_B);

    % rPB (O frame)
    rPB_O = [0; rCH; 0];

    % rPB (B frame)
    rPB_B = VectRotation(BcO, rPB_O);

    % leg hinge points from B (B frame)
    rh1B_B = sym('rh1B_B', [3, 1]);
    rh1B_B(1) = rh1B_x_B;
    rh1B_B(2) = rh1B_y_B;
    rh1B_B(3) = 0;
    
    rh2B_B = sym('rh2B_B', [3, 1]);
    rh2B_B(1) = rh2B_x_B;
    rh2B_B(2) = rh2B_y_B;
    rh2B_B(3) = 0;
    
    rh3B_B = sym('rh3B_B', [3, 1]);
    rh3B_B(1) = rh3B_x_B;
    rh3B_B(2) = rh3B_y_B;
    rh3B_B(3) = 0;
    
    rh4B_B = sym('rh4B_B', [3, 1]);
    rh4B_B(1) = rh4B_x_B;
    rh4B_B(2) = rh4B_y_B;
    rh4B_B(3) = 0;
    




    % leg CM points from hinge points (CM frames)
    rC1h1_C1 = sym('rC1h1_C1', [3, 1]);
    rC1h1_C1(1) = rC1h1_x_C1;
    rC1h1_C1(2) = rC1h1_y_C1;
    rC1h1_C1(3) = 0;
    
    rC2h2_C2 = sym('rC2h2_C2', [3, 1]);
    rC2h2_C2(1) = rC2h2_x_C2;
    rC2h2_C2(2) = rC2h2_y_C2;
    rC2h2_C2(3) = 0;
   
    rC3h3_C3 = sym('rC3h3_C3', [3, 1]);
    rC3h3_C3(1) = rC3h3_x_C3;
    rC3h3_C3(2) = rC3h3_y_C3;
    rC3h3_C3(3) = 0;
    
    rC4h4_C4 = sym('rC4h4_C4', [3, 1]);
    rC4h4_C4(1) = rC4h4_x_C4;
    rC4h4_C4(2) = rC4h4_y_C4;
    rC4h4_C4(3) = 0;
    
    
    
    
    
    
    

    % % leg CM points from B (CM frames)
    rC1B_B = rh1B_B + BcC1*rC1h1_C1;
    rC2B_B = rh2B_B + BcC2*rC2h2_C2;
    rC3B_B = rh3B_B + BcC3*rC3h3_C3;
    rC4B_B = rh4B_B + BcC4*rC4h4_C4;
    rC1B_C1 = C1cB*rC1B_B;
    rC2B_C2 = C2cB*rC2B_B;
    rC3B_C3 = C3cB*rC3B_B;
    rC4B_C4 = C4cB*rC4B_B;
    rC1O_B = rBO_B + rC1B_B;
    rC2O_B = rBO_B + rC2B_B;
    rC3O_B = rBO_B + rC3B_B;
    rC4O_B = rBO_B + rC4B_B;
    mT = 4*mL + mCH;
    rCMB_B = mL*(rC1B_B + rC2B_B + rC3B_B + rC4B_B)/(mT);

    disp('Position done');  
    % keyboard
    
    
    
    
    
    
    
    
%     % [OvBO] = Transport1(OwB, rBO, t)
    O_v_C1B_B = Transport1(O_omega_B_B, rC1B_B, t);
    O_v_C2B_B = Transport1(O_omega_B_B, rC2B_B, t);
    O_v_C3B_B = Transport1(O_omega_B_B, rC3B_B, t);
    O_v_C4B_B = Transport1(O_omega_B_B, rC4B_B, t);
    O_v_BO_B = Transport1(O_omega_B_B, rBO_B, t);
    O_v_CMB_B = Transport1(O_omega_B_B, rCMB_B, t);
    O_v_CMO_B = O_v_CMB_B + O_v_BO_B;
    
    % [OaBO] = Transport2(OwB, OaB, rBO, t)
    O_a_BO_B = Transport2(O_omega_B_B, O_alpha_B_B, rBO_B, t);
    O_a_CMB_B = Transport2(O_omega_B_B, O_alpha_B_B, rCMB_B, t);
    O_a_CMO_B = O_a_CMB_B + O_a_BO_B;
    
    disp('Acceleration done');  
    % keyboard
    
    
    
    
    
    % inertia tensors
    
    % Chassis MOI (B frame)
    I_Ch_B = sym('I_Ch_B', [3, 3]);
    I_Ch_B(1, 1) = 0;
    I_Ch_B(1, 2) = 0;
    I_Ch_B(1, 3) = 0;
    I_Ch_B(2, 1) = 0;
    I_Ch_B(2, 2) = 0;
    I_Ch_B(2, 3) = 0;
    I_Ch_B(3, 1) = 0;
    I_Ch_B(3, 2) = 0;
    I_Ch_B(3, 3) = I_Ch_zz;

    % Leg MOIs (CM frames)
    I_L1_C1 = sym('I_L1_C1', [3, 3]);
    I_L1_C1(1, 1) = 0;
    I_L1_C1(1, 2) = 0;
    I_L1_C1(1, 3) = 0;
    I_L1_C1(2, 1) = 0;
    I_L1_C1(2, 2) = 0;
    I_L1_C1(2, 3) = 0;
    I_L1_C1(3, 1) = 0;
    I_L1_C1(3, 2) = 0;
    I_L1_C1(3, 3) = I_L_zz;
    
    I_L2_C2 = I_L1_C1;
    I_L3_C3 = I_L1_C1;
    I_L4_C4 = I_L1_C1;
    
    % rotate
    I_L1_B = BcC1*I_L1_C1*C1cB; % B frame
    I_L2_B = BcC2*I_L2_C2*C2cB;
    I_L3_B = BcC3*I_L3_C3*C3cB;
    I_L4_B = BcC4*I_L4_C4*C4cB;
    
    
    
    
    
    
    % angular momentum
    OhOB_L1_B = CrossMe(rC1B_B, mL*O_v_C1B_B) + I_L1_B*O_omega_B_B;
    OhOB_L2_B = CrossMe(rC2B_B, mL*O_v_C2B_B) + I_L2_B*O_omega_B_B;
    OhOB_L3_B = CrossMe(rC3B_B, mL*O_v_C3B_B) + I_L3_B*O_omega_B_B;
    OhOB_L4_B = CrossMe(rC4B_B, mL*O_v_C4B_B) + I_L4_B*O_omega_B_B;
    OhOB_B_B = I_Ch_B*O_omega_B_B;
    OhOB_sys_B = OhOB_L1_B + OhOB_L2_B + OhOB_L3_B + OhOB_L4_B + OhOB_B_B;
    Oddt_OhOB_sys_B = Transport1(O_omega_B_B, OhOB_sys_B, t);
    VxV = CrossMe(O_v_BO_B, mT*O_v_CMO_B);
    disp('Momentum done');
    
    
    
    
    
    
    

    % sum forces
    Fg_CM_O = [0; mT*g; 0]; % weight
    Fn_P_O = [0; Fn; 0];    % normal
    Ffr_P_O = [Ffr; 0; 0];  % friction
    % Frr_P_O = [Frr; 0; 0];  % rolling resistance
    Frr_P_O = [-Crr*Fn; 0; 0];
    Fd_B_O = [Fd; 0; 0];    % disturbance
    eqn2 = BcO*(Fg_CM_O + Fn_P_O + Ffr_P_O + Frr_P_O + Fd_B_O) - mT*O_a_CMO_B;  % == 0, SumF = ma
    
    
    
    % sum torques
    Tg_CM_B = CrossMe(rCMB_B, BcO*Fg_CM_O); % weight
    Trr_P_B = BcO*CrossMe(rPB_O, Frr_P_O);
    eqn3 = Tg_CM_B + Trr_P_B - Oddt_OhOB_sys_B - VxV; % == 0, SumT = OddtOhOB + VxV
    disp('Forces and torques done');
    
    
    % solve and sub
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
        str2sym('mCH')];

    % vector of vars or vals to be substituted into their corresponding var in 
    % regVect
    subVect = [...
        str2sym('theB') str2sym('dtheB') str2sym('ddtheB') ...
        str2sym('gam1') str2sym('dgam1') str2sym('ddgam1') ...
        str2sym('gam2') str2sym('dgam2') str2sym('ddgam2') ...
        str2sym('gam3') str2sym('dgam3') str2sym('ddgam3') ...
        str2sym('gam4') str2sym('dgam4') str2sym('ddgam4') ...
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
        24];
    
    % % eqn1 = [0; 0; ###] so ignore first two elements
    % BIGEQNS(1) = simplify(subs(eqn1(3), regVect, subVect));
    
    % eqn2 = [##; ##; 0];
    % keyboard
    BIGEQNS(1) = simplify(subs(eqn2(1), regVect, subVect));
    BIGEQNS(2) = simplify(subs(eqn2(2), regVect, subVect));
    
    % eqn3 = [0; 0; ##]
    BIGEQNS(3) = simplify(subs(eqn3(3), regVect, subVect));
    
    % variables: thetadd, Fn, Ffr
    % equations: SumF (2x), SumT (1x)
    % keyboard
    
    % solve
    sobj = solve([BIGEQNS(1) == 0, BIGEQNS(2) == 0, BIGEQNS(3) == 0], ...
        [str2sym('ddtheB'), str2sym('Fn'), str2sym('Ffr')]);
    StateDeriv_orig = [simplify(sobj.ddtheB), simplify(sobj.Fn), simplify(sobj.Ffr)];
    
    % 
    
    
    
    
    
    
    % Substitute again and put into state space formulation
    % x(1) = thB        % x(5) = gam2       % x(9) = gam4
    % x(2) = dthB       % x(6) = dgam2      % x(10) = dgam4
    % x(3) = gam1       % x(7) = gam3       
    % x(4) = dgam1      % x(8) = dgam3      
    %*******************************************************

    regVect = [ ...
        str2sym('theB') str2sym('dtheB') ...
        str2sym('gam1') str2sym('dgam1') str2sym('ddgam1') ...
        str2sym('gam2') str2sym('dgam2') str2sym('ddgam2') ...
        str2sym('gam3') str2sym('dgam3') str2sym('ddgam3') ...
        str2sym('gam4') str2sym('dgam4') str2sym('ddgam4') ...
        ];

    subVect = [ ...
        str2sym('x(1)') str2sym('x(2)') ...
        str2sym('x(3)') str2sym('x(4)') str2sym('u(1)') ...
        str2sym('x(5)') str2sym('x(6)') str2sym('u(2)') ...
        str2sym('x(7)') str2sym('x(8)') str2sym('u(3)') ...
        str2sym('x(9)') str2sym('x(10)') str2sym('u(4)') ...
        ];

    EqVect_LHS = transpose([ ...
        str2sym('xdot(1)') str2sym('xdot(2)') str2sym('xdot(3)') ...
        str2sym('xdot(4)') str2sym('xdot(5)') str2sym('xdot(6)') ...
        str2sym('xdot(7)') str2sym('xdot(8)') str2sym('xdot(9)') ...
        str2sym('xdot(10)') ...
        ]);

    StateDeriv_SS = subs(StateDeriv_orig, regVect, subVect);

    % build the system of state derivative equations as they would go in an
    % m-file
    EQVect_RHS = [
    str2sym('x(2)')     % thetad = 
    StateDeriv_SS(1)    % thetadd = 
    str2sym('x(4)')     % gam1 d
    str2sym('u(1)')     % gam1 dd
    str2sym('x(6)')
    str2sym('u(2)')
    str2sym('x(8)')
    str2sym('u(3)')
    str2sym('x(10)')
    str2sym('u(4)')
    ];
    % keyboard
    EQString = append(string(EqVect_LHS), " = ", string(EQVect_RHS));


    % Finally, write the basic M-file
    % keyboard
    mFileName = [basename, '.m'];
    FID = fopen(mFileName, 'w');
    
    fprintf(FID,"function [xdot] = %s(t,x,Crr_nom,th_trig,g1,g2,g3,g4,flag)\n\n", basename);
    % fprintf(FID,"\tFd = getForces(t,x,u,Crr,thtrig);\n\n");
    
    fprintf(FID,"%% Variable schedule\n");
    fprintf(FID,"%% x(1) = thB        %% x(5) = gam2       %%  x(9) = gam4\n");
    fprintf(FID,"%% x(2) = dthB       %% x(6) = dgam2      %% x(10) = dgam4\n");
    fprintf(FID,"%% x(3) = gam1       %% x(7) = gam3\n");
    fprintf(FID,"%% x(4) = dgam1      %% x(8) = dgam3\n\n");  
    
    fprintf(FID, "%% Preallocate\n");
    fprintf(FID, "xdot = NaN(length(x), 1);\n\n");
    
    fprintf(FID, "%% Mass inputs\n");
    fprintf(FID, "u = MassSchedule(t, g1, g2, g3, g4);\n\n");
    
    fprintf(FID, "%% Coeffs:\n");
    fprintf(FID,"Crr=%s;\t%% Coefficient\n\n", "Crr_nom*tanh(x(2)/abs(th_trig))");
    
    fprintf(FID, "%% Forces:\n");
    fprintf(FID, "Fd=%s;\t%% Disturbance force\n", "GetDisturbance(t, x)");
    fprintf(FID, "Fn=%s;\t%% Normal force\n", string(StateDeriv_SS(2)));
    fprintf(FID, "Ffr=%s;\t%% Friction force\n", string(StateDeriv_SS(3)));
    fprintf(FID, "Fd=0;\t%% Disturbance force\n\n");  
    
    for i = 1:size(EQString,1)
        fprintf(FID,"%s;\n",EQString(i,:));
    end
    fprintf(FID,"\n\n");
    
    % code to get other parameters out
    fprintf(FID,"if flag == 1\n");
    fprintf(FID,"\txdn = xdot;\n");
    fprintf(FID,"\txdot = struct;\n");
    fprintf(FID,"\txdot.xd = xdn;\n");
    fprintf(FID,"\txdot.Frr = -Crr*Fn;\n");
    fprintf(FID,"\txdot.Fn = Fn;\n");
    fprintf(FID,"\txdot.Fd = Fd;\n");
    fprintf(FID,"\txdot.Ffr = Ffr;\n");
    fprintf(FID,"end\n\n");

    fprintf(FID,"end\n\n\n");
    fclose(FID);
    
    
end
