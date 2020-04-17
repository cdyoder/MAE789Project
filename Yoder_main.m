function Yoder_main

    % yoder code to run TRREx simulations
    clear
    close all
    clc
    fclose('all');

    % derive eoms
    matfile = DeriveEOMS;


end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function [xd] = states(t, x0, params)
    % code to calculate the state derivatives of a trrex model
    
    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function [eoms_mat] = DeriveEOMS

    % define symbolics
    syms t thetaB(t) gamma1(t) gamma2(t) gamma3(t) gamma4(t)...
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
    WX(t) WY(t) WZ(t) ...
    WXd WYd WZd ...
    ...
    Fn Ffr Frr Fd...        % force terms which can be calculated in a different
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
    O_omega_B_B(1) = WX;
    O_omega_B_B(2) = WY; 
    O_omega_B_B(3) = WZ;
    eqn1 = simplify([O_omega_B_x_B; O_omega_B_y_B; O_omega_B_z_B]) - O_omega_B_B; % == 0
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

    % rBO (O frame)
    rBO_O = VectRotation(OcB, rBO_B);

    % rBO (CM frames)
    rBO_C1 = VectRotation(C1cB, rBO_B);
    rBO_C2 = VectRotation(C2cB, rBO_B);
    rBO_C3 = VectRotation(C3cB, rBO_B);
    rBO_C4 = VectRotation(C4cB, rBO_B);

    % rPB (O frame)
    rPB_O = [0; rCH; 0];

    % rPB (B frame)
    rPB_B = VectRotation(BcO, rPB_O);

    % leg hinge points from B (B frame)
    rh1B_B = sym('rh1B_B', [3, 1]);
    rh1B_B(1) = rh1B_x_B;
    rh1B_B(2) = rh1B_y_B;
    rh1B_B(3) = rh1B_z_B;
    
    rh2B_B = sym('rh2B_B', [3, 1]);
    rh2B_B(1) = rh2B_x_B;
    rh2B_B(2) = rh2B_y_B;
    rh2B_B(3) = rh2B_z_B;
    
    rh3B_B = sym('rh3B_B', [3, 1]);
    rh3B_B(1) = rh3B_x_B;
    rh3B_B(2) = rh3B_y_B;
    rh3B_B(3) = rh3B_z_B;
    
    rh4B_B = sym('rh4B_B', [3, 1]);
    rh4B_B(1) = rh4B_x_B;
    rh4B_B(2) = rh4B_y_B;
    rh4B_B(3) = rh4B_z_B;
    


    % leg hinge points from B (CM frames)
    rh1B_C1 = VectRotation(C1cB,rh1B_B);
    rh2B_C2 = VectRotation(C2cB,rh2B_B);
    rh3B_C3 = VectRotation(C3cB,rh3B_B);
    rh4B_C4 = VectRotation(C4cB,rh4B_B);

    % leg CM points from hinge points (CM frames)
    rC1h1_C1 = sym('rC1h1_C1', [3, 1]);
    rC1h1_C1(1) = rC1h1_x_C1;
    rC1h1_C1(2) = rC1h1_y_C1;
    rC1h1_C1(3) = rC1h1_z_C1;
    
    rC2h2_C2 = sym('rC2h2_C2', [3, 1]);
    rC2h2_C2(1) = rC2h2_x_C2;
    rC2h2_C2(2) = rC2h2_y_C2;
    rC2h2_C2(3) = rC2h2_z_C2;
   
    rC3h3_C3 = sym('rC3h3_C3', [3, 1]);
    rC3h3_C3(1) = rC3h3_x_C3;
    rC3h3_C3(2) = rC3h3_y_C3;
    rC3h3_C3(3) = rC3h3_z_C3;
    
    rC4h4_C4 = sym('rC4h4_C4', [3, 1]);
    rC4h4_C4(1) = rC4h4_x_C4;
    rC4h4_C4(2) = rC4h4_y_C4;
    rC4h4_C4(3) = rC4h4_z_C4;
    
    
    
    
    
    
    

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
    
    disp('Position done');  
    % keyboard
    
    
    
    % TO DO
    % 1. Finish EOMs 
    % 2. Do simple test cases to prove functionality in 2D
    % 3. Use either time or angle dependant controller to demonstrate
    % openloop
    
    
    
    
    
%     % [OvBO] = Transport1(OwB, rBO, t)
%     O_v_C1O_C1 = Transport1(O_omega_C1_C1, rC1O_C1,t);
%     O_v_C2O_C2 = Transport1(O_omega_C2_C2, rC2O_C2, t);
%     O_v_C3O_C3 = Transport1(O_omega_C3_C3, rC3O_C3, t);
%     O_v_C4O_C4 = Transport1(O_omega_C4_C4, rC4O_C4, t);

    O_v_C1B_B = Transport1(O_omega_B_B, rC1B_B, t);
    O_v_C2B_B = Transport1(O_omega_B_B, rC2B_B, t);
    O_v_C3B_B = Transport1(O_omega_B_B, rC3B_B, t);
    O_v_C4B_B = Transport1(O_omega_B_B, rC4B_B, t);
    O_v_BO_B = Transport1(O_omega_B_B, rBO_B, t);

    
    % [OaBO] = Transport2(OwB, OaB, rBO, t)
    O_a_BO_B = Transport2(O_omega_B_B, O_alpha_B, rBO_B, t);

    O_a_C1B_C1 = Transport2(O_omega_C1_C1, O_alpha_C1, rC1B_C1, t);
    O_a_C2B_C2 = Transport2(O_omega_C2_C2, O_alpha_C2, rC2B_C2, t);
    O_a_C3B_C3 = Transport2(O_omega_C3_C3, O_alpha_C3, rC3B_C3, t);
    O_a_C4B_C4 = Transport2(O_omega_C4_C4, O_alpha_C4, rC4B_C4, t);

    O_a_C1B_B = VectRotation(BcC1, O_a_C1B_C1);
    O_a_C2B_B = VectRotation(BcC2, O_a_C2B_C2);
    O_a_C3B_B = VectRotation(BcC3, O_a_C3B_C3);
    O_a_C4B_B = VectRotation(BcC4, O_a_C4B_C4);
    
    disp('Acceleration done');  
    keyboard
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    syms t th(t) r1(t) r2(t) r3(t) r4(t) L(t)
    syms m Ixx g R M M_ch mu

    bco = [1 0 0;0 cos(th) sin(th);0 -sin(th) cos(th)];
    ocb = [1 0 0;0 cos(th) -sin(th);0 sin(th) cos(th)];
    ocdotb = diff(ocb,t);
    % Position and Angular Velocities 
    % B Frame
    rm1b = [0;r1;0];
    rm2b = [0;-r2;0];
    rm3b = [0;0;r3];
    rm4b = [0;0;-r4];

    % O Frame
    rbo_o = [0;L;0]; % Linear distance moved by TW from O frame. 
    rpb_o = [0;0;-R]; % Radius of TW (0.75m)

    % B Frame
    rbo = bco*rbo_o;
    rpb = bco*rpb_o;
    rcmb = M*(rm1b+rm2b+rm3b+rm4b)/(4*M + M_ch);
    rcmo = rbo + rcmb;

    wx = [0 0 1]*bco*ocdotb*[0;1;0];
    wy = [1 0 0]*bco*ocdotb*[0;0;1];
    wz = [0 1 0]*bco*ocdotb*[1;0;0];
    owb = simplify([wx;wy;wz]);

    % Linear velocities and accelerations


    % ovbo_o = diff(rbo_o,t); % io jo ko
    ovbo_o = cross(owb,rpb_o);


    ovbo = bco*ovbo_o;


    ovcmb = diff(rcmb, t) + cross(owb,rcmb);

    ovm1b = diff(rm1b,t)+cross(owb,rm1b);
    ovm2b = diff(rm2b,t)+cross(owb,rm2b);
    ovm3b = diff(rm3b,t)+cross(owb,rm3b);
    ovm4b = diff(rm4b,t)+cross(owb,rm4b);

    ovm1o = ovm1b + ovbo;
    ovm2o = ovm2b + ovbo;
    ovm3o = ovm3b + ovbo;
    ovm4o = ovm4b + ovbo;

    ovcmo = ovbo + ovcmb;  % ib jb kb
    ovcmo_o = ocb*ovcmo;

    % Acceleration
    oabo_o = diff(ovbo_o,t); %io jo ko
    oacmo = diff(ovcmo,t);
    oabo = diff(ovbo,t);

    %ang Momentum

    ohbody = Ixx*owb;
    ohm1b = M*cross(rm1b,(ovm1o));
    ohm2b = M*cross(rm2b,(ovm2o));
    ohm3b = M*cross(rm3b,(ovm3o));
    ohm4b = M*cross(rm4b,(ovm4o));

    ohbtotal = simplify(ohbody+ohm1b+ohm2b+ohm3b+ohm4b);

    %Torque

    tau_sys_o = simplify(diff(ohbtotal,t)+cross(owb,ohbtotal)+(4*M+M_ch)*cross(ovbo_o,ovcmo_o),50);

    % O_F_mass = [0;0;M*g];
    % taum1 = cross(rm1b,bco*O_F_mass);
    % taum2 = cross(rm2b,bco*O_F_mass);
    % taum3 = cross(rm3b,bco*O_F_mass);
    % taum4 = cross(rm4b,bco*O_F_mass);
    % 
    % tau_ext=taum1+taum2+taum3+taum4;

    O_F_mass = [0;0;-M*g];

    taum1 = cross(rm1b,bco*O_F_mass);
    taum2 = cross(rm2b,bco*O_F_mass);
    taum3 = cross(rm3b,bco*O_F_mass);
    taum4 = cross(rm4b,bco*O_F_mass);

    O_F_Kin_fric = [0;mu*(4*M+M_ch)*g;0]; % Fr = mu * N
    taufric = cross(rpb_o,O_F_Kin_fric);


    % O_F_Roll_resist = Nb/r
    % O_F_tot = O_F_mass + O_F_Kin_fric;

    syms c
    tau_c = [-c*diff(th(t), t); 0; 0];

    tau_ext = taum1 + taum2 + taum3 + taum4 + taufric + tau_c;
    tau_sys = bco*tau_sys_o;
    eq1=tau_sys==tau_ext;
    syms tta thdot thddot R1 R2 R3 R4 R1dot R2dot R3dot R4dot e1
    temp=subs(eq1,[th(t),diff(th(t), t),diff(th(t), t,t),r1(t),r2(t),r3(t),r4(t),diff(r1(t), t),diff(r2(t), t),diff(r3(t), t),diff(r4(t), t)],[tta,thdot,thddot,R1,R2,R3,R4,R1dot,R2dot,R3dot,R4dot])
    e1=solve(temp,thddot)

end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function [prod] = CrossMe(v1, v2)
    % function to perform the cross product on symbolics
    prod = sym('prod', [3, 1]);
    prod(1) = v1(2)*v2(3) - v1(3)*v2(2);
    prod(2) = v1(3)*v2(1) - v1(1)*v2(3);
    prod(3) = v1(1)*v2(2) - v1(2)*v2(1);
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function [OvBO] = Transport1(OwB, rBO, t)
    % first derivative transport theorem
    OvBO = simplify(diff(rBO, t) + CrossMe(OwB, rBO));
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function [OaBO] = Transport2(OwB, OaB, rBO, t)
    % second derivative transport theorem
    OaBO = simplify(diff(rBO, t, t) + 2*CrossMe(OwB, diff(rBO, t)) + ...
        CrossMe(OaB, rBO) + CrossMe(OwB, CrossMe(OwB, rBO)));
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function [rotatedVect] = VectRotation(rotationMat,vect2Rotate)
    rotatedVect = simplify(rotationMat * vect2Rotate);
end