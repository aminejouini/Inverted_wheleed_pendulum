%% pendulum model parameters
m_p = 0.267;        % full pendulum mass 
m_pen = 0.216;      % motor+wheel  
l_cm  = 0.111;      % center of mass distance (Solidworks)
J_w = 205971.26e-9      ;    % Wheel momentum of inertia (Solidworks)
g = 9.81;           % gravity acceleration
J_p =  1949047.64e-9 ;   % full pendulum momentum of inertia (Solidworks) m_p*l_cm^2 + 1153811.05e-9

%% motor parameters
c_m = 0.024;        % motor torque constant
c_e = 6.3/325;      % motor back EMF constant
R_a = 10.21;        % motor armature resistance

% controller sampling time
Ts = 50e-3;          % controller sampling time

% power supply
MAX_U = 12;        % V

%% linear system equations
A = [0          1       0;...
    g/l_cm      0       c_e*c_m/(J_p*R_a);...
    -g/l_cm     0       -c_m*c_e/(J_w*R_a)];

B = [0;...
    -c_m/(J_p*R_a);...
    c_m/(J_w*R_a)];

C = [1 0 0];

D = [0];

% discretization
G = ss(A,B,C,D);
Gd = c2d(G,Ts);

%% LQR controller
Qd = diag([1000,1,1]);
Rd = 1;

[Kr, P] = dlqr(Gd.a,Gd.b,Qd,Rd);
Kr =-Kr;  

% Compute closed-loop poles
Ac = Gd.a - Gd.b*Kr;
Bc = Gd.b;
Cc = Gd.c;
Dc = Gd.d;
sys_cl = ss(Ac,Bc,Cc,Dc);
cl_poles = pole(sys_cl);

% Plot closed-loop poles
figure;
pzmap(sys_cl);
title('Closed-loop pole locations');




% copy this to the code in the file pendulum.c
disp(strcat('/** LQR Matlab generated controller'))
disp(strcat('* Q = diag([',num2str(Qd(1,1)),',',num2str(Qd(2,2)),',',num2str(Qd(3,3)),']) R = [',num2str(Rd),']'))
disp(strcat('* Closed-loop poles: [',num2str(cl_poles.'),']'))
disp(strcat('*/'))
disp(strcat('float Kr[3] = {',num2str(Kr(1)),',',num2str(Kr(2)),',',num2str(Kr(3)),'};'))