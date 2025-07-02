%% LQR controller implementation for self balancing bot 

% Define system parameters
M = 1;    % Mass of the chassis (kg)
m = 0.5741;  % Mass of the Wheels and legs (kg)
b = 0.1;  % Coefficient of friction (Ns/m)
I_x = 0.0358; % Moment of Inertia (X-Axis) (kg*m^2)
I_y = 0.0197; % Moment of Inertia (Y-Axis) (kg*m^2)
I_z = 0.0268; % Moment of Inertia (Z-Axis) (kg*m^2)
I   = 0.0489; % Moment of Inertia (Total) (kg*m^2)
g = 9.8;  % Acceleration due to gravity (m/s^2)
l = 0.158;  % Length to wheel axis from the Center of Mass
p = I*(M+m)+M*m*l^2; % Denominator for the A and B matrices

% Define State-Space Matrices A and B
A = [0      1              0           0;
     0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
     0      0              0           1;
     0 -(m*l*b)/p       m*g*l*(M+m)/p  0];

B = [     0;
     (I+m*l^2)/p;
          0;
        m*l/p];

% Check Controllability Rank
a = rank(ctrb(A,B));
disp(['Rank of Controllability Matrix: ', num2str(a)]);
if a == 4
    fprintf('LQR can be applied, since Rank is equal to the Dimensions of Matrix A\n');
else
    fprintf('LQR cannot be applied, as Rank is not equal to the Dimensions of Matrix A\n');
end

% Define Output Matrices C and D
C = [1 0 0 0;
     0 0 1 0]; 
D=[0;0];

pos_0 = input('Enter initial position (m): ');
vel_0 = input('Enter initial velovity (m/s): ');
theta_0 = input('Enter the initial angle (rad): ');
angular_vel_0 = input('Enter initial angular velocity (rad/s): ');

% Define Initial State Vector
x0 = [pos_0; vel_0; theta_0; angular_vel_0];

% Define LQR Weighting Matrices Q and R

% Q penalizes state errors
Q = [1     0     0     0 ;       % Penalty on position error
     0     25    0     0 ;       % Penalty on velocity error
     0     0     25    0 ;       % Penalty on angle error
     0     0     0     25];     % Penalty on angular velocity error

% R penalizes control effort
R = 1;

% Calculate LQR Gain Matrix K
K = lqr(A,B,Q,R);
fprintf("Value of K:\n");
disp(K);

ltiSys = ss((A-B*K), B, C, D);

% Simulate the System Response
t = 0:0.01:50;
[y, t, x_states] = lsim(ltiSys, zeros(size(t)), t, x0);

% Display Poles of the System
fprintf("Poles of System (Must be negative real parts for stability):\n");
disp(pole(ltiSys));

% Plot all states versus time
figure;
subplot(4,1,1); 
plot(t, x_states(:,1)); 
ylabel('Position (m)');
title('System States over Time');
grid on;

subplot(4,1,2); 
plot(t, x_states(:,2)); 
ylabel('Velocity (m/s)');
grid on;

subplot(4,1,3);
plot(t, x_states(:,3)); 
ylabel('Angle (Radians)');
grid on;

subplot(4,1,4);
plot(t, x_states(:,4));
ylabel('Angular Velocity (rad/s)');
xlabel('Time (seconds)');
grid on;
