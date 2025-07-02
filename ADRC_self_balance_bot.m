%% ADRC control system implemented on self balance bot
% System Parameters
M = 10; % mass of body (kg)
k = 500; % equivalent spring constant of the damping spring (N/m)
c = 50; % equivalent damping coefficient (Ns/m)
y_eq = 0.3; % equilibrium vertical position (m)
g = 9.81; % gravity (m/s^2)

% Physical Constraints
y_min_physical = 0.2; % minimum height (m)
y_max_physical = 0.5; % maximum height (m)

% ADRC Parameters
b0 = 1 / M; % Control gain
omega_c = 3; % Controller bandwidth (rad/s)
omega_o = 10 * omega_c; % Observer bandwidth (rad/s)

% ESO gains
beta1 = 3 * omega_o;
beta2 = 3 * omega_o^2;
beta3 = omega_o^3;

% Controller gain for the error dynamics
K1 = omega_c^2;
K2 = 2 * omega_c;

% fal function parameters
alpha = 0.5;    % controls the nonlinearity
delta = 0.01;   % defines the linear region around zero error

% Simulation Parameters
T_sim = 5;
dt = 0.001;
t = 0:dt:T_sim;
N = length(t);

% Set Point
y_d_setpoint = input('Enter desired height (m): ');
yd_dot_setpoint = input('Enter desired vertical velocity (m/s): ');
yd_ddot_setpoint = input('Enter desired vertical accelaration (m/s^2): ');
y_d_vec = y_d_setpoint * ones(1, N);
yd_dot_vec = yd_dot_setpoint * ones(1, N);
yd_ddot_vec = yd_ddot_setpoint * ones(1, N);

% Initial Conditions
height_0 = input('Enter intial Height (m): ');
vel_0 = input('Enter initial vertical velocity (m/s): ');
x_sys = zeros(2, N); %[height , vertical velocity]
x_sys(:, 1) = [height_0; vel_0]; 

% Initialize the ESO states
x_eso = zeros(3, N); %[height , vertical velovity , total disturbance]
x_eso(:, 1) = [x_sys(1, 1); x_sys(2, 1); 0]; % ESO states initialized

% Initialize array to store the control input
u_control = zeros(1, N);

% Main Simulation Loop
%fal function definition
fal_func = @(e_val, a, d) (abs(e_val) > d) .* (sign(e_val) .* abs(e_val).^a) + (abs(e_val) <= d) .* (e_val / (d^(1-a)));

for i = 1:N-1
    % Getting current system output and ESO estimates
    y_current = x_sys(1, i); % current actual height

    z1 = x_eso(1, i); % current estimated height from ESO
    z2 = x_eso(2, i); % current estimated velocity from ESO
    z3 = x_eso(3, i); % current estimated total disturbance from ESO

    % Getting set point values for the current time step
    y_d = y_d_vec(i);
    yd_dot = yd_dot_vec(i);
    yd_ddot = yd_ddot_vec(i);

    % Calculate tracking errors
    e0 = y_d - z1;
    e1 = yd_dot - z2;

    % Apply fal function to errors to generate nonlinear error feedback
    f0 = fal_func(e0, alpha, delta);
    f1 = fal_func(e1, alpha, delta);

    % Calculate Control Input
    u = (K1 * f0 + K2 * f1 - z3) / b0;

    % Apply saturation to control input
    u_max = 500;
    u_min = -500;
    u = max(min(u, u_max), u_min); % Clamp the control input within limits

    u_control(i) = u; % Store the calculated control input

    % System Dynamics Updatew
    y_dot_current = x_sys(2, i); % current actual velocity

    % Calculate the net force on the system
    F_spring = -k * (y_current - y_eq); % Spring force
    F_gravity = -M * g;                 % Gravity force
    F_damping = -c * y_dot_current;     % Damping force
    F_net_other = F_spring + F_gravity + F_damping;

    % Calculate the actual acceleration of the system
    y_ddot_actual = (u + F_net_other) / M;

    % Integrate system states using Euler method
    x_sys(1, i+1) = x_sys(1, i) + x_sys(2, i) * dt;
    x_sys(2, i+1) = x_sys(2, i) + y_ddot_actual * dt;

    % Apply physical limits to height
    x_sys(1, i+1) = max(y_min_physical, x_sys(1, i+1));
    x_sys(1, i+1) = min(y_max_physical, x_sys(1, i+1));

    % ESO Dynamics Update
    % Observer error 
    eso_error = y_current - z1;

    % ESO state derivatives
    z1_dot = z2 + beta1 * eso_error;        % Derivative of estimated height
    z2_dot = z3 + beta2 * eso_error + b0 * u; % Derivative of estimated velocity
    z3_dot = beta3 * eso_error;             % Derivative of estimated total disturbance

    % Integrate ESO states using Euler method
    x_eso(1, i+1) = z1 + z1_dot * dt; % Update estimated height
    x_eso(2, i+1) = z2 + z2_dot * dt; % Update estimated velocity
    x_eso(3, i+1) = z3 + z3_dot * dt; % Update estimated total disturbance
end
fprintf('Simulation finished.\n');

% Plot Results
figure('Name', 'Simulation Results');

% Position Tracking
subplot(3,1,1);
plot(t, y_d_vec, 'k--', 'DisplayName', 'Desired Height'); hold on;
plot(t, x_sys(1, :), 'b', 'DisplayName', 'Actual Height');
plot(t, x_eso(1, :), 'r--', 'DisplayName', 'Estimated Height');
grid on;
legend('Location', 'northeastoutside');
xlabel('Time (s)');
ylabel('Position (m)');
title('Height Tracking');
box on;

% Estimated Velocity
subplot(3,1,2);
plot(t, x_sys(2, :), 'b', 'DisplayName', 'Actual Velocity');
hold on;
plot(t, x_eso(2, :), 'r--', 'DisplayName', 'Estimated Velocity');
grid on;
legend('Location', 'northeastoutside');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Estimated Velocity');
box on;

% Estimated Total Disturbance and Control Input
subplot(3,1,3);
yyaxis left;
plot(t, x_eso(3, :), 'g', 'DisplayName', 'Estimated Disturbance');
ylabel('Disturbance (N/kg)');

yyaxis right;
plot(t, u_control, 'm', 'DisplayName', 'Control Input');
ylabel('Force (N)');

grid on;
legend('Location', 'northeastoutside');
xlabel('Time (s)');
title('Estimated Disturbance and Control Input');
box on;

%final errors
final_error = abs(y_d_vec(end) - x_sys(1, end));
fprintf('Final steady-state position error: %.4f m\n', final_error);
