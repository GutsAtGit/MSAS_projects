% Modeling and Simulation of Aerospace systems
% Assignment #2
% Exercise #1
% Author: Tarun Singh

%% EX1_part1 (Causal Modeling)
clc; close all; clearvars;

% Struct for the parameters

params = struct;

% Physical constants
params.sigma = 5.67e-8;  % Stefan-Boltzmann constant [W/m^2 K^4]
params.T_ds = 3;         % Deep space temperature [K]
params.P_Sun = 1350;     % Solar constant [W/m^2]

% Characteristics of themal control system (DC motor and Radiators)
params.R = 0.1;        % Resistance [Ohm]
params.L = 0.001;      % Inductance [H]
params.k_m = 0.3;      % Motor Constant [Nm/A]
params.m_r = 0.2;      % Radiator mass [kg]
params.L_r = 0.5;      % Radiator length [m]

% Thermal parameters of the system
A_r1 = 2 * (1.5 * 0.5) + 2 * (1.49 * 0.5) + (0.5*0.5);  % Radiative area of Body1
A_r2 = (0.95 * 0.5) + 2 * (0.01 * 0.95) + (0.01 * 0.5); % = A_r3; Radiative area of Body 2 & 3
A_r4 = (0.5 * 0.5);                                     % Radiative area of Body 4 & 5
A_s1 = (0.5 * 0.5);                                     % Body 1 area exposed to sun
A_s2 = (0.95 * 0.5);                                    % Body 2 & 3 area exposed to sun
params.A_r = [A_r1, A_r2, A_r4];                % Areas corresponding to radiation(A1, A2 = A3, A4 = A5) [m^2]
params.A_a = [A_s1, A_s2];                      % Areas corresponding to absorption (A1, A2 = A3) [m^2]
params.C = [1.5e5, 1187.5, 1187.5, 30, 30];     % Heat Capacities (C1, C2, C3, C4, C5) [J/K]
params.G = [10, 10, 10, 10];                    % Thermal Conductance (G12, G13, G14, G15) [W/K]
params.alpha = [0.6, 0.78, 0.78];               % Absorptivity (alpha1,2,3) 
params.epsilon = [0.45,0.75, 0.75];             % Emissivity (epsilon1,2,3)
params.epsilon_min = 0.01;                      % Radiator Emissivity (min)
params.epsilon_max = 0.98;                      % Radiator Emissivity (max)

% Radiator control parameters
params.theta_min = -0.4 * pi;   % Minimum radiator angle (closed)
params.theta_max = 0;           % Maximum radiator angle (open)
params.T1_ref = 294.15;         % Body 1 reference temperature [K]
params.T1_min = 290;            % Body 1 Minimum allowable temperature [K]
params.T1_max = 300;            % Body 1 Maximum allowable temperature [K]

params.k_p = 4e-5;              % Gain constant


% Initial conditions
x0 = [298.15, 298.15, 298.15, 298.15, 298.15, -0.4 * pi, 0, 0];     % Initial temperatures (K), radiator angle (radians), angular velocity (rad/s), current (i)
tspan = [0 3600*50];                                                % Simulation time (50 hours)

% solve the ODE system
options = odeset('RelTol', 1e-11, 'AbsTol', 1e-12);
[t, x] = ode15s(@(t, x) odefun(t, x, params), tspan, x0, options);

% Extract results
T1 = x(:, 1); T2 = x(:, 2); T3 = x(:, 3); T4 = x(:, 4); T5 = x(:, 5);
theta = x(:, 6);
theta_dot = x(:, 7);
current = x(:, 8);

% Plot results
% plot T1
figure;
subplot(2, 2, 1);
plot(t / 3600, T1, 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Temperature (K)');
xline(10, '--c', 'LineWidth', 2);                                                   % t = 10 hrs
yline(params.T1_min, '--r', 'LineWidth', 2);                                        % limits
yline(params.T1_max, '--r', 'LineWidth', 2);
yline([293.85585 294.44415], '--m', {' ', '0.1% of T_1^{ref}'}, 'Linewidth', 2);    % 0.1% range of T1_ref
yline(294.15, '--k','Linewidth', 2);                                                % T1_ref
legend('Body 1 temperature');
grid on; title('Thermal Response');

% Plot Angle evolution
subplot(2, 2, 2);
plot(t / 3600, theta, 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Radiator Angle (rad)');
xline(10, '--c', 'LineWidth', 2);                               % t = 10 hrs
yline([0 -0.4*pi], '--r', 'LineWidth', 2);                      % limits
legend('\theta');
grid on; title('Control Effort');

% Plot Angular velocity
subplot(2, 2, 3);
plot(t / 3600, theta_dot, 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Radiator Angular velocity (rad/s)');
xline(10, '--c', 'LineWidth', 2);                               % t = 10 hrs
legend('$\dot{\theta}$', 'Interpreter','latex');
grid on; title('Control Speed');

% Plot Current
subplot(2, 2, 4);
plot(t / 3600, current, 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Current (Amp)');
xline(10, '--c', 'LineWidth', 2);                               % t = 10 hrs
lim = 0.8e-8; ylim([-lim, lim]); xlim([0, inf]);
legend('Current');
grid on; title('Current in E-M system');

% Plot all temperatures
figure;
hold on;
plot(t / 3600, T1, 'LineWidth', 2);
plot(t / 3600, T2, 'LineWidth', 2);
plot(t / 3600, T3, '--','LineWidth', 2);
plot(t / 3600, T4, 'LineWidth', 2);
plot(t / 3600, T5, '--', 'LineWidth', 2);
xlabel('Time (hours)'); ylabel('Temperature (K)');
xline(10, '--c', 'LineWidth', 1);                                                       % t = 10 hrs
yline([params.T1_min params.T1_max], '--r', {'T_{min}', 'T_{max}'}, 'LineWidth', 1);    % limits
yline([293.85585 294.44415], '--m', {' ', '0.1% of T_1^{ref}'}, 'Linewidth', 1);        % 0.1% range of T1_ref
yline(294.15, '--k','Linewidth', 1);                                                    % T1_ref
legend('T_1 (Main body)', 'T_2 (Solar Panel 1)', 'T_3 (Solar Panel 2)', 'T_4 (Radiator 1)', 'T_5 (Radiator 2)', 'Location', 'best');
hold off; grid on; title('Thermal Response');

% System of non-linear ODEs
function dxdt = odefun(t, x, params)

    % Extract parameters
    % Physical constants
    sigma = params.sigma;
    T_ds = params.T_ds;
    P_sun = params.P_Sun;
    
    % Control parameters
    R = params.R;
    L = params.L;
    k_m = params.k_m;
    m_r = params.m_r;
    L_r = params.L_r;
    
    % Thermal paramters
    A_r = params.A_r;
    A_a = params.A_a;
    C = params.C;
    G = params.G;
    alpha = params.alpha;
    epsilon = params.epsilon;
    epsilon_min = params.epsilon_min;
    epsilon_max = params.epsilon_max;
    theta_max = params.theta_max;
    theta_min = params.theta_min;
    
    % Radiator control
    T1_ref = params.T1_ref;
    k_p = params.k_p;


    % Radiator moment of inertia
    J_r = (1/3) * m_r * L_r^2;

    % State variables
    % Temperature variables
    T1 = x(1); T2 = x(2); T3 = x(3); T4 = x(4); T5 = x(5);

    % theta variable
    theta = x(6); theta_dot = x(7);

    % current variable
    i = x(8);

    % Input voltage and theta limitation relation
    V_in = k_p*(T1 - T1_ref);

    % Control strategy for Voltage input
    if (theta > theta_max && V_in > 0)
        V_in = 0;                                           % stop power supply
        theta = theta_max;                                  % maximum angle acheieved
    elseif (theta < theta_min && V_in < 0)
        V_in = 0;                                           % stop Power supply
        theta = theta_min;                                  % minimum angle acheieved
    end

    % Radiator emissivity
    epsilon_theta = epsilon_min + ((epsilon_max - epsilon_min) / (0.4*pi)) *(theta + 0.4*pi);

    % Heat from sun power
    Q_sun1 = P_sun * A_a(1) * alpha(1);
    Q_sun2 = P_sun * A_a(2) * alpha(2);      % = Q_sun3

    % Radiation to deep space
    Q_rad1 = epsilon(1) * sigma * A_r(1) * (T1^4 - T_ds^4);
    Q_rad2 = epsilon(2) * sigma * A_r(2) * (T2^4 - T_ds^4);
    Q_rad3 = epsilon(3) * sigma * A_r(2) * (T3^4 - T_ds^4);
    Q_rad4 = epsilon_theta * sigma * A_r(3) * (T4^4 - T_ds^4);
    Q_rad5 = epsilon_theta * sigma * A_r(3) * (T5^4 - T_ds^4);

    % Define ODEs
    % Tempearature ODEs
    dT1dt = (G(1) * (T2 - T1) + G(2) * (T3 - T1) + G(3) * (T4 - T1) + G(4) * (T5 - T1) ...
        + Q_sun1 - Q_rad1) / C(1);
    dT2dt = (G(1) * (T1 - T2) + Q_sun2 - Q_rad2) / C(2);
    dT3dt = (G(2) * (T1 - T3) + Q_sun2 - Q_rad3) / C(3);
    dT4dt = (G(3) * (T1 - T4) - Q_rad4) / C(4);
    dT5dt = (G(4) * (T1 - T5) - Q_rad5) / C(5);

    % Angle ODEs (radiator dynamics)
    dtheta_dt = theta_dot;
    dtheta_dot_dt = (k_m * i)/J_r;

    % Current ODE
    di_dt = (V_in - i*R - k_m*theta_dot)/L;
    
    % Collect derivatives
    dxdt = [dT1dt; dT2dt; dT3dt; dT4dt; dT5dt; dtheta_dt; dtheta_dot_dt; di_dt];
end