% Modeling and Simulation of Aerospace systems
% Assignment #2
% Exercise #2
% Author: Tarun Singh

%% EX2_part1
clc; close all; clearvars;

% Struct for the parameters

params = struct;

% Accelerometer parameters
params.M_sc = 300;      % Spacecraft mass [kg]
params.m_a = 0.32;      % Seismic mass [kg]
params.b_a = 1.5e3;     % Accelerometer damper [Ns/m]; range = (1.5e3 - 2e4)
params.k_a = 5e-5;      % Accelerometer spring [N/m]; range = (5e-5 - 3e-3)
params.K_acc = 1;       % Accelerometer proportional coefficient [Vs/m]

% Amplifier parameters
params.R_in = 0.1;      % Inverting resistance [ohm]; range = (0.1 - 10)
params.R_f = 8e4;       % Feedback resistance [ohm]; range = (1e4 - 8e4)

% Solenoid Valve parameters
params.m_v = 0.1;       % Spool mass [kg]
params.k_v = 1e3;       % Valve spring [N/m]
params.b_v = 1e3;       % Valve damper [Ns/m]
params.alpha = 2.1e-2;  % Solenoid constant [1/H]
params.beta = -60;      % Solenoid gain [1/Hm]
params.A_0 = 4.7e-12;   % Minimum aperture area [m^2]
params.l = 1e-5;        % Width of the flap [m]
params.x_vmax = 1e-5;   % Maximum spool extension [m]

% Thrusters parameters
params.k = 1.66;            % Heat ratio []
params.p_T = 2e5;           % Tank pressure [Pa]
params.T_T = 240;           % Tank temperature [K]
params.R_Xe = 63.32754;     % Xenon gas constant [J/kg-K]
params.q = 1.6e-19;         % Xenon ion charge [C]
params.delta_V = 2000;      % Xenon ion voltage [V]
params.m_i = 2.188e-25;     % Xenon ion mass [kg]

% Drag parameters
params.omega_s = 1.658226e-6; % Secular pulsation [rad/s]
params.omega_0 = 1.160758e-3; % Orbital pulsation [rad/s] 

% initial Conditions
x0 = [0, 0, 0, 0, 0];
Torb_0 = 2*pi/params.omega_0;
tspan = linspace(0,3*Torb_0, 500);
%tspan = [0 3*Torb_0];


% Solve ode system
% solve the ODE system
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);
[t, x] = ode15s(@(t, x) odefun(t, x, params), tspan, x0, options);
% Extract results

x_a = x(2:end, 1); v_a = x(2:end, 2);                   % Accelerometer Dynamics
x_v = x(2:end, 3); v_v = x(2:end, 4); I = x(2:end, 5);  % Solenoid Dynamics
V_out = params.K_acc*v_a;                               % V_out (Op-Amp input)
Vhat_out = -(params.R_f/params.R_in)*V_out;             % Vhat_out (Op-Amp output)

% Drag
D = 2.2 - cos(params.omega_s .* tspan(2:end)) + 1.2 * sin(params.omega_0 .* tspan(2:end)) .* cos(params.omega_0 .* tspan(2:end));

% Thrust calculation from state variables
% l = params.x_vmax;
A_v = params.A_0 + params.l .* (params.x_vmax - x_v);
rho_Xe = params.p_T / (params.R_Xe * params.T_T);
k = params.k;
mdot_Xe = A_v .* sqrt(k .* rho_Xe .* params.p_T .* (2 / (k + 1))^((k + 1) / (k - 1)));
v_exit = sqrt(2 * params.q * params.delta_V / params.m_i);
T = mdot_Xe .* v_exit;

% Drag and Thrust difference
DminusT = D' - T*1e3;

% Plot results
% plot x_a
figure;
subplot(2, 3, 1);
plot(t(2:end) / Torb_0, x_a, 'LineWidth', 1.5);
xlabel('Time (Orbits)'); ylabel('Seismic Position (m)');
legend('x_a');
grid on; title('Accelerometer Dynamics');

% plot v_a
subplot(2, 3, 2);
plot(t(2:end) / Torb_0, v_a, 'LineWidth', 1.5);
xlabel('Time (Orbits)'); ylabel('Seismic Velocity (m/s)');
legend('v_a');
grid on; title('Accelerometer Dynamics');

% plot x_v
subplot(2, 3, 3);
plot(t(2:end) / Torb_0, x_v, 'LineWidth', 1.5);
xlabel('Time (Orbits)'); ylabel('Spool Position (m)');
legend('x_v');
grid on; title('Solenoid Dynamics');

% plot v_v
subplot(2, 3, 4);
plot(t(2:end) / Torb_0, v_v, 'LineWidth', 1.5);
xlabel('Time (Orbits)'); ylabel('Spool velocity (m/s)');
legend('v_v');
grid on; title('Solenoid Dynamics');

% plot I
subplot(2, 3, 5);
plot(t(2:end) / Torb_0, I, 'LineWidth', 1.5);
xlabel('Time (Orbits)'); ylabel('Solenoid Current (Amp)');
legend('I');
grid on; title('Solenoid Dynamics');

% plot Vhat_out (Op-Amp output)
subplot(2, 3, 6);
plot(t(2:end) / Torb_0, Vhat_out, 'LineWidth', 1.5);
xlabel('Time (Orbits)'); ylabel('Solenoid Voltage (V)');
legend('Op-Amp output');
grid on; title('Solenoid Dynamics');

% Drag vs Thrust
figure;
subplot(1, 2, 1)
hold on;
plot(t(2:end) / Torb_0, D, 'LineWidth', 1.5);
plot(t(2:end) / Torb_0, T*1e3, 'LineWidth', 1.5);
xlabel('Time (Orbits)'); ylabel('Drag vs Thrust (mN)');
legend('Drag', 'Thrust');
hold off; grid on; title('Drag Compensation');

subplot(1, 2, 2)
plot(t(2:end) / Torb_0, DminusT, 'LineWidth', 1.5);
xlabel('Time (Orbits)'); ylabel('Drag vs Thrust (mN)');
legend('Drag - Thrust');
grid on; title('Drag and Thrust difference');

% System of non-linear ODEs
function dxdt = odefun(t, x, params)

    % Extract parameters
    % Accelerometer
    M_sc = params.M_sc;
    m_a = params.m_a;
    b_a = params.b_a;
    k_a = params.k_a;
    K_acc = params.K_acc;
    
    % Amplifier
    R_in = params.R_in;
    R_f = params.R_f;
    
    % Solenoid Valve
    m_v = params.m_v;
    k_v = params.k_v;
    b_v = params.b_v;
    alpha = params.alpha;
    beta = params.beta;
    A_0 = params.A_0;
    x_vmax = params.x_vmax;
    
    % Thruster
    k = params.k;
    p_T = params.p_T;
    T_T = params.T_T;
    R_Xe = params.R_Xe;
    q = params.q;
    delta_V = params.delta_V;
    m_i = params.m_i;
    
    % Drag
    omega_s = params.omega_s;
    omega_0 = params.omega_0;

    % State Variables
    % Accelerometer
    x_a = x(1);     % Seismic mass position
    v_a = x(2);     % Seismic mass velocity

    % Solenoid valve
    x_v = x(3);     % valve spool position
    v_v = x(4);     % Valve spool velocity
    I = x(5);       % Solenoid current
    
    % Drag force
    D = 2.2 - cos(omega_s * t) + 1.2 * sin(omega_0 * t) * cos(omega_0 * t);
    
    % Voltage modulation
    V_out = K_acc * v_a;
    V_hat_out = -R_f / R_in * V_out;
    
    % Solenoid Inductance, Force and aperture area
    l = x_vmax;
    L = 1 / (alpha + beta * x_v);
    dL_dxv = -beta / (alpha + beta * x_v)^2;
    A_v = A_0 + l * (x_vmax - x_v);

    % Xenon mass flow rate
    rho_Xe = p_T / (R_Xe * T_T);
    mdot_Xe = A_v * sqrt(k * rho_Xe * p_T * (2 / (k + 1))^((k + 1) / (k - 1)));

    % Thrust
    v_exit = sqrt(2 * q * delta_V / m_i);
    T = mdot_Xe * v_exit;

    % Accelerometer dynamics
    dx_adt = v_a; 
    dv_adt = (T - D * 1e-3) / M_sc - (1 / m_a) * (b_a * v_a + k_a * x_a);

    % Solenoid valve dynamics
    dx_vdt = v_v;
    dv_vdt = (1 / m_v) * (-k_v * x_v - b_v * v_v + 0.5 * I^2 * dL_dxv);

    % Current
    dIdt = V_hat_out / L;

    % Collect derivatives
    dxdt = [dx_adt; dv_adt; dx_vdt; dv_vdt; dIdt];
end
