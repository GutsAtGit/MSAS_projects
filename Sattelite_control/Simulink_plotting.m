% Plot Similink-Simscape results
clc; close all;
% Parameters
omega_0 = 1.160758e-3;
Torb_0 = 2*pi/omega_0;

% Extract simulink data
% time 
t = ss_x_v.Time;
% accelerometer
x_a = ss_x_a.Data;
v_a = ss_v_a.Data;
% solenoid
x_v = ss_x_v.Data;
v_v = ss_v_v.Data;
Vhat_out = ss_Vhat_out.Data;
% drag and thrust
D = ss_D.Data;
T = ss_T.Data;

% Drag and Thrust diff
DminusT = D - T;


% Plot results
% plot x_a
figure;
subplot(2, 3, 1);
plot(t(2:end) / Torb_0, x_a(2:end), 'LineWidth', 1.5);
xlabel('Time (Orbits)'); ylabel('Seismic Position (mm)');
legend('x_a');
grid on; title('Accelerometer Dynamics');

% plot v_a
subplot(2, 3, 2);
plot(t(2:end) / Torb_0, v_a(2:end), 'LineWidth', 1.5);
xlabel('Time (Orbits)'); ylabel('Seismic Velocity (mm/s)');
lim_va = 1e-7; ylim([-lim_va, lim_va]); xlim([0, inf]);
legend('v_a');
grid on; title('Accelerometer Dynamics');

% plot x_v
subplot(2, 3, 3);
plot(t(2:end) / Torb_0, x_v(2:end), 'LineWidth', 1.5);
xlabel('Time (Orbits)'); ylabel('Spool Position (mm)');
legend('x_v');
grid on; title('Solenoid Dynamics');

% plot v_v
subplot(2, 3, 4);
plot(t(2:end) / Torb_0, v_v(2:end), 'LineWidth', 1.5);
lim_vv = 3e-6; ylim([-lim_vv, lim_vv]); xlim([0, inf]);
xlabel('Time (Orbits)'); ylabel('Spool velocity (mm/s)');
legend('v_v');
grid on; title('Solenoid Dynamics');

% plot Vhat_out (Op-Amp output)
subplot(2, 3, 5);
plot(t(2:end) / Torb_0, Vhat_out(2:end), 'LineWidth', 1.5);
xlabel('Time (Orbits)'); ylabel('Solenoid Voltage (mV)');
lim_Vol = 0.1; ylim([-lim_Vol, lim_Vol]); xlim([0, inf]);
legend('Op-Amp output');
grid on; title('Solenoid Dynamics');

% Drag vs Thrust
figure;
subplot(1, 2, 1)
hold on;
plot(t(2:end) / Torb_0, D(2:end), 'LineWidth', 1.5);
plot(t(2:end) / Torb_0, T(2:end), 'LineWidth', 1.5);
xlabel('Time (Orbits)'); ylabel('Drag vs Thrust (mN)');
ylim([0.4, 2]);
legend('Drag', 'Thrust');
hold off; grid on; title('Drag Compensation');

subplot(1, 2, 2)
plot(t(2:end) / Torb_0, DminusT(2:end), 'LineWidth', 1.5);
xlabel('Time (Orbits)'); ylabel('Drag vs Thrust (mN)');
lim_T = 1.5e-4; ylim([-lim_T, lim_T]); xlim([0, inf]);
legend('Drag - Thrust');
grid on; title('Drag and Thrust difference');