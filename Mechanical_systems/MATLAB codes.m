% Modeling and Simulation of Aerospace Systems (2024-2025)
% Assingment #1
% Author: Tarun Singh; 10894156

%% EX 1
clearvars; close all; clc; clear;


a1 = 10; a2 = 13; a3 = 8; a4 = 10;          % Rod lengths    
beta = linspace(0, 2*pi/3, 69);             % beta                                         
C = (a1^2 + a2^2 - a3^2 + a4^2)/(2*a2*a4);  % Constant 'C'
h_val = 1e-6;                               % step size for Finite Difference

% function and its analytical derivative
falpha = @(alpha, beta) (a1/a2)*cos(beta) - (a1/a4)*cos(alpha) - cos(beta - alpha) + C;
dfalpha = @(alpha, beta) (a1/a4)*sin(alpha) - sin(beta - alpha);

% Initializing vectors
fzero_alpha = zeros(2, length(beta));
alpha = zeros(2, length(beta));
iterations = zeros(2, length(beta));
success = false(2, length(beta));

% Initializing Vectors for Finite Difference
alpha_FD = zeros(2, length(beta));
iterations_FD = zeros(2, length(beta));
success_FD = false(2, length(beta));

% Solving 
alpha0 = [-0.1, 2*pi/3];    % Initial guesses
max_iter = 20;              % number of iterations
tol = 1e-5;                 % set the tolerance 

% running the loop
for i = 1:length(beta)
    
    % initial guess alpha = -0.1
    [root, iter, flag] = Newton_Method(falpha, dfalpha, alpha0(1), tol, max_iter, beta(i), false, h_val);
    alpha(1, i) = root;
    iterations(1, i) = iter;
    success(1, i) = flag;
    
    % initial guess alpha = 2*pi/3
    [root, iter, flag] = Newton_Method(falpha, dfalpha, alpha0(2), tol, max_iter, beta(i), false, h_val);
    alpha(2, i) = root;
    iterations(2, i) = iter;
    success(2, i) = flag;

    % defining function again for fzero
    fun = @(x) (a1/a2)*cos(beta(i)) - (a1/a4)*cos(x) - cos(beta(i) - x) + C;
    
    % initial guess alpha = -0.1
    fzero_alpha(1, i) = fzero(fun, alpha0(1));
    
    % initial guess alpha = 2*pi/3
    fzero_alpha(2, i) = fzero(fun, alpha0(2));

    % Finite Difference with finite diff
    % initial guess alpha = -0.1
    [root, iter, flag] = Newton_Method(falpha, dfalpha, alpha0(1), tol, max_iter, beta(i), true, h_val);
    alpha_FD(1, i) = root;
    iterations_FD(1, i) = iter;
    success_FD(1, i) = flag;
    
    % initial guess alpha = 2*pi/3
    [root, iter, flag] = Newton_Method(falpha, dfalpha, alpha0(2), tol, max_iter, beta(i), true, h_val);
    alpha_FD(2, i) = root;
    iterations_FD(2, i) = iter;
    success_FD(2, i) = flag;   
end

% Errors for NS wrt analytical sol (-0.1 & 2*pi/3)
error_NS1 = abs(fzero_alpha(1,:) - alpha(1,:));         % alpha = -0.1
error_NS2 = abs(fzero_alpha(2,:) - alpha(2,:));         % alpha = 2*pi/3

error_NS_FD1 = abs(fzero_alpha(1,:) - alpha_FD(1,:));   % alpha = -0.1
error_NS_FD2 = abs(fzero_alpha(2,:) - alpha_FD(2,:));   % alpha = 2*pi/3

% Display results
disp(table(beta', alpha', iterations', success', fzero_alpha', alpha_FD', iterations_FD', success_FD', ...
     'VariableNames', {'beta', 'alpha_NS(-0.1 and 2pi/3)', 'Iterations_NS', 'Converged_NS', ...
     'alpha_fzero(-0.1 and 2pi/3)', 'alpha_FD(-0.1 and 2pi/3)', 'Iterations_NS_FD', 'Converged_NS_FD'}));

% Plotting NS vs fzero
figure(1);
hold on
grid on
plot(beta, alpha(1,:),'-o');
plot(beta, fzero_alpha(1,:),'-+');
plot(beta, alpha(2,:), '-square');
plot(beta, fzero_alpha(2,:),'-x');
title('Newton Solver vs fzero; validation');
xlabel('\beta (rad)', 'FontWeight', 'bold');
ylabel('\alpha (rad)', 'FontWeight', 'bold');
legend('Newton Solver; \alpha_0 = -0.1','fzero function; \alpha_0 = -0.1', 'Newton Solver; \alpha_0 = 2\pi/3', ...
    'fzero function; \alpha_0 = 2\pi/3', 'Location', 'northwest');
hold off

% Plotting NS vs NS FD
figure(2);
hold on
grid on
plot(beta, alpha(1,:),'-o');
plot(beta, alpha_FD(1,:),'-+');
plot(beta, alpha(2,:), '-square');
plot(beta, alpha_FD(2,:),'-x');
title('NS (Analytical) vs NS (Finite Diff; h = 1e-6)');
xlabel('\beta (rad)', 'FontWeight', 'bold');
ylabel('\alpha (rad)', 'FontWeight', 'bold');
legend('Newton Solver; \alpha_0 = -0.1','NS with FD; \alpha_0 = -0.1', 'Newton Solver; \alpha_0 = 2\pi/3', ...
    'NS with FD; \alpha_0 = 2\pi/3','Location', 'northwest');
hold off

% Plotting error in NS and NS_FD wrt analytical solution
figure(3);
hold on
grid on
plot(beta, error_NS1,'-o',"LineWidth", 2);              % alpha = -0.1 and NS
plot(beta, error_NS_FD1,'-+', "LineWidth", 1);          % alpha = -0.1 and NS_FD
plot(beta, error_NS2, '-square', "LineWidth", 2);       % alpha = 2*pi/3 and NS
plot(beta, error_NS_FD2,'-x', "LineWidth", 1);          % alpha = 2*pi/3 and NS_FD
title('NS and NS_{FD} vs fzero(analytical solution); validation');
xlabel('\beta (rad)', 'FontWeight', 'bold');
ylabel('Absolute Error (rad)', 'FontWeight', 'bold');
legend('NS; \alpha_0 = -0.1','NS_{FD}; \alpha_0 = -0.1', 'NS; \alpha_0 = 2\pi/3', ...
    'NS_{FD}; \alpha_0 = 2\pi/3', 'Location', 'best');
hold off

%% EX 2
clearvars; close all; clc; clear;

% Numerical Solution
% Define constants and initial values
v0 = 0;                 % [m/s]
m0 = 20;                % [kg]
c_m = 0.1;              % [kg/s]
f_t = 1;                % [N]
alpha = 0.01;           % [Ns/m]
rho = [0, 900];         % [kg/m^3]
Cd = 2.05;              % drag coefficient
Am = 1;                 % [m^2]

tspan = [0, 160];       % initial and final time [s]
h = [50, 20, 10, 1];    % step-sizes

N = 100;                            % number of CPU times
elap2 = zeros(N, length(h));        % to store the elapsed time RK2
elap4 = zeros(N, length(h));        % to store the elapsed time RK4
integ_err2 = zeros(size(h));        % to store max error RK2
integ_err4 = zeros(size(h));        % to store max error RK4

% Define the ODE equation, rho = 0
f_0 = @ (t, v) f_t/(m0 - c_m*t) - alpha * v/(m0 - c_m*t); % f = dv/dt 

% Analytical Solution
t_ana1 = tspan(1):1:tspan(2);   % corresponding to h = 1 
v_ana1 = f_t/alpha - (f_t/alpha - v0)*(1 - (c_m/m0).*t_ana1).^(alpha/c_m);

t_ana10 = tspan(1):10:tspan(2); % corresponding to h = 10
v_ana10 = f_t/alpha - (f_t/alpha - v0)*(1 - (c_m/m0).*t_ana10).^(alpha/c_m);

t_ana20 = tspan(1):20:tspan(2); % corresponding to h = 20
v_ana20 = f_t/alpha - (f_t/alpha - v0)*(1 - (c_m/m0).*t_ana20).^(alpha/c_m);

t_ana50 = tspan(1):50:tspan(2); % corresponding to h = 50
v_ana50 = f_t/alpha - (f_t/alpha - v0)*(1 - (c_m/m0).*t_ana50).^(alpha/c_m);

% Heun's Method
for i = 1:N                     % running through loop to store CPU time
    % for h = 50;
    tic;
    [v_sol50_RK2, t_sol50_RK2] = heuns_method(f_0, tspan, v0, h(1));
    elap2(i, 1) = toc;
    
    % for h = 20; 
    tic;
    [v_sol20_RK2, t_sol20_RK2] = heuns_method(f_0, tspan, v0, h(2));
    elap2(i, 2) = toc;
    
    % for h = 10; 
    tic;
    [v_sol10_RK2, t_sol10_RK2] = heuns_method(f_0, tspan, v0, h(3));
    elap2(i, 3) = toc;
    
    % for h = 1; 
    tic;
    [v_sol1_RK2, t_sol1_RK2] = heuns_method(f_0, tspan, v0, h(4));
    elap2(i, 4) = toc;
end

elap2_avg = mean(elap2);        % average values of elapsed time (RK2)

% Error plotting
% error at specific time points
integ_err2(1) = max(abs(v_ana50 - v_sol50_RK2));
integ_err2(2) = max(abs(v_ana20 - v_sol20_RK2));
integ_err2(3) = max(abs(v_ana10 - v_sol10_RK2));
integ_err2(4) = max(abs(v_ana1  - v_sol1_RK2 ));

% Interpolate to match exact time points
v_sol50_RK2_interp = interp1(t_sol50_RK2, v_sol50_RK2, t_ana1);
v_sol20_RK2_interp = interp1(t_sol20_RK2, v_sol20_RK2, t_ana1);
v_sol10_RK2_interp = interp1(t_sol10_RK2, v_sol10_RK2, t_ana1);
v_sol1_RK2_interp = interp1(t_sol1_RK2, v_sol1_RK2, t_ana1);

% Calculate errors
error_h50_RK2 = abs(v_ana1 - v_sol50_RK2_interp);
error_h20_RK2 = abs(v_ana1 - v_sol20_RK2_interp);
error_h10_RK2 = abs(v_ana1 - v_sol10_RK2_interp);
error_h1_RK2 = abs(v_ana1 - v_sol1_RK2_interp);

% Plot results RK2
figure(1);
hold on;
grid on;
plot(t_ana1, v_ana1, 'k-', "LineWidth", 1.5);
plot(t_sol1_RK2, v_sol1_RK2, 'square', "LineWidth", 1);
plot(t_sol10_RK2, v_sol10_RK2, 'b:pentagram', "LineWidth", 1);
plot(t_sol20_RK2, v_sol20_RK2, 'm:^', "LineWidth", 1);
plot(t_sol50_RK2, v_sol50_RK2,'g--o', "LineWidth", 1);
xlim([0 170])
title('Heuns Method (RK2) with various time steps vs Analytical solution');
xlabel('Time [s]', 'FontWeight', 'bold');
ylabel('Velocity [m/s]', 'FontWeight', 'bold');
legend('Analytical Solution','h = 1', 'h = 10', ...
    'h = 20', 'h = 50', 'Location', 'southeast');
hold off;

% Plot errors RK2
figure(2);
hold on; grid on;
plot(t_ana1, error_h50_RK2, "LineWidth", 1.5);
plot(t_ana1, error_h20_RK2, "LineWidth", 1.5);
plot(t_ana1, error_h10_RK2, "LineWidth", 1.5);
plot(t_ana1, error_h1_RK2, "LineWidth", 1.5);
xlim([0 170])
title('Heuns Method (RK2) error (h = 50, 20, 10 , 1) vs Analytical solution');
xlabel('Time [s]', 'FontWeight', 'bold');
ylabel('Velocity error [m/s]', 'FontWeight', 'bold');
legend('h = 50', 'h = 20', ...
    'h = 10', 'h = 1', 'Location', 'best');
hold off;

% CPU time vs max error RK2
figure(3);
semilogy(h, integ_err2, '-o', "LineWidth", 1.5);
hold on;
grid on;
semilogy(h, elap2_avg, '-o', "LineWidth", 1.5);
title('Error (max) & CPU time Trade off for RK2', 'FontWeight', 'bold');
xlabel('Step size (h) [s]', 'FontWeight', 'bold');
ylabel('Metric [dimensions in legend]', 'FontWeight', 'bold');
legend('Max Err [m/s]','CPU time [s]', 'Location','best');
hold off

% RK4
for i = 1:N
    % for h = 50;
    tic;
    [v_sol50_RK4, t_sol50_RK4] = RK4(f_0, tspan, v0, h(1));
    elap4(i, 1) = toc;
    
    % for h = 20; 
    tic;
    [v_sol20_RK4, t_sol20_RK4] = RK4(f_0, tspan, v0, h(2));
    elap4(i, 2) = toc;
    
    % for h = 10;
    tic;
    [v_sol10_RK4, t_sol10_RK4] = RK4(f_0, tspan, v0, h(3));
    elap4(i, 3) = toc;
    
    % for h = 1;
    tic;
    [v_sol1_RK4, t_sol1_RK4] = RK4(f_0, tspan, v0, h(4));
    elap4(i, 4) = toc;
end

elap4_avg = mean(elap4);        % average values of elapsed time (RK4)

% Error plotting
% error at specific time points
integ_err4(1) = max(abs(v_ana50 - v_sol50_RK4));
integ_err4(2) = max(abs(v_ana20 - v_sol20_RK4));
integ_err4(3) = max(abs(v_ana10 - v_sol10_RK4));
integ_err4(4) = max(abs(v_ana1 - v_sol1_RK4));

% Interpolate to match exact time points
v_sol50_RK4_interp = interp1(t_sol50_RK4, v_sol50_RK4, t_ana1);
v_sol20_RK4_interp = interp1(t_sol20_RK4, v_sol20_RK4, t_ana1);
v_sol10_RK4_interp = interp1(t_sol10_RK4, v_sol10_RK4, t_ana1);
v_sol1_RK4_interp = interp1(t_sol1_RK4, v_sol1_RK4, t_ana1);

% Calculate errors
error_h50_RK4 = abs(v_ana1 - v_sol50_RK4_interp);
error_h20_RK4 = abs(v_ana1 - v_sol20_RK4_interp);
error_h10_RK4 = abs(v_ana1 - v_sol10_RK4_interp);
error_h1_RK4 = abs(v_ana1 - v_sol1_RK4_interp);

% Plot results RK4
figure(5);
hold on
grid on
plot(t_ana1, v_ana1, 'k-', "LineWidth", 1.7);
plot(t_sol1_RK4, v_sol1_RK4, 'square', "LineWidth", 1);
plot(t_sol10_RK4, v_sol10_RK4, 'b:pentagram', "LineWidth", 1);
plot(t_sol20_RK4, v_sol20_RK4, 'm:^', "LineWidth", 1);
plot(t_sol50_RK4, v_sol50_RK4,'g--o', "LineWidth", 1);
xlim([0 170])
title('RK4 with various time step sizes vs Analytical solution');
xlabel('Time [s]', 'FontWeight', 'bold');
ylabel('Velocity [m/s]', 'FontWeight', 'bold');
legend('Analytical Solution','h = 1', 'h = 10', ...
    'h = 20', 'h = 50', 'Location', 'southeast');
hold off

% Plot errors RK4
figure(6);
hold on; grid on;
plot(t_ana1, error_h50_RK4, "LineWidth", 1.5);
plot(t_ana1, error_h20_RK4, "LineWidth", 1.5);
plot(t_ana1, error_h10_RK4, "LineWidth", 1.5);
plot(t_ana1, error_h1_RK4, "LineWidth", 1.5);
xlim([0 170])
title('RK4 error (h = 50, 20, 10 , 1) vs Analytical solution');
xlabel('Time [s]', 'FontWeight', 'bold');
ylabel('Velocity error [m/s]', 'FontWeight', 'bold');
legend('h = 50', 'h = 20', ...
    'h = 10', 'h = 1', 'Location', 'best');
hold off

% CPU time vs max error
figure(4);
semilogy(h, integ_err4, '-o', "LineWidth", 1.5);
hold on; grid on;
semilogy(h, elap4_avg, '-o', "LineWidth", 1.5);
title('Error (max) & CPU time Trade off for RK4', 'FontWeight', 'bold');
xlabel('Step size (h) [s]', 'FontWeight', 'bold');
ylabel('Metric [dimensions in legend]', 'FontWeight', 'bold');
legend('Max Err [m/s]','CPU time [s]', 'Location','best');
hold off

% Define the ODE equations, rho = 900
f_900 = @ (t, v) (-0.5*rho(2)*Cd*Am*v^2)/(m0 - c_m*t) + f_t/(m0 - c_m*t) - alpha * v/(m0 - c_m*t); % f = dv/dt

% RK2 with h = 1
[v_sol1_900, t_sol1_900] = heuns_method(f_900, tspan, v0, h(4));
%[v_sol1_900, t_sol1_900] = RK4(f_900, tspan, v0, 0.1); % RK4 is able to
%solve this for 0.1 step size

% Plot results
figure(7);
hold on
grid on
plot(t_sol1_900, v_sol1_900, '--*', "LineWidth", 1);
title('\rho = 900; RK2 with h = 1');
xlabel('Time [s]', 'FontWeight', 'bold');
ylabel('Velocity [m/s]', 'FontWeight', 'bold');
hold off

%Solving with an ODE
F = ode;
F.ODEFcn = f_900;
F.InitialValue = 0;
F.SelectedSolver
h = solve(F,0,160);

% Plot results
figure(8);
hold on
grid on
plot(h.Time, h.Solution, '--o', "LineWidth", 0.7);
title('ode45 for modified problem');
xlabel('Time [s]', 'FontWeight', 'bold');
ylabel('Velocity [m/s]', 'FontWeight', 'bold');
legend('ODE45', 'Location','southeast');
xlim([-10 170]);
hold off;
% ylim([-1e6 10]);



%% EX 3
% RK and BI methods
%% %%%%%%%%% Runge-Kutta (EX3) %%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all; clc; clear;

% alpha values
alpha = linspace(0, pi, 100);

% Define parameters and initial values
h_02 = 1.2;
h_04 = [1.6 2.5];               % two different initial guesses
h2 = zeros(size(alpha));        % for RK2
h4 = zeros(2, length(alpha));   % step size corresponding to different initial values
eigA = zeros(2, length(alpha)); % to store eigenvalues of A matrix

% solving for alpha = [0, pi]
for i = 1:length(alpha)

    % Matrix A
    A = [0, 1; -1, 2*cos(alpha(i))];
    eigA(:, i) = eig(A);
    
    % RK2 
    eig_RK2 = @(h) max(abs(eig((eye(2) + h.*A + 0.5*(h^2.*A^2))))) - 1; % eigenvalue of F_RK2
    try
        h2(i) = fzero(eig_RK2, h_02);
    catch
        h2(i) = NaN;
    end
    
    % RK4 
    eig_RK4 = @(h) max(abs(eig((eye(2) + h.*A + 0.5*(h^2.*A^2) + ...
        (1/6)*(h^3.*A^3) + (1/24)*(h^4.*A^4))))) - 1;                   % eigenvalue of F_RK4
    try
        h4(1, i) = fzero(eig_RK4, h_04(1));
        h4(2, i) = fzero(eig_RK4, h_04(2));
    catch
        h4(1, i) = NaN;
        h4(2, i) = NaN;
    end
end

% Plot the results 

fprintf('The value of h for RK2 and alpha = pi is: %.5f\n', h2(100));
fprintf('The value of h for RK2 and alpha = pi is: %.5f\n', h4(2, 100));

% %%%% h vs alpha %%%%%
% figure(1);
% hold on;
% plot(alpha, h2,'c.', 'LineWidth', 1.5);
% plot(alpha, h4(1,:), 'r.', 'LineWidth', 1.5);
% plot(alpha, h4(2,:), 'b.', 'LineWidth', 1.5);
% xlabel('\alpha (rad)');
% ylabel('h');
% grid on;
% hold off;
% axis equal;

%%%%%% h-lambda plane %%%%%%%%
% for plotting purpose
real2 = real(eigA.*h2);
imag2 = imag(eigA.*h2);
real41 = real(eigA.*h4(1, :));
imag41 = imag(eigA.*h4(1, :));
real42 = real(eigA.*h4(2, :));
imag42 = imag(eigA.*h4(2, :));

figure(2);
hold on;
% plot(real(eigA.*h2), imag(eigA.*h2), 'ro', 'LineWidth', 1.5);
% plot(real(eigA.*h4(1,:)), imag(eigA.*h4(1,:)), 'co', 'LineWidth', 1.5); % for guess-1
% plot(real(eigA.*h4(2,:)), imag(eigA.*h4(2,:)), 'co', 'LineWidth', 1.5); % for guess-2
plot(real2(1,:), imag2(1,:), 'ro', 'LineWidth', 2);     % upper half
plot(real2(2,:), imag2(2,:), 'ro', 'LineWidth', 2);     % lower half
plot(real41(1,:), imag41(1,:), 'co', 'LineWidth', 2);   % for guess-1
plot(real41(2,:), imag41(2,:), 'co', 'LineWidth', 2);
plot(real42(1,:), imag42(1,:), 'co', 'LineWidth', 2);   % for guess-2
plot(real42(2,:), imag42(2,:), 'co', 'LineWidth', 2);
xlabel('Re(h\lambda)');
ylabel('Im(h\lambda)');
title('Stability domain of RK2 and RK4 on (h\lambda)-plane');
legend('RK2','', 'RK4', 'Location', 'best');
grid on;
hold off;
axis equal;

%% %%%%%%%%% Backward Interpolation %%%%%%%%%%%%%%%%%%%
clearvars; close all; clc; clear;

% alpha values
alpha = linspace(0, pi, 100);

% Define parameters
theta = [0.2 0.3 0.4 0.6 0.8];        % values of theta
h0 = [5 5 6 6 5];                     % Initial values of h

% initialize vectors
h = zeros(5, length(alpha));
eigA = zeros(2, length(alpha));

% solving for alpha = [0, pi]
for i = 1:length(alpha)

    % Matrix A
    A = [0, 1; -1, 2*cos(alpha(i))];
    eigA(:, i) = eig(A);
    
    for j = 1:length(theta)

        thetaF = theta(j);     % for RK2 step (forward)
        thetaB = 1 - theta(j); % for BRK2 step (backward)

        eigF = @(h) max(abs(eig((eye(2) - thetaB*h.*A + 0.5*((thetaB*h)^2.*A^2))\ ...
        (eye(2) + thetaF*h.*A + 0.5*((thetaF*h)^2.*A^2))))) - 1; % eigenvalues of F_BI2_theta
        
        try
            h(j, i) = fzero(eigF, h0(j));
        catch
            h(j, i) = NaN;
        end
    end
end

% Plot the results 
%fprintf('The value of h for RK2 and alpha = pi is: %.5f\n', h(1, 250));

% %%%% h vs alpha %%%%%
% figure(1);
% hold on;
% plot(alpha, h ,'r.', 'LineWidth', 1.5);
% %plot(alpha, h(2, :) ,'r.', 'LineWidth', 1.5);
% xlabel('\alpha (rad)');
% ylabel('h');
% grid on;
% hold off;
% axis equal;

%%%%%% h-lambda plane %%%%%%%%
figure(2);
hold on;
%plot(real(eigA.*h), imag(eigA.*h), 'r.', 'LineWidth', 1.5);
plot(real(eigA.*h(2, :)), imag(eigA.*h(2, :)), 'ro', 'LineWidth', 1.5);
xlabel('Re(h\lambda)');
ylabel('Im(h\lambda)');
title('Stability domain of BI2 for \theta = 0.3 in (h\lambda)-plane');
legend('\theta = 0.3');
grid on;
hold off;
axis equal;

% for plotting purposes
real1 = real(eigA.*h(1, :));
imag1 = imag(eigA.*h(1, :));
real3 = real(eigA.*h(3, :));
imag3 = imag(eigA.*h(3, :));
real4 = real(eigA.*h(4, :));
imag4 = imag(eigA.*h(4, :));
real5 = real(eigA.*h(5, :));
imag5 = imag(eigA.*h(5, :));

figure(3);
hold on;
%plot(real(eigA.*h), imag(eigA.*h), 'r.', 'LineWidth', 1.5);
plot(real1(1,:), imag1(1,:), 'ro', 'LineWidth', 2); % upper half 
plot(real1(2,:), imag1(2,:), 'ro', 'LineWidth', 2); % lower half
plot(real3(1,:), imag3(1,:), 'co', 'LineWidth', 2);
plot(real3(2,:), imag3(2,:), 'co', 'LineWidth', 2);
plot(real4(1,:), imag4(1,:), 'bo', 'LineWidth', 2);
plot(real4(2,:), imag4(2,:), 'bo', 'LineWidth', 2);
plot(real5(1,:), imag5(1,:), 'go', 'LineWidth', 2);
plot(real5(2,:), imag5(2,:), 'go', 'LineWidth', 2);
xlabel('Re(h\lambda)');
ylabel('Im(h\lambda)');
title('Stability domain of BI2 for \theta = [0.2, 0.4, 0.6, 0.8] in (h\lambda)-plane');
legend('\theta = 0.2', '', '\theta = 0.4', '', '\theta = 0.6', '', '\theta = 0.8', 'Location', 'best');
grid on;
hold off;
axis equal;

%% EX 4
clearvars; close all; clc; clear;

% Define constants
K_c = 0.0042;   % Convective coefficient [J/(sK)]
K_r = 6.15e-11; % Radiation coefficient [J/sK^4]
T_a = 277;      % Air temperature [K]
C = 45;         % Mass thermal capicity [J/K]
T0 = 555;       % Initial mass temperature [K]
t0 = 0;         % Initial time [s]
tf = 36000;     % Final time [s] (arbitrarily selected; simple multiple of 720 and 1440)
h = [720 1440]; % step sizes for RK2 and RK4, recpectively

% mathematical model of the problem
f = @ (t, T) (-1/C) * (K_c * (T - T_a) + K_r * (T^4 - T_a^4));

% using general purpose RK2 (Heun's)
[T_RK2, tspan] = heuns_method(f, [t0 tf], T0, h(1));

figure(1);
hold on
grid on
plot(tspan, T_RK2, '-o', "LineWidth",1);
xlabel('Time');
ylabel('Temperature (K)');
title('RK2 Solution');
hold off

% Using general pupose RK4
[T_RK4, tspan_RK4] = RK4(f, [t0 tf], T0, h(2));

figure(2);
hold on
grid on
plot(tspan_RK4, T_RK4, '-o', "LineWidth",1);
xlabel('Time');
ylabel('Temperature (K)');
title('RK4 Solution');
hold off

% Solve with ODE45 for exact solution
[t_exact, T_exact] = ode45(f, [t0 tf], T0);

figure(3);
hold on
grid on
plot(tspan, T_RK2,'-o', "LineWidth",1);
plot(tspan_RK4, T_RK4,'-^', "LineWidth",1);
plot(t_exact, T_exact,'-square', "LineWidth",1);
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Temperature (K)', 'FontWeight', 'bold');
title('RK2(h = 720) vs RK4(h = 1440) vs ODE45(Exact Solution)');
legend('RK2 (h=720)', 'RK4(h = 1440)', 'ode45 (exact)', 'Location', 'best');
hold off

% Interpolate RK2 and RK4 to match exact time points
T_RK2_interp = interp1(tspan, T_RK2, t_exact);
T_RK4_interp = interp1(tspan_RK4, T_RK4, t_exact);

% Calculate errors
error_RK2 = abs(T_exact - T_RK2_interp);
error_RK4 = abs(T_exact - T_RK4_interp);

% Plot errors
figure(4);
plot(t_exact, error_RK2, '-o', 'LineWidth', 1); 
hold on;
plot(t_exact, error_RK4, '-^', 'LineWidth', 1);
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Absolute Error (K)', 'FontWeight', 'bold');
title('Error Analysis for RK2 and RK4', 'FontWeight', 'bold');
legend('RK2 Error', 'RK4 Error', 'Location', 'best');
grid on;

%% EX 5
clearvars; close all; clc; clear;

% Given circuit parameters
R = 25;      % Resistance (Ohms)
L = 20e-3;   % Inductance (H)
C = 200e-3;  % Capacitance (F)
v0 = 12;     % Volts (V)

% System matrix
A = [0, 1; -1/(L*C), -R/L];

x0 = [v0; 0];                       % Initial condition: charge on capacitor and current in circuit
h_IEX4 = 0.1;                       % Step size for IEX4
t_end = 30;                         % End time >= 5*tau = 5*R*C
t_IEX4 = 0:h_IEX4:t_end;            % Time vector
x_IEX4 = zeros(2, length(t_IEX4));  % State variable storage
x_IEX4(:, 1) = x0;                  % Intial vector

% Compute eigenvalues
eigenvalues2 = eig(A);

% Plot eigenvalues in the complex plane
% figure(1);
% plot(real(eigenvalues), imag(eigenvalues), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
% xlabel('Real Part');
% ylabel('Imaginary Part');
% title('Eigenvalues of the System');
% grid on;

% IEX4 algorithm (10 fcn evaluation)
for i = 1:length(t_IEX4)-1
    % 1st predictor
    F1 = @(x_next) x_IEX4(:, i) + h_IEX4 * (A * x_next) - x_next; % Implicit equation
    k1 = fsolve(F1, x_IEX4(:, 1)); % Solve for k1

    % 2nd predictor
    F2a = @(x_next) x_IEX4(:, i) + (h_IEX4/2) * (A * x_next) - x_next; % Implicit equation
    k2a = fsolve(F2a, x_IEX4(:, 1)); 

    F2 = @(x_next) k2a + (h_IEX4/2) * (A * x_next) - x_next; % Implicit equation
    k2 = fsolve(F2, x_IEX4(:, 1)); % solve for k2

    % 3rd predictor
    F3a = @(x_next) x_IEX4(:, i) + (h_IEX4/3) * (A * x_next) - x_next; % Implicit equation
    k3a = fsolve(F3a, x_IEX4(:, 1));

    F3b = @(x_next) k3a + (h_IEX4/3) * (A * x_next) - x_next; % Implicit equation
    k3b = fsolve(F3b, x_IEX4(:, 1));

    F3 = @(x_next) k3b + (h_IEX4/3) * (A * x_next) - x_next; % Implicit equation
    k3 = fsolve(F3, x_IEX4(:, 1)); % Solve for k3

    % 4th predictor
    F4a = @(x_next) x_IEX4(:, i) + (h_IEX4/4) * (A * x_next) - x_next; % Implicit equation
    k4a = fsolve(F4a, x_IEX4(:, 1));

    F4b = @(x_next) k4a + (h_IEX4/4) * (A * x_next) - x_next; % Implicit equation
    k4b = fsolve(F4b, x_IEX4(:, 1));

    F4c = @(x_next) k4b + (h_IEX4/4) * (A * x_next) - x_next; % Implicit equation
    k4c = fsolve(F4c, x_IEX4(:, 1));

    F4 = @(x_next) k4c + (h_IEX4/4) * (A * x_next) - x_next; % Implicit equation
    k4 = fsolve(F4, x_IEX4(:, 1)); % Solve for k4

    % Corrector
    x_IEX4(:, i+1) = (-1/6) * k1 + 4 * k2 - (27/2) * k3 + (32/3) * k4; 
end

% RK2 implementation
h_RK2 = 2/1400;                     % step size for RK2 
t_RK2 = 0:h_RK2:t_end;              % Time vector
x_RK2 = zeros(2, length(t_RK2));    % State variable storage
x_RK2(:, 1) = x0;                   % Intial vector

% Solving RK2
for i = 1:length(t_RK2)-1
    % Predictor step
    xp = x_RK2(:, i) + h_RK2 * A * x_RK2(:, i);

    % Corrector step
    x_RK2(:, i+1) = x_RK2(:, i) + (h_RK2/2) * (A * x_RK2(:, i) + A * xp);
end

% Plotting eigenvalues in RK2 and IEX4 stability domain

% for RK2
% Grid in the complex plane
[x1, y1] = meshgrid(-4:0.1:2, -3:0.1:3);
z1 = x1 + 1i*y1;
% Stability function for RK2
R_RK2 = 1 + z1 + z1.^2 / 2;
% Magnitude of R(z)
R_mag_RK2 = abs(R_RK2);

figure(1);
contour(x1, y1, R_mag_RK2, [1, 1], 'LineWidth', 1.5); hold on; % RK2
plot(real(eigenvalues2*h_RK2), imag(eigenvalues2*h_RK2), 'rx', 'MarkerSize', 10, 'LineWidth', 1.5);
xlabel('Re(h\lambda)');
ylabel('Im(h\lambda)');
title('(h\lambda)-plane for RK2');
legend('Stability domain of RK2 inside curve', 'Eigenvalues', 'Location', 'northeast');
grid on;
axis equal;

% for IEX4
% Grid in the complex plane
[x2, y2] = meshgrid(-5:0.1:15, -10:0.1:10);
z2 = x2 + 1i*y2;
% Stability function for IEX4
R_IEX4 = z2-12.5/2;
% Magnitude of R(z)
R_mag_IEX4 = abs(R_IEX4);

figure(2);
contour(x2, y2, R_mag_IEX4, [12.5/2, 12.5/2], 'LineWidth', 1.5); hold on; % IEX4
plot(real(eigenvalues2*0.01), imag(eigenvalues2*0.01), 'rx', 'MarkerSize', 10, 'LineWidth', 1.5);
xlabel('Re(h\lambda)');
ylabel('Im(h\lambda)');
title('(h\lambda)-plane for IEX4');
legend('Stability domain of IEX4 outside curve', 'Eigenvalues', 'Location', 'northeast');
grid on;
axis equal;

% Plotting IEX4 results
% Charge vector
figure(3);
plot(t_RK2, x_RK2(1, :), '-', 'LineWidth', 1.5); hold on;
plot(t_IEX4, x_IEX4(1, :), '-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Charge (Q)');
legend('RK2', 'IEX4');
title('IEX4 (10 Function Evaluations) and RK2 for electrical Circuit');
grid on;

% Current vector (On capacitor)
figure(4);
plot(t_RK2, x_RK2(2, :), '-', 'LineWidth', 1.5); hold on;
plot(t_IEX4, x_IEX4(2, :), '-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Current (Amp)');
legend('RK2', 'IEX4');
title('IEX4 (10 Function Evaluations) and RK2 for electrical Circuit');
grid on;

%Errors
% Interpolate RK2 and RK4 to match exact time points
current_RK2_interp = interp1(t_RK2, x_RK2(2, :), t_IEX4);

% Calculate errors
error_RK2 = abs(x_IEX4(2, :) - current_RK2_interp);

% Plot errors
figure(5);
plot(t_IEX4, error_RK2, '-', 'LineWidth', 1.5); 
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Absolute Error (K)', 'FontWeight', 'bold');
title('Error Analysis for RK2', 'FontWeight', 'bold');
legend('RK2 Error', 'Location', 'best');
grid on;

%% EX 6
clearvars; close all; clc; clear;

% Define parameters
k = 0.9;                % attenuation factor
rho = [1.225, 15, 60];  % Density vector [kg/m3]
m_b = 1;                % Ball mass [kg]
A_b = 0.07;             % Ball area [m2]
V_b = 0.014;            % Ball volume [m3]
C_d = 1.17;             % Drag coeff
g = 9.81;               % Gravitional coefficient [m/s2]

% Initialize storage for results
t_all = [];
x_all = [];
v_all = [];

% Initial condition
x0 = 10;                % Initial height (m)
v0 = 0;                 % Initial velocity (m/s)

% Simulation loop for multiple bounces
t_start = 0;            % Start time
t_final = 10;           % Final time
x = x0;                 % Initial position
v = v0;                 % Initial velocity

% Define ODE systems
dxdt_1 = @(t, X) [X(2); -g + (-0.5*rho(1)*C_d*A_b*X(2)*abs(X(2)))/m_b + (rho(1)*V_b*g)/m_b]; % for rho = 1.225
dxdt_2 = @(t, X) [X(2); -g + (-0.5*rho(2)*C_d*A_b*X(2)*abs(X(2)))/m_b + (rho(2)*V_b*g)/m_b]; % for rho = 15
dxdt_3 = @(t, X) [X(2); -g + (-0.5*rho(3)*C_d*A_b*X(2)*abs(X(2)))/m_b + (rho(3)*V_b*g)/m_b]; % for rho = 60

% Looooooping for every event detection
while t_start < t_final
    
    % Options for event detection
    options = odeset('Events', @(t, X) bounceBall(t, X, k));

    % Solve using ode45
    [t, X, te, Xe, ie] = ode45(dxdt_1, [t_start, t_final], [x; v], options); % edit dxdt_# for different fluid densities

    % Store results
    t_all = [t_all; t]; 
    x_all = [x_all; X(:, 1)]; 
    v_all = [v_all; X(:, 2)]; 
    
    % Check if the event (ground contact) occurred
    if isempty(te)
        break;              % No further events detected
    end
    
    % Update initial conditions for the next bounce
    t_start = te;           % Update start time
    x = 0;                  % Ball position resets to the ground
    v = -k * Xe(2);         % Reverse and attenuate velocity
end

% Plotting results
figure(1);
plot(t_all, x_all, '-', 'LineWidth', 1.5); 
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Position (m)', 'FontWeight', 'bold');
title('Bouncing ball (ode45) position for \rho = 1.225 kg/m^3', 'FontWeight', 'bold');
grid on;

figure(2);
plot(t_all, v_all, '-', 'LineWidth', 1.5); 
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Velocity (m)', 'FontWeight', 'bold');
title('Bouncing ball (ode45) velocity for \rho = 1.225 kg/m^3', 'FontWeight', 'bold');
grid on;

% ------------------------------------------------------
% iMPLETEMENT rk4
%-------------------------------------------------------

h = 0.65;                       % Step size (s)

% Initialize variables
time = 0:h:t_final;             % Time vector
pos = zeros(size(time));        % Height array
vel = zeros(size(time));        % Velocity array
pos(1) = x0;
vel(1) = v0;

% Define the system of equations
f1 = @(x, v) v;                 % dx/dt = v
f2 = @(x, v) -g + (-0.5*rho(3)*C_d*A_b.*v.*abs(v))/m_b + (rho(3)*V_b*g)/m_b; % dv/dt

% RK4 loop
for i = 1:length(time)-1

    
    % % Calculate abbreviates k1, k2, k3, and k4 for position and velocity
    k1_pos = f1(pos(i), vel(i));
    k1_vel = f2(pos(i), vel(i));
    
    k2_pos = f1(pos(i) + h/2 * k1_pos, vel(i) + h/2 * k1_vel);
    k2_vel = f2(pos(i) + h/2 * k1_pos, vel(i) + h/2 * k1_vel);
    
    k3_pos = f1(pos(i) + h/2 * k2_pos, vel(i) + h/2 * k2_vel);
    k3_vel = f2(pos(i) + h/2 * k2_pos, vel(i) + h/2 * k2_vel);
    
    k4_pos = f1(pos(i) + h * k3_pos, vel(i) + h * k3_vel);
    k4_vel = f2(pos(i) + h * k3_pos, vel(i) + h * k3_vel);
    
    % Update position and velocity
    pos(i+1) = pos(i) + h/6 * (k1_pos + 2*k2_pos + 2*k3_pos + k4_pos);
    vel(i+1) = vel(i) + h/6 * (k1_vel + 2*k2_vel + 2*k3_vel + k4_vel);
    
    % Check for bounce
    if pos(i+1) < 0
        pos(i+1) = 0;               % Reset position to ground
        vel(i+1) = -k * vel(i+1);   % Reverse and attenuate velocity
    end
end

% Plotting results
figure(3);
plot(time, pos, '-', 'LineWidth', 1.5); 
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Position (m)', 'FontWeight', 'bold');
title('Bouncing ball (RK4) position for \rho = 60 kg/m^3', 'FontWeight', 'bold');
grid on;

figure(4);
plot(time, vel, '-', 'LineWidth', 1.5); 
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Velocity (m)', 'FontWeight', 'bold');
title('Bouncing ball (RK4) velocity for \rho = 60 kg/m^3', 'FontWeight', 'bold');
grid on;

%% Functions

% Implementation of general purpose Newton Solver
function [root, iter, success] = Newton_Method(f, df, x0, tol, max_iter, beta, finite_diff, h)
    
    % Inputs
    % f: the equation to be solved
    % df: analytical derivative of equation
    % x0: initial value
    % tol: accepted tolerance
    % max_iter: number of maximum iterations allowed
    % beta: input angle
    % finite_diff: condition for finite difference method
    % h: step-size in FD method to estimate derivative of equation
    
    % Outputs
    % root: soltion to the equation
    % iter: number of iterations taken
    % success: condition to check if solution converged in max_iter

    % Initialize
    success = false;
    x = x0;
    b = beta;

    % loop for iterations
    for i = 1:max_iter
        fx = f(x,b);
        if ~finite_diff
            % Analytical f'(x)
            dfx = df(x,b);
        else
            % Finite Difference (Central) f'(x)
            dfx = (f(x + h, beta) - f(x - h, beta)) / (2 * h);
        end

    % Check if derivative is zero (avoid division by zero)
        if abs(dfx) < eps
            warning('Derivative near zero; stopping iterations.');
            return;
        end

    % compute next guess
        x_new = x - fx/dfx;

    % when converge
        if abs(fx) < tol
            success = true;
            break;
        end

    % update variable
        x = x_new;
    end

    % Assigning final value
    root = x;
    iter = i;

    % Check if maximum iterations were reached
    if ~success
        warning('Newton solver did not converge within the maximum number of iterations.');
    end
end

% Implementation of general purpose Heun's Method (RK2)
function [x, t] = heuns_method(f, tspan, x0, h)
    
    % Inputs
    % f: f(x,t) = dx/dt
    % tspan: time span for which the problem is required to be solved [t0, tf]
    % x0: value of x at t0
    % h: step-size

    % Outputs
    % x: solution vector
    % t: time vector (distributed based on step size)

    % Initialize
    t = tspan(1):h:tspan(2);    % time vector
    x = zeros(1, length(t));    % Solution vector
    x(1) = x0;                  % defining initial condition

    % Iterations
    for i = 1:length(t)-1
        % Predictor step
        xp = x(i) + h * f(t(i), x(i));

        % Corrector step
        x(i+1) = x(i) + (h/2) * (f(t(i), x(i)) + f(t(i+1), xp));
    end
end

% Implementation of general purpose RK4
function [x, t] = RK4(f, tspan, x0, h)

    % Inputs - same as Heusn's method
    % Outputs - Solution and time vectors

    % Initialize
    t = tspan(1):h:tspan(2);    % time vector
    x = zeros(1, length(t));    % Solution vector
    x(1) = x0;                  % defining initial condition

    % Iterations
    for i = 1:length(t) - 1
        % Calculate abbreviates k1, k2, k3, and k4
        k1 = f(t(i), x(i));
        k2 = f(t(i) + h/2, x(i) + (h/2)*k1);
        k3 = f(t(i) + h/2, x(i) + (h/2)*k2);
        k4 = f(t(i) + h, x(i) + h*k3);
    
    % update the solution
        x(i+1) = x(i) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end

% Bouncing ball event detection function in odeset
function [value, isterminal, direction] = bounceBall(t, x, k)
    value = x(1);       % Check ball height = 0
    isterminal = 1;     % Stop the integration
    direction = -1;     % Detect downward crossing
    if value < 0
        x(2) = -k*x(2); % Update velocity after impact
    end
end