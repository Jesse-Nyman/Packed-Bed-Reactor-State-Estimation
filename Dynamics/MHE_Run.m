% Full MHE and Simulation Code in One File

clear all; clc; close all
addpath('casadi-3.6.6\')
import casadi.*

%% Model parameters
k0 = 3.16 * 10^-2;      % (kmol/(kg cat * s * Pa))   Pre-exponential factor for hydrogenation
K0 = 3.16 * 10^-13;     % (1/Pa)                     Adsorption constant for benzene
Q = -6.89 * 10^7;       % (kJ/kmol)                  Benzene adsorption activation energy
E = 5.76 * 10^7;        % (kJ/kmol)                  Activation energy of hydrogenation
rhoCp = 7.322 * 10^5;   % (J/(m^3*C))                Effective heat capacity, of reactor
negDeltaH = 2.09 * 10^8;% (J/kmol)                   Heat of reaction
rhoB = 4.14 * 10^2;     % (kg/m^3)                   Catalyst bulk density
epz = 0.58;             %                            Void fraction
MT = 4.05 * 10^-4;      % (kg/kmol)                  Poison adsorption capacity of catalyst
k0d = 1.80 * 10^-4;     % ((Pa*s)^-1)                Pre-exponential factor for poisoning
Ed = 4.52 * 10^6;       % (J/kmol)                   Activation energy of poisoning
P = 1.01 * 10^5;        % (Pa or N/m^2)              Pressure
L = 0.118 * 5;          % (m)                        Catalyst's bed length
r = 0.0078;             % (m)                        Reactor's radius
DeB = 4.5 * 10^-5;      % Diffusion coefficient
nstages = 100;
sensorNumber = 5;
sensorPosition = round(linspace(20, 80, sensorNumber));
nstates = 6;  % Number of states per stage (e.g., XB, XH, XP, XT, T, Activity)
global theta;
theta = {nstages, k0, K0, Q, E, rhoCp, negDeltaH, rhoB, epz, MT, k0d, Ed, P, L, r, DeB, sensorPosition};
A = pi * r^2;           % (m^2)                      Reactor's cross-section area

% Initial conditions for state variables
xBzt = 1.5E-2 * ones(nstages, 1); % Initial concentration for XB
xPzt = zeros(nstages, 1);         % Initial concentration for XP
xTzt = zeros(nstages, 1);         % Initial concentration for XT
xHzt = 1.0 - xBzt - xPzt;
Tzt = 350 * ones(nstages, 1);     % Initial temperature
Azt = zeros(nstages, 1);          % Initial catalyst activity
Azt(20:80) = 1;                   % Catalyst bed position
Xz0 = [xBzt; xHzt; xPzt; xTzt; Tzt; Azt];

%% Inputs
Feed = 2.5 * 10^-5;   % (m^3 / sec) Volumetric flowrate
ut = Feed / (A * epz) * ones(1, 100); % Linear velocity
xB0t = xBzt(1) * ones(1, 100); % Feed concentration XB
xP0t = xPzt(1) * ones(1, 100); % Feed concentration XP
xT0t = 6.36E-4 + xTzt(1) * ones(1, 100); % Feed concentration XT
xH0t = 1 - xB0t - xT0t;
T0t = Tzt(1) * ones(1, 100); % Feed temperature

% Combine inputs
Ut = [xB0t', xH0t', xP0t', xT0t', T0t', ut'];

% Integrator
stepSize = 1; % Simulation step size in seconds
F = buildIntegrator(@modelCD, [nstates * nstages, 6], stepSize);

% Simulate the system to obtain yout and xout
[yout, ~, xout] = simulateSystem(F, @modelCD_out, Xz0, Ut);

%% MHE Parameters
nx = nstates * nstages; % State dimension
ny = length(sensorPosition); % Number of sensors
nu = 6; % Control dimension
nw = 5; % Disturbance dimension
N_mhe = 20; % Horizon length for MHE

% Initialize variables
x0_mhe = zeros(nx, 1);
u_mhe = zeros(nu, N_mhe);
w_mhe = zeros(nw, N_mhe);
Q = diag(ones(nx, 1) * 1e-4);
R = diag(ones(ny, 1) * 1e-2);
x_est = x0_mhe;

%% MHE Setup
x = MX.sym('x', nx, N_mhe + 1);
u = MX.sym('u', nu, N_mhe);
w = MX.sym('w', nw, N_mhe);
y_meas = MX.sym('y_meas', ny, N_mhe);

cost = 0;
g = [];

for k = 1:N_mhe
    x_next = F('x0', x(:, k), 'p', [u(:, k); w(:, k)]);
    x_pred = x_next.x_next;
    g = [g; x(:, k+1) - x_pred];
    y_pred = modelCD_out(0, x(:, k+1), theta);
    cost = cost + (y_meas(:, k) - y_pred)' * R * (y_meas(:, k) - y_pred);
end
x0_est = MX.sym('x0_est', nx);
cost = cost + (x(:, 1) - x0_est)' * Q * (x(:, 1) - x0_est);

opts = struct('ipopt', struct('print_level', 0, 'tol', 1e-6));
nlp = struct('x', [x(:); u(:); w(:)], 'f', cost, 'g', g);
solver = nlpsol('solver', 'ipopt', nlp, opts);

%% MHE Simulation Loop
for t = 1:100 - N_mhe
    y_current = yout(:, t:t+N_mhe-1);
    sol = solver('x0', [x_est(:); u_mhe(:); w_mhe(:)], 'p', y_current, ...
                 'lbg', zeros(nx * N_mhe, 1), 'ubg', zeros(nx * N_mhe, 1), ...
                 'x0_est', x_est);
    x_est = reshape(full(sol.x(1:nx)), nx, 1);
    x_estimates(:, t) = x_est;
end

% Plot the results
figure;
for i = 1:nstates
    subplot(3, 2, i);
    plot(1:100 - N_mhe, x_estimates((i-1)*nstages+1:i*nstages, :));
    title(['Estimated State ' num2str(i)]);
end

%% Define Functions

function F = buildIntegrator(f, d, stepSize)
    import casadi.*
    x = MX.sym('x', d(1));
    u = MX.sym('u', d(2));
    DAE = struct('x', x, 'p', u, 'ode', f(0, x, u, theta));
    options = struct('tf', stepSize);
    F = integrator('F', 'cvodes', DAE, options);
end

function xdot = modelCD(~, x, u, theta)
    % Insert the full modelCD code here
    % Use the theta parameters and compute xdot based on your equations
    xdot = zeros(length(x), 1); % Placeholder, replace with actual equations
end

function y = modelCD_out(~, x, theta)
    % Insert the full modelCD_out code here
    % Simulate measurements based on sensor positions
    nstages = theta{1};
    sensorPosition = theta{17};
    T = x(4*nstages+1 : 5*nstages);
    y = T(sensorPosition);
end

function [yout, tout, xout] = simulateSystem(F, g, X0, Ut)
    nSteps = size(Ut, 1);
    xout = zeros(length(X0), nSteps);
    yout = zeros(length(g(X0)), nSteps);
    x = X0;
    for i = 1:nSteps
        xout(:, i) = x;
        yout(:, i) = g(0, x, theta);
        res = F('x0', x, 'p', Ut(i, :)');
        x = full(res.xf);
    end
    tout = (1:nSteps) * stepSize;
end
