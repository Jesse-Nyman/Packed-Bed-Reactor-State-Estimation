%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% catalystDecay_simulate.m

% After the nonlinear simulation:
% - At the start of each segment, linearize around the nonlinear state.
% - simulate the linear model forward for that segment using the 
%   newly obtained linearization.
% This creates a linear time-varying approximation that stays closer 
% to the nonlinear trajectory over the full simulation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all

%% ADD OWN PATH
addpath('C:\Users\Jesse_N\OneDrive - Aalto University\Documents\Masters Thesis\starting code\casadi-3.6.6')
import casadi.*; 

%% Model parameters (check dynamic sim for comments)
k0        = 3.16e-2;       
K0        = 3.16e-13;      
Q         = -6.89e7;       
E         = 5.76e7;        
rhoCp     = 7.322e5;       
negDeltaH = 2.09e8;        
rhoB      = 4.14e2;        
epz       = 0.58;          
MT        = 4.05e-4;       
k0d       = 1.80e-4;       
Ed        = 4.52e6;        
P         = 1.01e5;        
L         = 0.118*5;       
r         = 0.0078;        
DeB       = 4.5e-5;        

% Spatial discretization (vary as wanted)
nstages   = 100;
front = 0.2;
back  = 0.8;
anteCatalyst = 1:floor(front*nstages);
postCatalyst = floor(back*nstages)+1:nstages;
catalystBed  = anteCatalyst(end)+1 : postCatalyst(1)-1;
sensorNumber   = 5;
sensorPosition = round(linspace(catalystBed(1),catalystBed(end),sensorNumber));

theta     = {nstages, k0, K0, Q, E, rhoCp, negDeltaH, rhoB, epz, MT, k0d, Ed, P, L, r, DeB, sensorPosition};
A         = pi*r^2;  

%% Initial conditions
xBzt = 1.5E-2 * ones(nstages,1);
xPzt = 0.0     * ones(nstages,1);
xTzt = 0.0     * ones(nstages,1);
xHzt = 1.0 - xBzt - xPzt - xTzt; 
Tzt  = 350 * ones(nstages,1);
Azt  = zeros(nstages,1);
Azt(catalystBed)=1;
Xz0 = [xBzt; xHzt; xPzt; xTzt; Tzt; Azt];

% Inputs (constant)
Feed = 2.5e-5;
epz = 0.58; 
ut   = (Feed/(A*epz));
xB0 = xBzt(1);
xP0 = xPzt(1);
xT0 = 6.36e-4 + xTzt(1);
xH0 = 1 - xB0 - xT0 - xP0;
T0_init  = Tzt(1);
U = [xB0; xH0; xP0; xT0; T0_init; ut]; % 6x1

f = @(x,u) modelCD(0, x, u, theta);
g = @(x)   modelCD_out(0, x, theta);

nstates = 6; 
Nx = nstates*nstages;
Ny = size(g(Xz0),1);
Nu = 6;

simulationTime = 180;      % 3 hours
stepSize       = 1/60;     % 1 second steps in minutes
nSteps         = round(simulationTime / stepSize);
stepSize_in_seconds = stepSize * 60;

F = buildIntegrator(f, [Nx,Nu], stepSize_in_seconds);

%% Nonlinear Simulation
Ut = repmat(U',nSteps,1);  % nSteps x 6
xout = zeros(Nx,nSteps);
yout = zeros(Ny,nSteps);
x_current = Xz0;

for k=1:nSteps
    x_current = full(F(x_current, Ut(k,:)'));
    xout(:,k) = x_current;
    yout(:,k) = g(x_current);
    if mod(k,100)==0
        fprintf('Nonlinear simulation progress: %d out of %d steps\n',k,nSteps);
    end
end

timeVec = (1:nSteps)*stepSize; % in minutes
z = linspace(0,L,nstages);
stateVars = {'xB','xH2','xP','xT','T','Activity'};

%% Piecewise-Linear Approximation with Re-Linearization Every Step
chunkSize = 1; % crazy diverging numbers when using 100, relinerize every time step
nChunks = ceil(nSteps/chunkSize);

x_lin = zeros(Nx,nSteps);
x_lin(:,1)=xout(:,1); % start from the same initial condition
dx = zeros(Nx,1);

for c=1:nChunks
    startIdx = (c-1)*chunkSize + 1;
    endIdx = min(c*chunkSize, nSteps);
    
    % Linearize around the nonlinear state at start of this chunk
    x_bar = xout(:,startIdx);
    u_bar = Ut(startIdx,:)';    
    [A_lin, B_lin, C_lin, D_lin] = linearizeSystem(x_bar, u_bar, theta);

    fprintf('Linearization at step %d: chunk %d of %d\n', startIdx, c, nChunks);
    
    % Reset dx for this chunk (we define x - x_bar)
    dx = x_lin(:,startIdx)-x_bar;
    
    for k = startIdx:endIdx-1
        dx = dx + A_lin*dx*stepSize_in_seconds; % du=0 since input constant
        x_lin(:,k+1) = x_bar + dx;
    end
end

%% Plot Comparison: Nonlinear vs Piecewise-Linear
timeInstances = [1, round(nSteps/6), round(nSteps/3), round(nSteps/2), nSteps];
timeLabels = arrayfun(@(t) sprintf('%.1f min', (t-1)*stepSize), timeInstances,'UniformOutput',false);
colors = lines(length(timeInstances));

figure('Name','State Variables Along Reactor Length','NumberTitle','off','Color','w');
for iState = 1:nstates
    subplot(3,2,iState);
    hold on;
    % Nonlinear: solid
    for idx = 1:length(timeInstances)
        tIdx = timeInstances(idx);
        plot(z, xout((iState-1)*nstages+(1:nstages), tIdx), 'Color', colors(idx,:), 'LineWidth', 2);
    end
    % Linear: dashed
    for idx = 1:length(timeInstances)
        tIdx = timeInstances(idx);
        plot(z, x_lin((iState-1)*nstages+(1:nstages), tIdx), '--', 'Color', colors(idx,:), 'LineWidth', 1.5);
    end

    xlabel('Reactor Length [m]','FontSize',12);
    ylabel(stateVars{iState},'FontSize',12);
    title(['Spatial Profile of ', stateVars{iState}], 'FontSize',14);
    grid on; set(gca,'FontSize',12);

    if iState == 1
        legend([timeLabels, strcat(timeLabels,' (Linear)')],'Location','Best','FontSize',10);
    end
end
sgtitle('State Variables Along Reactor Length: Nonlinear (solid) vs Piecewise-Linear (dashed)','FontSize',16);

% Temporal profiles at selected positions
spaceStep4plot = round(0.1 * nstages);
spaceIndices = 1:spaceStep4plot:nstages;
spaceLabels = arrayfun(@(x) sprintf('%.2f m', z(x)), spaceIndices,'UniformOutput',false);

figure('Name','State Variables Over Time','NumberTitle','off','Color','w');
for iState = 1:nstates
    subplot(3,2,iState);
    hold on;
    colorIdx = 1;
    for sIdx = spaceIndices
        % Nonlinear
        plot(timeVec, xout((iState-1)*nstages + sIdx, :), 'Color', colors(colorIdx, :), 'LineWidth', 2);
        % Piecewise-Linear
        plot(timeVec, x_lin((iState-1)*nstages + sIdx, :), '--', 'Color', colors(colorIdx, :), 'LineWidth',1.5);
        colorIdx = colorIdx + 1;
    end
    xlabel('Time [min]','FontSize',12);
    ylabel(stateVars{iState},'FontSize',12);
    title(['Temporal Profile of ', stateVars{iState}], 'FontSize',14);
    grid on; set(gca,'FontSize',12);
    if iState == 1
        legend([spaceLabels, strcat(spaceLabels,' (Linear)')], 'Location','best','FontSize',10);
    end
end
sgtitle('State Variables Over Time at Positions: Nonlinear (solid) vs Piecewise-Linear (dashed)','FontSize',16);

% Contour plot for Temperature comparison
[Z_grid, T_grid] = meshgrid(z, timeVec);
TemperatureData_NL = xout(4*nstages+(1:nstages),:)';  
TemperatureData_LIN = x_lin(4*nstages+(1:nstages),:)';

figure('Name','Temperature Contour Comparison','NumberTitle','off','Color','w');
subplot(1,2,1);
contourf(Z_grid, T_grid, TemperatureData_NL, 50, 'LineStyle', 'none');
colorbar;
xlabel('Reactor Length [m]','FontSize',12);
ylabel('Time [min]','FontSize',12);
title('Nonlinear Temperature','FontSize',14);
set(gca,'FontSize',12);

subplot(1,2,2);
contourf(Z_grid, T_grid, TemperatureData_LIN, 50, 'LineStyle', 'none');
colorbar;
xlabel('Reactor Length [m]','FontSize',12);
ylabel('Time [min]','FontSize',12);
title('Piecewise-Linear Temperature','FontSize',14);
set(gca,'FontSize',12);

sgtitle('Temperature Distribution: Nonlinear vs Piecewise-Linear (relinearized every 100 steps)','FontSize',16);

disp('Piecewise-linear approximation complete. A,B,C,D matrices available at each linearization point.');
disp('Simulation and plotting finished.');
