clear all; clc; close all

%% ADD CASADI PATH IF NEEDED
addpath('C:\Users\Jesse_N\OneDrive - Aalto University\Documents\Masters Thesis\starting code\casadi-3.6.6')
import casadi.*; 
%% Model parameters

k0        =  3.16 * 10^-2;     % (kmol/(kg cat * s * Pa))   Pre exponential factor for hydrogenation   !
K0        =  3.16 * 10^-13;    % (1/Pa)                     Adsorption constant for benzene            !
Q         = -6.89 * 10^7;      % (kJ/kmol)                  Benzene adsorption activation energy       !
E         =  5.76 * 10^7;      % (kJ/kmol)                  Activation energy of hydrogenation         !

rhoCp     =  7.322 * 10^5;     % (J/(m^3*C))                Effective heat capacity, of reactor        !
                               %                             (First value in the paper, Table 1 and 2)
negDeltaH =  2.09 * 10^8;      % (J/kmol)                   Heat of reaction                           !

rhoB      =  4.14 * 10^2;      % (kg/m^3)                   Catalyst bulk density                      ! 
                               %                             (First value in the paper, also Table 1)  
epz       =  0.58;             %                            Void fraction                              !

MT        =  4.05 * 10^-4;     % (kg/kmol)                  Poison adsorption capacity of catalyst     !
                               %                             (First value in the paper, Table 1)

k0d       =  1.80 * 10^-4;     % ((Pa*s)^-1)                Pre exponential factor for poisoning       !
Ed        =  4.52 * 10^6;      % (J/kmol)                   Activation energy of poisoning             !

P         =  1.01*10^5;        % (Pa or N/m^2)              Pressure (Tab.1: Example)                  ! 

L         =  0.118 * 5;        % (m)                        Catalysist's bed length                    ! <= Must be fixed to account
                               %                             (First value in Table 1)                  !       for the actual length
r         =  0.0078;           % (m)                        Reactor's radius                           !
DeB       =  4.5 * 10^-5;      % 

nstages   = 100; 

front = 0.2;
back  = 0.8;
anteCatalyst = 1:floor(front*nstages);
postCatalyst = floor(back*nstages)+1:nstages;
catalystBed  = anteCatalyst(end)+1 : postCatalyst(1)-1;

sensorNumber   = 5;
sensorPosition = round(linspace(catalystBed(1),catalystBed(end),sensorNumber));

theta     = {nstages, k0, K0, Q, E, rhoCp, negDeltaH, rhoB, epz, MT, k0d, Ed, P, L, r, DeB, sensorPosition};
%            1        2   3   4  5  6      7          8     9    10  11   12  13 14 15 16   17

A         =  pi*r^2;          % (m^2)                      Reactor's cross-section area               !

%% Model controls and disturbances 

% xB = 1.42 * 10^-2;                                % (kmol/kmol)       Benzene mole fraction at inlet (disturbance)
% xT = 6.38 * 10^-4;                                % (kmol/kmol)       Tiophenol mole fraction at inlet (disturbance)
% T0 = 50;                                          % (C)               Inlet temperature (control)

%% Model equations
f = @(x,u) modelCD(0, x, u, theta);                 
g = @(x)   modelCD_out(0, x, theta);                
                                                    % temperatures as measurements
                                                    
%% Integrator

% Variables
nstates  = 6;

Nx = nstates * nstages;                             % [xb, xh, xp, xt, T, theta]
Ny = 1;                                             % To be defined

Nu = 6;                                             % 1 (or 2) control [T0 (and U)]
Nw = 5;                                             % 5 disturbances [XB0, XH, XT0, XP0, T0, U] 

% Times (model is in seconds) 

simulationTime = 3 * 60;                            % Minutes (set to 3h)
stepSize       = 1 / 60;                            % Minutes (set to 1s)
nSteps         = round(simulationTime / stepSize);  % 

simulationTime_in_seconds = simulationTime * 60;    %
stepSize_in_seconds       = stepSize * 60;          %

F = buildIntegrator(f, [Nx,Nu], stepSize_in_seconds); %F = buildIntegrator(f, [Nx,Nu]);

%% Simulate    

% Boundary and initial conditions
xBzt = 1.5E-2 * ones(nstages,nSteps); 
xPzt = 0.0    * ones(nstages,nSteps);
xTzt = 0.0    * ones(nstages,nSteps);
xHzt = 1.0 - xBzt - xPzt;

Tzt = 350 * ones(nstages,nSteps);

Azt =  zeros(nstages,nSteps); 
Azt(catalystBed,1:nSteps) = 1;

Xz0 = [xBzt(:,1);                                  % NSTAGE-long vector of XBs
       xHzt(:,1);                                  % NSTAGE-long vector of XHs
       xPzt(:,1);                                  % NSTAGE-long vector of XPs
       xTzt(:,1);                                  % NSTAGE-long vector of xTs
       Tzt(:,1);                                   % NSTAGE-long vector of tTs
       Azt(:,1)                                    % NSTAGE-long vector of xAs
      ];
  
% Inputs 
Feed = 2.5 * 10^-5;                                % (m^3 / sec) Volumetric flowrate
ut   = (Feed/(A*epz))*ones(1,nSteps);              % (m / sec)   Linear velocity    TF-long vector of Us   (disturbance/control)

xB0t = 1 * xBzt(1,1) * ones(1,nSteps);             % ()          Feed xB            TF-long vector of xB0  (disturbance)
xP0t = 1 * xPzt(1,1) * ones(1,nSteps);             % ()          Feed xP            TF-long vector of xP0  (disturbance)
xT0t = 6.36E-4 + xTzt(1,1) * ones(1,nSteps);       % ()          Feed xT            TF-long vector of xT0  (disturbance)
xH0t = 1 - xB0t - xT0t;                            % ()          Feed xH            TF long vector of XH0  (disturbance)

T0t  = 1 * Tzt(1,1)  * ones(1,nSteps);             % (C)         Feed T             TF-long vector of T0s  (control)

Ut   = [xB0t', xH0t', xP0t', xT0t', T0t', ut'];    %             Put stuff together

% Simulation  
[yout,~,xout] = simulateSystem(F, g, Xz0, Ut);
 
%% Plotting

stateVars = {'xB','xH2','xP','xT','T','Activity'};

% Time step plotting interval
timeStep4plot = round(0.1 * nSteps); % Every 10% of the total steps
% Correct legend based on total simulation time of 180 minutes
timeLabels = arrayfun(@(x) sprintf('%.1f min', (x-1) * 180 / (nSteps/timeStep4plot)), 1:(nSteps/timeStep4plot), 'UniformOutput', false);

% Use a colormap with distinct colors for each line
colors = lines(length(1:timeStep4plot:nSteps));

%% Figure 1:
figure(1); 
for iState = 1:nstates
    subplot(3, 2, iState); 
    hold on;
    % Adjust the color index to match the number of time steps
    colorIdx = 1;
    for tIdx = 1:timeStep4plot:nSteps
        plot(linspace(0, 1, nstages), xout((iState-1)*nstages + 1 : iState*nstages, tIdx), 'Color', colors(colorIdx, :), 'LineWidth', 1.5);
        colorIdx = colorIdx + 1; % Increment color index
    end
    if iState == 6
        xlabel('Reactor length [%]');
    end
    axis tight; 
    title(stateVars{iState}); 
    set(gca, 'FontSize', 18);
end

% Single legend for the figure (Correct Time Legend in Minutes)
legend(timeLabels, 'Location', 'BestOutside');

%% Figure 2: Along Time with Correct Reactor Position Legend (in reactor %)

% Space step plotting interval
spaceStep4plot = round(0.1 * nstages); % Every 10% of the reactor length
spaceLabels = arrayfun(@(x) sprintf('Reactor Position: %.1f %%', (x-1) * 100 / nstages), 1:spaceStep4plot:nstages, 'UniformOutput', false);

% Use a colormap with distinct colors for each line
colors = lines(length(1:spaceStep4plot:nstages));

figure(2); 
for iState = 1:nstates
    subplot(3, 2, iState); 
    hold on;
    % Adjust the color index to match the number of space steps
    colorIdx = 1;
    for sIdx = 1:spaceStep4plot:nstages
        % Correct index for spatial points, no tIdx involved here
        stairs((1:nSteps) * stepSize, xout((iState-1) * nstages + sIdx, 1:nSteps), 'Color', colors(colorIdx, :), 'LineWidth', 1.5);
        colorIdx = colorIdx + 1; % Increment color index
    end
    if iState == 6
        xlabel('Time [min]');
    end
    axis tight; 
    title(stateVars{iState});
    set(gca, 'FontSize', 18);
end

% Single legend for the figure (Reactor Position Legend)
legend(spaceLabels, 'Location', 'BestOutside');

%% Figure 1

% Define reactor length vector
z = linspace(0, L, nstages); % Reactor length from 0 to L

% Time steps to plot (e.g., every 30 minutes)
timeInstances = [1, round(nSteps/6), round(nSteps/3), round(nSteps/2), nSteps];
timeLabels = arrayfun(@(t) sprintf('%.1f min', (t-1) * stepSize), timeInstances, 'UniformOutput', false);

% Colors for the lines
colors = lines(length(timeInstances));

% Create Figure 1
figure('Name', 'State Variables Along Reactor Length', 'NumberTitle', 'off', 'Color', 'w');
for iState = 1:nstates
    subplot(3, 2, iState);
    hold on;
    for idx = 1:length(timeInstances)
        tIdx = timeInstances(idx);
        plot(z, xout((iState-1)*nstages + (1:nstages), tIdx), 'Color', colors(idx, :), 'LineWidth', 1.5);
    end
    xlabel('Reactor Length [m]', 'FontSize', 12);
    ylabel(stateVars{iState}, 'FontSize', 12);
    title(['Spatial Profile of ', stateVars{iState}], 'FontSize', 14);
    grid on;
    set(gca, 'FontSize', 12);
    if iState == 1
        legend(timeLabels, 'Location', 'Best', 'FontSize', 10);
    end
end
sgtitle('State Variables Along Reactor Length at Different Time Steps', 'FontSize', 16);

%% FIGURE 2

% Reactor positions to plot
positionIndices = [1, round(nstages/4), round(nstages/2), round(3*nstages/4), nstages];
positionLabels = arrayfun(@(pos) sprintf('%.2f m', z(pos)), positionIndices, 'UniformOutput', false);

% Colors for the lines
colors = lines(length(positionIndices));

% Create Figure 2
figure('Name', 'State Variables Over Time', 'NumberTitle', 'off', 'Color', 'w');
for iState = 1:nstates
    subplot(3, 2, iState);
    hold on;
    for idx = 1:length(positionIndices)
        sIdx = positionIndices(idx);
        plot((1:nSteps) * stepSize, xout((iState-1)*nstages + sIdx, :), 'Color', colors(idx, :), 'LineWidth', 1.5);
    end
    xlabel('Time [min]', 'FontSize', 12);
    ylabel(stateVars{iState}, 'FontSize', 12);
    title(['Temporal Profile of ', stateVars{iState}], 'FontSize', 14);
    grid on;
    set(gca, 'FontSize', 12);
    if iState == 1
        legend(positionLabels, 'Location', 'Best', 'FontSize', 10);
    end
end
sgtitle('State Variables Over Time at Different Reactor Positions', 'FontSize', 16);

%% FIGURE 3

% Time instances to plot
timeInstances = [1, round(nSteps/3), round(2*nSteps/3), nSteps];
timeLabels = arrayfun(@(t) sprintf('%.1f min', (t-1) * stepSize), timeInstances, 'UniformOutput', false);

% Colors for the lines
colors = lines(length(timeInstances));

% Create Figure 3
figure('Name', 'Catalyst Activity Profile', 'NumberTitle', 'off', 'Color', 'w');
hold on;
for idx = 1:length(timeInstances)
    tIdx = timeInstances(idx);
    plot(z, xout(5*nstages + (1:nstages), tIdx), 'LineWidth', 2, 'Color', colors(idx, :));
end
xlabel('Reactor Length [m]', 'FontSize', 12);
ylabel('Catalyst Activity', 'FontSize', 12);
title('Catalyst Activity Profile Along Reactor Length', 'FontSize', 16);
legend(timeLabels, 'Location', 'Best', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 12);


%% FIGURE 4

% Prepare data for contour plot
[Z, T_grid] = meshgrid(z, (1:nSteps) * stepSize);
TemperatureData = xout(4*nstages + (1:nstages), :)';

% Create Figure 4
figure('Name', 'Temperature Contour Plot', 'NumberTitle', 'off', 'Color', 'w');
contourf(Z, T_grid, TemperatureData, 50, 'LineStyle', 'none');
colorbar;
xlabel('Reactor Length [m]', 'FontSize', 12);
ylabel('Time [min]', 'FontSize', 12);
title('Temperature Distribution Over Time and Reactor Length', 'FontSize', 16);
set(gca, 'FontSize', 12);

