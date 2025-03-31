%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% directFeedback_control.m
%
% This script runs the reactor simulation with direct feedback control,
% produces plots and a video of the reactor states.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all

%% ADD OWN CASADI PATH
import casadi.*

%% Model parameters (unchanged)
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
A         = pi*r^2;

nstages = 50; %can modify -- depends on discretization goal
front = 0.2;
back  = 0.8;
anteCatalyst = 1:floor(front*nstages);
postCatalyst = floor(back*nstages)+1:nstages;
catalystBed  = anteCatalyst(end)+1 : postCatalyst(1)-1;
sensorNumber = 5;
sensorPosition = round(linspace(catalystBed(1),catalystBed(end),sensorNumber));

theta = {nstages, k0, K0, Q, E, rhoCp, negDeltaH, rhoB, epz, MT, k0d, Ed, P, L, r, DeB, sensorPosition};

% Initial conditions
xBzt = 1.5e-2*ones(nstages,1);
xPzt = zeros(nstages,1);
xTzt = zeros(nstages,1);
xHzt = 1.0 - xBzt - xPzt - xTzt; 
Tzt  = 350*ones(nstages,1);
Azt  = zeros(nstages,1);
Azt(catalystBed)=1;
X0 = [xBzt; xHzt; xPzt; xTzt; Tzt; Azt];

Feed = 2.5e-5;
ut   = (Feed/(A*0.58));
xB_in = xBzt(1);  % inlet benzene fraction reference

xB0 = xBzt(1); xP0=xPzt(1); xT0=6.36e-4+xTzt(1); xH0=1 - xB0 - xT0 - xP0;
T0_in=350; 
U0=[xB0,xH0,xP0,xT0,T0_in,ut];

f = @(x,u) modelCD(0,x,u,theta);
g = @(x) modelCD_out(0,x,theta);

nstates = 6; Nx=nstates*nstages; Nu=6;

% Simulation settings
simulationTime=180;   % run for 3 hours = 180 min
stepSize=1/600;       % 1-second step in minutes
nSteps=round(simulationTime/stepSize);
stepSize_s=stepSize*60;

F=buildIntegrator(f,[Nx,Nu],stepSize_s);

%% Control Settings
conv_target = 0.9;  
Kp = 5;  % proportional gain

% Arrays for storing data
X_nl=zeros(Nx,nSteps);X_nl(:,1)=X0;
Y_meas=zeros(length(g(X0)),nSteps);Y_meas(:,1)=g(X0);
U_hist=zeros(nSteps,Nu);U_hist(1,:)=U0;

T_in_hist=zeros(1,nSteps);T_in_hist(1)=U0(5);
cumulative_product=zeros(1,nSteps);

fprintf('Starting direct feedback control simulation...\n');
tic;

for k=1:nSteps-1
    x_next=full(F(X_nl(:,k),U_hist(k,:)'));
    if any(isnan(x_next))||any(isinf(x_next))
        x_next = max(x_next,1e-12);
    end
    X_nl(:,k+1)=x_next;
    Y_meas(:,k+1)=g(X_nl(:,k+1));

    % Compute conversion at outlet
    xB_block = X_nl(1:nstages,k+1);
    xB_out = xB_block(end);
    conversion = (xB_in - xB_out)/xB_in;

    % Compute error and adjust T_in
    err = conversion - conv_target;
    delta_T_in = -Kp * err;
    T_in_new = U_hist(k,5) + delta_T_in;
    
    % Clamp T_in_new to avoid unrealistic values
    T_in_new = min(max(T_in_new,300),500);

    % Update input
    U_hist(k+1,:)=U_hist(k,:);
    U_hist(k+1,5)=T_in_new;
    
    T_in_hist(k+1)=T_in_new;
    cumulative_product(k+1)=cumulative_product(k)+X_nl(3*nstages,k)*Feed*stepSize_s;

    if mod(k,round(nSteps/10))==0
        fprintf('Progress: %d/%d steps done.\n',k,nSteps);
    end
end

total_time=toc;
fprintf('Simulation completed in %.3f s.\n',total_time);
fprintf('Final Cumulative Product: %.6f\n', cumulative_product(end));

timeVec=(0:nSteps-1)*stepSize; 
z=linspace(0,L,nstages);

% Compute conversion over time
conversion_time = zeros(1,nSteps);
for i=1:nSteps
    xB_out_i = X_nl(nstages,i);
    conversion_time(i)=(xB_in - xB_out_i)/xB_in;
end

% Extract T-block for visualization (temperature)
T_block  = X_nl(4*nstages+1:5*nstages,:);

%% Plotting

% Plot Inlet Temperature over Time
figure('Name','Inlet Temperature');
plot(timeVec,T_in_hist,'b-','LineWidth',2);
xlabel('Time [min]','FontSize',12);
ylabel('Inlet Temperature [°C]','FontSize',12);
title('Inlet Temperature over Time','FontSize',14);
grid on;

% Plot Conversion vs Target
figure('Name','Conversion vs Target');
plot(timeVec,conversion_time,'b','LineWidth',2); hold on;
plot(timeVec,conv_target*ones(size(timeVec)),'r--','LineWidth',2);
xlabel('Time [min]','FontSize',12);
ylabel('Conversion (dimensionless)','FontSize',12);
legend('Conversion','Target','Location','best');
title('Outlet Conversion vs Target','FontSize',14);
grid on;

% Cumulative Product
figure('Name','Cumulative Product');
plot(timeVec,cumulative_product,'b-','LineWidth',2);
xlabel('Time [min]','FontSize',12);ylabel('Cumulative Product','FontSize',12);
title('Cumulative Product Over Time','FontSize',14);
grid on;

% Temperature Contour
[Z_grid,T_grid]=meshgrid(z,timeVec);
TemperatureData=T_block';
figure('Name','Temperature Contour');
contourf(Z_grid,T_grid,TemperatureData,50,'LineStyle','none');
colorbar;
xlabel('Reactor Length [m]','FontSize',12);
ylabel('Time [min]','FontSize',12);
title('Temperature Distribution Over Time and Length','FontSize',14);
grid on;

% Select snapshot times - ensure they are within simulationTime
snapshot_times = [0,60,120,180]; 

% Find snapshot indices with a for loop to handle empties
snapshot_indices = zeros(size(snapshot_times));
for ii = 1:length(snapshot_times)
    idx = find(timeVec >= snapshot_times(ii), 1, 'first');
    if isempty(idx)
        idx = length(timeVec); % fallback to the last time step
    end
    snapshot_indices(ii) = idx;
end

xH_block = X_nl(nstages+1:2*nstages,:);
xP_block = X_nl(2*nstages+1:3*nstages,:);
xT_block = X_nl(3*nstages+1:4*nstages,:);
A_block  = X_nl(5*nstages+1:6*nstages,:);

figure('Name','State Spatial Profiles','Color','w');

subplot(3,2,1); hold on;
for i=1:length(snapshot_indices)
    idx=snapshot_indices(i);
    plot(z,X_nl(1:nstages,idx),'LineWidth',2,'DisplayName',sprintf('%d min',snapshot_times(i)));
end
xlabel('Length [m]','FontSize',12);ylabel('x_B','FontSize',12);
title('Benzene (xB)','FontSize',14);
legend('Location','best');grid on;

subplot(3,2,2); hold on;
for i=1:length(snapshot_indices)
    idx=snapshot_indices(i);
    plot(z,xH_block(:,idx),'LineWidth',2,'DisplayName',sprintf('%d min',snapshot_times(i)));
end
xlabel('Length [m]','FontSize',12);ylabel('x_H','FontSize',12);
title('Hydrogen (xH)','FontSize',14);
legend('Location','best');grid on;

subplot(3,2,3); hold on;
for i=1:length(snapshot_indices)
    idx=snapshot_indices(i);
    plot(z,xP_block(:,idx),'LineWidth',2,'DisplayName',sprintf('%d min',snapshot_times(i)));
end
xlabel('Length [m]','FontSize',12);ylabel('x_P','FontSize',12);
title('Product (xP)','FontSize',14);
legend('Location','best');grid on;

subplot(3,2,4); hold on;
for i=1:length(snapshot_indices)
    idx=snapshot_indices(i);
    plot(z,xT_block(:,idx),'LineWidth',2,'DisplayName',sprintf('%d min',snapshot_times(i)));
end
xlabel('Length [m]','FontSize',12);ylabel('x_T','FontSize',12);
title('Poison (xT)','FontSize',14);
legend('Location','best');grid on;

subplot(3,2,5); hold on;
for i=1:length(snapshot_indices)
    idx=snapshot_indices(i);
    plot(z,T_block(:,idx),'LineWidth',2,'DisplayName',sprintf('%d min',snapshot_times(i)));
end
xlabel('Length [m]','FontSize',12);ylabel('T [°C]','FontSize',12);
title('Temperature (T)','FontSize',14);
legend('Location','best');grid on;

subplot(3,2,6); hold on;
A_block(A_block<1e-12)=1e-12;
for i=1:length(snapshot_indices)
    idx=snapshot_indices(i);
    plot(z,A_block(:,idx),'LineWidth',2,'DisplayName',sprintf('%d min',snapshot_times(i)));
end
xlabel('Length [m]','FontSize',12);ylabel('Activity','FontSize',12);
title('Catalyst Activity (A)','FontSize',14);
legend('Location','best');grid on;

sgtitle('State Spatial Profiles at Hourly Intervals','FontSize',16);

%% After the simulation and plotting of static figures:

% Assume the simulation has run, and we have:
% timeVec, conversion_time, z, T_block, and variables defined.

% Speed up the video: For example, if stepSize is 1/600 min, each step = 1 sec.
% To show 1 minute per frame, capture every 600 steps.

stride = 600; % adjust as needed for desired speed-up

videoFilename = 'TemperatureEvolution_Fast.mp4';
v = VideoWriter(videoFilename,'MPEG-4');
v.Quality = 100;
v.FrameRate = 20;
open(v);

% Create a fixed-size figure for consistent frame dimensions
hfig = figure('Color','w','Name','Temperature Evolution Video',...
    'Units','pixels','Position',[100 100 560 420]); % fixed size 560x420 px

% Precompute max temperature for annotation
max_T = max(T_block(:));

for k = 1:stride:length(timeVec)
    plot(z, T_block(:,k), 'LineWidth', 2);
    title('Temperature Profile Over Time','FontSize',14);
    xlabel('Reactor Length [m]','FontSize',12); 
    ylabel('T [°C]','FontSize',12);
    grid on;

    % Add text annotation with current time and conversion
    txt = sprintf('Time: %.2f min\nConversion: %.3f', timeVec(k), conversion_time(k));
    text(z(end)*0.1, max_T*0.9, txt, 'FontSize',12, 'BackgroundColor','w','EdgeColor','k');

    drawnow;

    % Use getframe with a specified rectangle matching the figure size
    % [0 0 width height] relative to the figure's drawable area
    frame = getframe(hfig,[0 0 560 420]);

    writeVideo(v, frame);
end

close(v);
disp(['Video saved as ', videoFilename]);

