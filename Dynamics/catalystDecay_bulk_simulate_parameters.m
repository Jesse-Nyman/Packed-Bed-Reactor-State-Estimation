clear all; close all; clc

addpath('casadi-3.6.6\')
import casadi.*; 

%% Model parameters

k0        =  3.16 * 10^-2;     % (kmol/(kg cat * s * Pa))   Pre exponential factor for hydrogenation   !
K0        =  3.16 * 10^-13;    % (1/Pa)                     Adsorption constant for benzene            !
Q         = -6.89 * 10^7;      % (kJ/kmol)                  Benzene adsorption activation energy       !
E         =  5.76 * 10^7;      % (kJ/kmol)                  Activation energy of hydrogenation         !
rhoCp     =  7.322 * 10^5;     % (J/(m^3*C))                Effective heat capacity, of reactor        !
                               %                            (First value in the paper, Table 1 and 2)
negDeltaH =  2.09 * 10^8;      % (J/kmol)                   Heat of reaction                           !
rhoB      =  4.14 * 10^2;      % (kg/m^3)                   Catalyst bulk density                      ! 
                               %                             (First value in the paper, also Table 1)  
epz       =  0.58;             %                            Void fraction                              !
MT        =  4.05 * 10^-4;     % (kg/kmol)                  Poison adsorption capacity of catalyst     !
                               %                             (First value in the paper, Table 1)
k0d       =  1.80 * 10^-4;     % ((Pa*s)^-1)                Pre exponential factor for poisoning       !
Ed        =  4.52 * 10^6;      % (J/kmol)                   Activation energy of poisoning             !
P         =  1.01*10^5;        % (Pa or N/m^2)              Pressure (Tab.1: Example)                  ! 
L         =  0.118 *10/9;      % (m)                        Catalysist's bed length                    ! <= Must be fixed to account
                               %                             (First value in Table 1)                  !       for the actual length
r         =  0.0078;           % (m)                        Reactor's radius                           !
DeB       =  4.5 * 10^-5;      % 
A         =  pi*r^2;           % (m^2)                      Reactor's cross-section area               !

nstages   = 100; 

front = 0.05;
back  = 0.95;
anteCatalyst = 1:round(front*nstages);
postCatalyst = round(back*nstages)+1:nstages;
catalystBed  = anteCatalyst(end)+1 : postCatalyst(1)-1;

sensorNumber   = 5;
sensorPosition = round(linspace(catalystBed(1),catalystBed(end),sensorNumber)); 
%% Configure integrator

% Variables
nstates  = 6;

Nx = nstates * nstages;                             % [xb, xh, xp, xt, T, theta]
Ny = 1;                                             % To be defined

Nu = 6;                                             % 1 (or 2) control [T0 (and U)]
Nw = 5;                                             % 5 disturbances [XB0, XH, XT0, XP0, T0, U] 

% Times (model is in seconds) | Do not fully trust these time units

simulationTime = 2 * 60;                            % Minutes (set to 3h)
stepSize       = 1 / 60;                            % Minutes (set to 1s)
nSteps         = round(simulationTime / stepSize);  % 

simulationTime_in_seconds = simulationTime * 60;    %
stepSize_in_seconds       = stepSize * 60;          %
%% Inputs
Feed = 2.5 * 10^-5;                                 % (m^3 / sec) Volumetric flowrate

% Boundary and initial conditions
xBzt = 1.5E-2 * ones(nstages,nSteps);
xPzt = 0.0    * ones(nstages,nSteps);
xTzt = 0.0    * ones(nstages,nSteps);
xHzt = 1.0 - xBzt - xPzt;

Tzt = 350 * ones(nstages,nSteps);

Azt =  zeros(nstages,nSteps); 
Azt(catalystBed,1:nSteps) = 1;

Xz0 = [
       xBzt(:,1);                                  % NSTAGE-long vector of XBs
       xHzt(:,1);                                  % NSTAGE-long vector of XHs
       xPzt(:,1);                                  % NSTAGE-long vector of XPs
       xTzt(:,1);                                  % NSTAGE-long vector of xTs
       Tzt(:,1);                                   % NSTAGE-long vector of tTs
       Azt(:,1)                                    % NSTAGE-long vector of xAs
      ]; 

xB0t = 1 * xBzt(1,1) * ones(1,nSteps);             % ()          Feed xB            NSTEP-long vector of xB0  (disturbance)
xP0t = 1 * xPzt(1,1) * ones(1,nSteps);             % ()          Feed xP            NSTEP-long vector of xP0  (disturbance)
xT0t = 6.36E-4 + xTzt(1,1) * ones(1,nSteps);       % ()          Feed xT            NSTEP-long vector of xT0  (disturbance)
xH0t = 1 - xB0t - xT0t;                            % ()          Feed xH            NSTEP long vector of XH0  (disturbance)

T0t  = 1 * Tzt(1,1)  * ones(1,nSteps);             % (C)         Feed T             NSTEP-long vector of T0s  (control)

for j = 1:10

 %% Disturbances 
    prcnt1 = 5*10^-2;                                                       %B
    prcnt2 = 5*10^-2;                                                       %T
    prcnt3 = 5*10^-2;                                                       %P
    prcnt4 = 5*10^-2;                                                       % temperature
    prcnt5 = 5*10^-2;                                                       % Feed
    p = 0.05;                                                               %Probability of success
    T0t = temperature_disturbance(nSteps,T0t,Tzt(1,1),prcnt4,p);
    ut   = Feed_Disturbances(nSteps,Feed,A,epz,prcnt5,p);              % (m / sec)   Linear velocity    NSTEP-long vector of Us   (disturbance/control)
    

    Feed_C = [xB0t', xH0t', xP0t', xT0t'];
    Feed_C  = Concentration_Disturbances(nSteps,xB0t(1,1), xT0t(1,1), xP0t(1,1),prcnt1,prcnt2,prcnt3,p, Feed_C);
    Ut = [Feed_C(:,1),Feed_C(:,2),Feed_C(:,3),Feed_C(:,4) T0t', ut'];
   
    %% Simulation loop
    for n = 1:10
        
        %% Model controls and disturbances 

        % xB = 1.42 * 10^-2;                                % (kmol/kmol)       Benzene mole fraction at inlet (disturbance)
        % xT = 6.38 * 10^-4;                                % (kmol/kmol)       Tiophenol mole fraction at inlet (disturbance)
        % T0 = 50;                                          % (C)               Inlet temperature (control)

        %% Sample parameters, uniform +/- some percent (<= 1) around value in paper

        prcnt       = 0.10 * 10^-2;

         k0_          = unifrnd((1-prcnt)*k0        , (1+prcnt)*k0       , 1);  
         K0_          = unifrnd((1-prcnt)*K0        , (1+prcnt)*K0       , 1);  
         Q_           = unifrnd((1+prcnt)*Q         , (1-prcnt)*Q        );  
         E_           = unifrnd((1-prcnt)*E         , (1+prcnt)*E        , 1);  
         rhoCp_       = unifrnd((1-prcnt)*rhoCp     , (1+prcnt)*rhoCp    , 1); 
         negDeltaH_   = unifrnd((1-prcnt)*negDeltaH , (1+prcnt)*negDeltaH, 1); 
         rhoB_        = unifrnd((1-prcnt)*rhoB      , (1+prcnt)*rhoB     , 1); 
        % epz_         = unifrnd((1-prcnt)*epz       , (1+prcnt)*epz      , 1);
         MT_          = unifrnd((1-prcnt)*MT        , (1+prcnt)*MT       , 1);
         k0d_         = unifrnd((1-prcnt)*k0d       , (1+prcnt)*k0d      , 1);
         Ed_          = unifrnd((1-prcnt)*Ed        , (1+prcnt)*Ed       , 1);
         P_           = unifrnd((1-prcnt)*P         , (1+prcnt)*P        , 1);
        % L_           = unifrnd((1-prcnt)*L         , (1+prcnt)*L        , 1);  
        % r_           = unifrnd((1-prcnt)*r         , (1+prcnt)*r        , 1);
         DeB_         = unifrnd((1-prcnt)*DeB       , (1+prcnt)*DeB      , 1);

        theta    = {nstages, k0_, K0_, Q_, E_, rhoCp_, negDeltaH_, rhoB_, epz, MT_, k0d_, Ed_, P_, L, r, DeB_, sensorPosition}; % needs to be made compatible with parfor
        %            1        2    3    4   5   6        7            8     9    10  11    12    13 14 15 16    17

        %% Define model equations and integrator
        f = @(x,u) modelCD( 0, x, u, theta);                % Must fix this to include disturbances as inputs (?)
        g = @(x)   modelCD_out(0, x, theta);                % Must finish this to extract some side 
                                                            % temperatures as measurements

        F = buildIntegrator(f, [Nx,Nu], stepSize_in_seconds); %F = buildIntegrator(f, [Nx,Nu]);


      
        %% Simulation  
        [yout,~,xout] = simulateSystem(F, g, Xz0, Ut);

        %% Storing data in txt file
        %
% 
%         dlmwrite(sprintf('%s%d_%d_%s','Data\',j,n,'xout.txt'),xout);
%         dlmwrite(sprintf('%s%d_%d_%s','Data\',j,n,'yout.txt'),yout);
%         dlmwrite(sprintf('%s%d_%d_%s','Data\',j,n,'Ut.txt'),Ut);
%         dlmwrite(sprintf('%s%d_%d_%s','Data\',j,n,'theta.txt'),theta);

        %% Store data in mat file
        Parsave(j,n,xout',yout',Ut,theta);
%         save(sprintf('%s%d%d_%s','Data\',j,n,'xout.mat'),'xout');
%         save(sprintf('%s%d%d_%s','Data\',j,n,'yout.mat'),'yout');
%         save(sprintf('%s%d%d_%s','Data\',j,n,'Ut.mat'),'Ut');
%         save(sprintf('%s%d%d_%s','Data\',j,n,'theta.mat'),'theta');
    end
end

%% Plotting

stateVars = {'xB','xH2','xP','xT','T','Activity'};

timeStep4plot = round(0.1*nSteps);
timeColor = gray(numel(1:timeStep4plot:nSteps));

spaceStep4plot = round(0.1*nstages);
spaceColor = gray(numel(1:spaceStep4plot:nstages));

figure(); % Along space
for iState = 1:nstates
 subplot(3,2,iState); 
    plot(linspace(0,1,nstages),xout((iState-1)*nstages + 1 : iState*nstages, 1:timeStep4plot:nSteps));
    if iState == 6; xlabel('Reactor length [%]'); end
    axis tight; title(stateVars(iState))
    set(gca,'FontSize',18)
end

figure() % Along time
for iState = 1:nstates
   subplot(3,2,iState); 
    stairs((1:nSteps)*stepSize,xout((iState-1)*nstages+1 : 10 : iState*nstages, 1:nSteps)'); 
    if iState == 6; xlabel('Time [min]'); end
    axis tight; title(stateVars(iState))  
    set(gca,'FontSize',18)
end

figure() % Input-Output
plot(1:nSteps,ut)