clear all; import casadi.*; clc; clf

%% Model parameters

k0        =  3.16 * 10^-2;     % (kmol/(kg cat * s * Pa))   Pre exponential factor for hydrogenation   !
K0        =  3.16 * 10^-13;    % (1/Pa)                     Adsorption constant for benzene            !
Q         = -6.89 * 10^7;      % (kJ/kmol)                  Benzene adsorption activation energy       !
E         =  5.76 * 10^7;      % (kJ/kmol)                  Activation energy of hydrogenation         !

rhoCp     =  7.322 * 10^5;     % (J/(m^3*C))                Effective heat capacity, of reactor        !
                               %                            (First value in the paper, Table 1 and 2)
negDeltaH =  2.09 * 10^8;      % (J/kmol)                   Heat of reaction                           !

rhoB      =  4.14 * 10^2;     % (kg/m^3)                   Catalyst bulk density                      ! <=
                               %                            (First value in the paper, also Table 1)  
epz       =  0.58;             %                            Void fraction                              !

MT        =  4.05 * 10^-4;     % (kg/kmol)                  Poison adsorption capacity of catalyst     !
                               %                            (First value in the paper, Table 1)

k0d       =  1.80 * 10^-4;     % ((Pa*s)^-1)                Pre exponential factor for poisoning       !
Ed        =  4.52 * 10^6;      % (J/kmol)                   Activation energy of poisoning             !

P         =  1.01*10^5;         % (Pa or N/m^2)              Pressure (Tab.1: Example)                  ! 

L         =  0.118 /(3/5);            % (m)                        Catalysist's bed length                    !
                               %                            (First value in Table 1) 
r         =  0.0078;            % (m)                        Reactor's radius                           !
DeB       =  4.5*10^-5;                           
nstages = 50;
M1 = 1; 
M2 = nstages*1/5+1;
M3 = nstages*2/5+1;
M4 = nstages*3/5+1;
M5 = nstages*4/5;
M6 = nstages;
theta = [nstages, k0, K0, Q, E, rhoCp, negDeltaH, rhoB, epz, MT, k0d, Ed, P, L, r, DeB, M1,M2,M3,M4,M5,M6];

A         =  pi*r^2;           % (m^2)                      Reactor's cross-section area               !
 
% mwB       =  78.11;          % (kg/kmol)                  Benzene molar mass (Google)                ! <=
% Brho      =  876;            % (kg/m^3)                   Benzene density    (Google)                ! <=

%% Model controls and disturbances 

% xB = 1.42 * 10^-2;             % (kmol/kmol)       Benzene mole fraction at inlet (disturbance)
% xT = 6.38 * 10^-4;             % (kmol/kmol)       Tiophenol mole fraction at inlet (disturbance)
% T0 = 50;                       % (C)               Inlet temperature (control)

%% Model equations
f = @(x,u) modelCD(0, x, u, theta);                 % Must fix this to include disturbances as inputs (?)
g = @(x)   modelCD_out(0, x, theta);                % Must finish this to extract some side 
                                                    % temperatures as measurements

%% Build integrator
nstates  = 6;

Nx = nstates*nstages;           % [xb, xh, xp, xt, T, theta]
Ny = 1;                         % To be defined

Nu = 6;                         % 1 (or 2) control [T0 (and U)]
Nw = 5;                         % 5 disturbances [XB0, XH, XT0, XP0, T0, U] 

F = buildIntegrator(f,[Nx,Nu]);

%% Simulate model    

% Final time (units, minutes?)
finalTime  =6000; 

% Boundary and initial conditions
xBzt = 1.5E-2 * ones(nstages,finalTime);
xPzt = 0.0    * ones(nstages,finalTime);
xTzt = 0.0    * ones(nstages,finalTime);
xHzt = 1.0 - xBzt - xPzt;

Tzt = 350 * ones(nstages,finalTime);

Azt =  zeros(nstages,finalTime);
Azt((1/5*nstages+1):(4/5*nstages),1:finalTime) = 1;
Xz0 = [xBzt(:,1);    % NSTAGE-long vector of XBs
       xHzt(:,1);    % NSTAGE-long vector of XHs
       xPzt(:,1);    % NSTAGE-long vector of XPs
       xTzt(:,1);    % NSTAGE-long vector of xTs
       Tzt(:,1);     % NSTAGE-long vector of tTs
       Azt(:,1)      % NSTAGE-long vector of xAs
      ];
  
% Inputs 
Feed = 2.5 * 10^-5;                              % (m^3 / sec) Volumetric flowrate
ut   = (Feed/(A*epz))*ones(1,finalTime);         % (m/sec)     Linear velocity    TF-long vector of Us   (disturbance/control)

xB0t = 1 * xBzt(1,1) * ones(1,finalTime);        % ()          Feed xB            TF-long vector of xB0  (disturbance)
xP0t = 1 * xPzt(1,1) * ones(1,finalTime);        % ()          Feed xP            TF-long vector of xP0  (disturbance)
xT0t = 6.36E-4 + xTzt(1,1) * ones(1,finalTime);        % ()          Feed xT            TF-long vector of xT0  (disturbance)
xH0t = 1 - xB0t - xT0t;                   % ()          Feed xH            TF long vector of XH0  (disturbance)

T0t  = 1 * Tzt(1,1)  * ones(1,finalTime);        % (C)         Feed T             TF-long vector of T0s  (control)
Ut   = [xB0t', xH0t', xP0t', xT0t', T0t', ut'];  %             Put stuff together
  
[yout,~,xout] = simulateSystem(F, g, Xz0, Ut);
    
vars = {'xB','xH2','xP','xT','Temp','Theta'};
for ip = 1:6
 subplot(3,2,ip); plot(linspace(0,1,nstages),(xout((ip-1)*nstages+1:ip*nstages,1:500:end))); 
 if ip == 6; xlabel('length'); end
 axis tight; title(vars(ip))
end
figure
for ip = 1:6
    subplot(3,2,ip); plot((yout(ip,:))'); 
    if ip == 8; xlabel('time'); end
    if ip == 1
        axis tight; title('Inlet temperature')
    elseif ip == 6
        axis tight; title('Exit temperature')
    else
        axis tight; title(sprintf('Catalyst section %d of 4', (ip-1)))
    end
end




