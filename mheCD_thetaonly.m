function mheCD_thetaOnly_cleanedUp()
    % mheCD_thetaOnly_cleanedUp
    %
    % Approach #1: "Theta-only" MHE
    %
    % - 5 PDE stages, 6 states each => 30 total states (xB,xH,xP,xT,T,Th)
    % - We measure xB,xH,xP,xT,T in each stage => treat them as "known/experimental"
    %   from data (with ±2% noise).
    %
    % Plots:
    %   1) PDE_True_Dynamics: entire PDE simulation for all 6 states
    %   2) MHE_ThetaOnly_main: main figure with xB..T subplots (true vs. noisy measurements)
    %   PLot 2 doesn't work/don't need.
    %   3) MHE_ThetaOnly_Catalyst: separate large figure for Theta only (true vs. MHE).

    % 1) Add your CasADi path
    addpath('/Users/jessenyman/Documents/casadi-3')

    import casadi.*

    %% A) Setup
    nstages   = 5;                 % PDE segments
    Nx_full   = 6*nstages;         % 6 states/stage => 30 total
    Nsteps    = 10800;             % main time steps (could reduce for testing)
    dt_main   = 1.0;
    nSubSteps = 10;
    dt_sub    = dt_main / nSubSteps;

    timeVec   = 0:Nsteps;

    % PDE parameters
    k0        = 3.16e-2;
    K0        = 3.16e-13;
    Q         = -6.89e7;
    E         =  5.76e7;
    rhoCp     =  7.322e5;
    negDeltaH =  2.09e8;
    rhoB      =  4.14e2;
    epz       =  0.58;
    MT        =  4.05e-4;
    k0d       =  1.80e-4;
    Ed_       =  4.52e6;
    P         =  1.01e5;
    L         =  0.118*5;
    r_unused  =  0.0078; 
    DeB       =  4.5e-5;

    paramVal_num = [k0;K0;Q;E;rhoCp;negDeltaH;rhoB;epz;MT; ...
                    k0d;Ed_;P;L;r_unused;DeB];

    % PDE function for "true" sim
    x_sym = SX.sym('x', Nx_full,1);
    u_sym = SX.sym('u', 6,1);
    p_sym = SX.sym('p',15,1);
    xdot_expr = modelCD_pde_symbolic(x_sym, u_sym, p_sym, nstages);
    fPDE      = Function('fPDE',{x_sym,u_sym,p_sym},{xdot_expr});

    % fix constant inputs; could be variable in the future
    u_val = [0.02; 0.02; 0; 1e-4; 380; 0.1];

    %% B) Generate "True" PDE + Synthetic Data
    xTrue_0 = zeros(Nx_full,1);
    for i=1:nstages
        xTrue_0(i)                 = 0.01;            % xB
        xTrue_0(nstages + i)       = 0.01;            % xH
        xTrue_0(2*nstages + i)     = 0.01;            % xP
        xTrue_0(3*nstages + i)     = 0.01;            % xT
        xTrue_0(4*nstages + i)     = 400 + 5*(i-1);   % T
        xTrue_0(5*nstages + i)     = 0.95;            % Th
    end

    xTrue_all = zeros(Nx_full, Nsteps+1);
    xTrue_all(:,1) = xTrue_0;

    % measure xB,xH,xP,xT,T => 5 per stage => 25 total
    nMeas   = 5*nstages;
    meas_all= zeros(nMeas, Nsteps+1);

    noiseFrac = 0.02;  % ±2% noise

    for k=1:Nsteps
        xk = xTrue_all(:,k);
        xTemp = xk;
        % sub-step PDE
        for ss=1:nSubSteps
            xdot_num= full(fPDE(xTemp, u_val, paramVal_num));
            xTemp   = xTemp + dt_sub*xdot_num;
            xTemp   = clampStateWithTh(xTemp, nstages);
        end
        xTrue_next= xTemp;
        xTrue_all(:, k+1)= xTrue_next;

        % build measurements => xB,xH,xP,xT,T only
        idxM=1;
        for st=1:nstages
            xB_st= xTrue_next(st);
            xH_st= xTrue_next(nstages + st);
            xP_st= xTrue_next(2*nstages + st);
            xT_st= xTrue_next(3*nstages + st);
            TT_st= xTrue_next(4*nstages + st);

            meas_all(idxM,   k+1)= xB_st*(1+ noiseFrac*(2*rand-1));
            meas_all(idxM+1, k+1)= xH_st*(1+ noiseFrac*(2*rand-1));
            meas_all(idxM+2, k+1)= xP_st*(1+ noiseFrac*(2*rand-1));
            meas_all(idxM+3, k+1)= xT_st*(1+ noiseFrac*(2*rand-1));
            meas_all(idxM+4, k+1)= TT_st*(1+ noiseFrac*(2*rand-1));
            idxM= idxM+5;
        end
    end

    %% B1) Optional Plot => PDE_True_Dynamics (all states)
    % PDE evolution for xB..Th
    plotPDEonlyDynamics(xTrue_all, nstages, timeVec);

    %% C) "Theta-only" MHE
    Nx_Th= nstages;
    ThEst_all= zeros(Nx_Th, Nsteps+1);
    % initial guess
    ThEst_all(:,1)= xTrue_0(5*nstages+(1:nstages));

    for k=1:Nsteps
        Th_prev= ThEst_all(:,k);

        % read the measured states for time k+1
        measVec_k= meas_all(:,k+1);
        [xB_, xH_, xP_, xT_, T_] = extractMeasuredStates(measVec_k, nstages);

        % PDE mismatch => Th_next ~ Th_prev + dt*( -rDeact( xT, T, Th_prev ) )
        Th_next_sym= SX.sym('Th_next', Nx_Th,1);
        cost_k= SX(0);

        dt_approx= 1.0;    % each main step is 1s
        wPDE    = 1000;    % PDE mismatch weight
        wPrior  = 10;      % prior cost on Th

        for st=1:nstages
            rD_val= rxnDeact( xT_(st), T_(st), Th_prev(st), paramVal_num );
            Th_pred= Th_prev(st) + dt_approx*( -rD_val );
            mismatch= Th_next_sym(st) - Th_pred;
            cost_k= cost_k + wPDE*(mismatch^2);
        end

        % prior => stay near Th_prev
        cost_k= cost_k + wPrior*sumsqr( Th_next_sym - Th_prev );

        nlp_struct= struct('x', Th_next_sym, 'f', cost_k, 'g',[], 'p',[]);
        solver_opts= struct; solver_opts.ipopt.print_level=0;
        solver_k= nlpsol('ThetaMHE','ipopt', nlp_struct, solver_opts);

        solK= solver_k('x0', Th_prev);
        Th_next= full(solK.x);

        ThEst_all(:,k+1)= Th_next;
    end

    %% D) Plot main figure => measured states xB..T in subplots (true vs. data).
    figure('Name','MHE_ThetaOnly_main','Color','white');
    stageColors= defineCustomColors(nstages);
    stride=500;  % large step downsample for 10800 steps

    subplot(2,3,1); hold on; grid on;
    plotState_measVsTrue('x_B', 1, xTrue_all, meas_all, timeVec, nstages, stageColors, stride);

    subplot(2,3,2); hold on; grid on;
    plotState_measVsTrue('x_H', 2, xTrue_all, meas_all, timeVec, nstages, stageColors, stride);

    subplot(2,3,3); hold on; grid on;
    plotState_measVsTrue('x_P', 3, xTrue_all, meas_all, timeVec, nstages, stageColors, stride);

    subplot(2,3,4); hold on; grid on;
    plotState_measVsTrue('x_T', 4, xTrue_all, meas_all, timeVec, nstages, stageColors, stride);

    subplot(2,3,5); hold on; grid on;
    plotState_measVsTrue('T',   5, xTrue_all, meas_all, timeVec, nstages, stageColors, stride);

    subplot(2,3,6); axis off;
    text(0.1,0.5,'\Theta is unmeasured. See separate figure','FontSize',10);

    %% E) Separate figure => Catalyst (Theta) alone
    figure('Name','MHE_ThetaOnly_Catalyst','Color','white');
    hold on; grid on;
    idxPlot= 1:stride:length(timeVec);

    for i=1:nstages
        ThTrue_i= xTrue_all(5*nstages + i,:);
        ThEst_i = ThEst_all(i,:);
        plot(timeVec(idxPlot), ThTrue_i(idxPlot), 'o-', ...
            'Color', stageColors(i,:), 'LineWidth',1.5, 'MarkerSize',5,...
            'DisplayName', sprintf('\\Theta True (stage %d)', i));
        plot(timeVec(idxPlot), ThEst_i(idxPlot), 'x--', ...
            'Color', stageColors(i,:), 'HandleVisibility','off');
    end
    xlabel('Time'); ylabel('\Theta');
    legend('Location','best','FontSize',10);
    title('\Theta (unmeasured) : Approach #1 MHE', 'FontWeight','bold');
end

%% ======================================================================
%% PDE model for "true" sim
function xdot = modelCD_pde_symbolic(x, u, param, nstages)
    import casadi.*

    k0        = param(1);
    K0        = param(2);
    Q         = param(3);
    E         = param(4);
    rhoCp     = param(5);
    negDeltaH = param(6);
    rhoB      = param(7);
    epz       = param(8);
    MT        = param(9);
    k0d       = param(10);
    Ed_       = param(11);
    P         = param(12);
    L         = param(13);
   
    DeB       = param(15);

    R  = 8.314e3;
    dz = L / nstages;

    xB_ = x(1:nstages);
    xH_ = x(nstages+1 : 2*nstages);
    xP_ = x(2*nstages+1 : 3*nstages);
    xT_ = x(3*nstages+1 : 4*nstages);
    T_  = x(4*nstages+1 : 5*nstages);
    Th_ = x(5*nstages+1 : 6*nstages);

    xB_in = u(1);
    xH_in = u(2);
    xP_in = u(3);
    xT_in = u(4);
    T_in  = u(5);
    flow  = u(6);

    xdot = SX.zeros(6*nstages,1);

    Tpad= SX.zeros(nstages+1,1);
    Tpad(1)= T_in;
    for j=1:nstages
        Tpad(j+1)= T_(j);
    end

    for i=1:nstages
        rH_ = rxnHydro(i);
        rD_ = rxnDeact(i);

        [dB,aB] = diffAdv(i, xB_, xB_in);
        xdot(i) = dB - aB - (rhoB/epz)*rH_;

        iH= nstages + i;
        [dH,aH] = diffAdv(i, xH_, xH_in);
        xdot(iH)= dH - aH - (rhoB/epz)*rH_;

        iP= 2*nstages + i;
        [dP,aP] = diffAdv(i, xP_, xP_in);
        xdot(iP)= dP - aP + (rhoB/epz)*rH_;

        iT_=3*nstages + i;
        [dT,aT] = diffAdv(i, xT_, xT_in);
        consumeT= (rhoB/epz)*MT*rD_;
        xdot(iT_)= dT - aT - consumeT;

        iT= 4*nstages + i;
        T_up  = Tpad(i);
        T_cur = Tpad(i+1);
        cMixVal= (2.902*xH_(i)+9.686*xB_(i))*1e4;
        advEn  = - flow*(P/(R*T_cur))/rhoCp * cMixVal*(T_cur - T_up)/dz;
        heatTerm= (negDeltaH/rhoCp)*(P/(R*T_cur))*rH_*(rhoB/epz);
        xdot(iT)= advEn + heatTerm;

        iTh= 5*nstages + i;
        xdot(iTh)= -rD_;
    end

    function val= clampExp(e)
        val= if_else(e>50, 50, e);
        val= if_else(val<-50,-50,val);
    end
    function rHval= rxnHydro(idx)
        argQE= clampExp(-(Q+E)/(R*T_(idx)));
        argQ = clampExp(-Q/(R*T_(idx)));
        denom= 1 + K0*exp(argQ)*(P*xB_(idx));
        rHval= k0*K0*exp(argQE)*(P*xB_(idx))*(P*xH_(idx))/denom * Th_(idx);
    end
    function rDval= rxnDeact(idx)
        argD= clampExp(-Ed_/(R*T_(idx)));
        rDval= k0d*exp(argD)*(P*xT_(idx))*Th_(idx);
    end
    function [dTerm,aTerm] = diffAdv(i_, xVec, xInVal)
        if nstages>1
            if i_==1
                dTerm= DeB/(dz^2)*( xVec(2)-2*xVec(1)+ xInVal );
                aTerm= (flow/dz)*( xVec(1)- xInVal );
            elseif i_==nstages
                dTerm= DeB/(dz^2)*( xVec(nstages-1)-2*xVec(nstages)+ xVec(nstages-1) );
                aTerm= (flow/dz)*( xVec(nstages)- xVec(nstages-1) );
            else
                dTerm= DeB/(dz^2)*( xVec(i_+1)-2*xVec(i_)+ xVec(i_-1) );
                aTerm= (flow/dz)*( xVec(i_)- xVec(i_-1) );
            end
        else
            dTerm= SX(0);
            aTerm= (flow/dz)*( xVec(i_)- xInVal );
        end
    end
end

%% ======================================================================
%% PDEonly dynamics => plot xB,xH,xP,xT,T,Th subplots
function plotPDEonlyDynamics(xTrue_all, nstages, timeVec)
    figure('Name','PDE_True_Dynamics','Color','white');
    stageColors= defineCustomColors(nstages);
    stride= max(round(length(timeVec)/200),1);
    idxPlot= 1:stride:length(timeVec);

    % xB
    subplot(2,3,1); hold on; grid on;
    plotVar_simple(xTrue_all(1:nstages,:), timeVec, idxPlot, stageColors, 'x_B');

    % xH
    subplot(2,3,2); hold on; grid on;
    plotVar_simple(xTrue_all(nstages+(1:nstages),:), timeVec, idxPlot, stageColors, 'x_H');

    % xP
    subplot(2,3,3); hold on; grid on;
    plotVar_simple(xTrue_all(2*nstages+(1:nstages),:), timeVec, idxPlot, stageColors, 'x_P');

    % xT
    subplot(2,3,4); hold on; grid on;
    plotVar_simple(xTrue_all(3*nstages+(1:nstages),:), timeVec, idxPlot, stageColors, 'x_T');

    % T
    subplot(2,3,5); hold on; grid on;
    plotVar_simple(xTrue_all(4*nstages+(1:nstages),:), timeVec, idxPlot, stageColors, 'T');

    % Th
    subplot(2,3,6); hold on; grid on;
    plotVar_simple(xTrue_all(5*nstages+(1:nstages),:), timeVec, idxPlot, stageColors, '\Theta');
end

function plotVar_simple(xMat, timeVec, idxPlot, colors, labelStr)
    nstages= size(xMat,1);
    for i=1:nstages
        plot(timeVec(idxPlot), xMat(i, idxPlot), 'o-',...
            'MarkerSize',4,'LineWidth',1,'Color',colors(i,:),...
            'DisplayName', sprintf('%s True (st%d)', labelStr, i));
    end
    xlabel('Time'); ylabel(labelStr);
    legend('Location','best','FontSize',8);
    title(sprintf('%s in %d stages', labelStr, nstages));
end

%% ======================================================================
%% Helper to define custom stage colors
function cMap = defineCustomColors(nstages)
    baseColors = [
        0.00,0.45,0.70;  % dark blue
        0.80,0.00,0.00;  % red
        0.47,0.67,0.19;  % green
        0.93,0.69,0.13;  % orange
        0.49,0.18,0.56;  % purple
        0.30,0.75,0.93;  % skyblue
        0.64,0.08,0.08;  % dark red
        0.66,0.66,0.66;  % grey
    ];
    cMap= baseColors(1:nstages,:);
end

%% ======================================================================
%% clamp
function xC= clampStateWithTh(xRaw, nstages)
    xC= xRaw;
    for j=1:4*nstages
        if xC(j)<1e-10, xC(j)=1e-10; end
        if xC(j)>1, xC(j)=1; end
    end
    for j=4*nstages+1 : 5*nstages
        if xC(j)<320, xC(j)=320; end
    end
    for j=5*nstages+1 : 6*nstages
        if xC(j)<0, xC(j)=0; end
        if xC(j)>1, xC(j)=1; end
    end
end

%% ======================================================================
%% read out xB,xH,xP,xT,T from measVec
function [xB_, xH_, xP_, xT_, T_] = extractMeasuredStates(measVec, nstages)
    xB_ = zeros(nstages,1);
    xH_ = zeros(nstages,1);
    xP_ = zeros(nstages,1);
    xT_ = zeros(nstages,1);
    T_  = zeros(nstages,1);

    idx=1;
    for s=1:nstages
        xB_(s)= measVec(idx);   idx=idx+1;
        xH_(s)= measVec(idx);   idx=idx+1;
        xP_(s)= measVec(idx);   idx=idx+1;
        xT_(s)= measVec(idx);   idx=idx+1;
        T_(s) = measVec(idx);   idx=idx+1;
    end
end

%% ======================================================================
%% local rD => -rDeact
function rVal= rxnDeact(xTval, Tval, Thval, paramVal)
    % from PDE => rD= k0d*exp( -Ed/(R*T) )*(P*xT)* Th
    k0d= paramVal(10);
    Ed_= paramVal(11);
    P  = paramVal(12);
    R  = 8.314e3;

    argD= -(Ed_)/(R*Tval);
    if argD>50, argD=50; elseif argD<-50, argD=-50; end
    rVal= k0d*exp(argD)*(P*xTval)* Thval;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotState_measVsTrue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotState_measVsTrue(labelStr, indexBlock, xTrue_all, meas_all, ...
                              timeVec, nstages, stageColors, stride)

    idxPlot = 1:stride:length(timeVec);
    hold on;
    for st = 1:nstages
        % PDE "true" state for this variable in stage st
        xTrueVals = xTrue_all((indexBlock - 1)*nstages + st, :);

        % The measurement data is in meas_all:
        % each stage has 5 measurements => (st-1)*5 + indexBlock
        measIdx  = (st - 1)*5 + indexBlock;
        measVals = meas_all(measIdx, :);

        % Plot PDE "true" with circle markers
        plot(timeVec(idxPlot), xTrueVals(idxPlot), 'o-', ...
            'LineWidth', 1, ...
            'MarkerSize', 4, ...
            'Color', stageColors(st,:), ...
            'DisplayName', sprintf('%s True (stage %d)', labelStr, st));

        % Plot measurement data with square markers
        plot(timeVec(idxPlot), measVals(idxPlot), 's', ...
            'Color', stageColors(st,:), ...
            'HandleVisibility', 'off');
    end

    xlabel('Time'); 
    ylabel(labelStr);
    legend('Location','best','FontSize',8);
    title(sprintf('%s in %d stages', labelStr, nstages));
    grid on;
end
