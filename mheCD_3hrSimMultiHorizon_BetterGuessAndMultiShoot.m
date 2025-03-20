function mheCD_3hrSimMultiHorizon_BetterGuessAndMultiShoot()
    %
    % Demonstrates MHE for a PDE with:
    %  - PDE states: xB, xH, xP, xT, T, Theta in 5 stages  => Nx=30
    %  - Only T in each stage + xB(final stage) measured
    %  - Horizon = 30 steps, each with 10 substeps (multi-shooting)
    %  - Forward-sim "open-loop" guess for the entire horizon
    %  - Optionally add node priors to anchor each horizon state
    
    %Add Casadi Path as needed
    import casadi.*
    
    %% 1) Basic PDE/Simulation Setup
    nstages   = 5;             
    Nx        = 6 * nstages;   % xB,xH,xP,xT,T,Theta per stage
    Nsteps    = 10800;         % 3 hours
    dt_main   = 1.0;           
    nSubSteps = 10;            % 10 sub-steps => matches the "true" PDE simulation
    dt_sub    = dt_main / nSubSteps;

    timeVec   = 0:Nsteps;

    %% 2) PDE Parameters (unchanged from your script)
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

    paramVal_num = [k0; K0; Q; E; rhoCp; negDeltaH; ...
                    rhoB; epz; MT; k0d; Ed_; P; L; r_unused; DeB];

    %% 3) Build symbolic PDE function
    x_sym = SX.sym('x', Nx,1);
    u_sym = SX.sym('u', 6,1);     % (xB_in, xH_in, xP_in, xT_in, T_in, flow)
    p_sym = SX.sym('p',15,1);

    xdot_expr = modelCD_pde_symbolic(x_sym, u_sym, p_sym, nstages);
    fPDE      = Function('fPDE',{x_sym,u_sym,p_sym},{xdot_expr});

    %% 4) "True" PDE Simulation (with sub-stepping)
    % constant input
    u_val = [0.02; 0.02; 0; 1e-4; 380; 0.1];

    % Initial condition
    xTrue_0 = zeros(Nx,1);
    for i=1:nstages
        xTrue_0(i)                = 0.01;          % xB
        xTrue_0(nstages + i)      = 0.01;          % xH
        xTrue_0(2*nstages + i)    = 0.01;          % xP
        xTrue_0(3*nstages + i)    = 0.01;          % xT
        xTrue_0(4*nstages + i)    = 400 + 5*(i-1); % T
        xTrue_0(5*nstages + i)    = 0.95;          % Theta. Can change to 1, depends on what you consider.
    end

    xTrue_all = zeros(Nx, Nsteps+1);
    xTrue_all(:,1) = xTrue_0;

    % Generate measurements: T in each stage + xB in final stage
    nMeas = nstages + 1; 
    meas_all = zeros(nMeas, Nsteps+1);

    noise_std_T   = 1.0;   % temperature noise std
    noise_frac_xB = 0.03;  % 3% noise for xB final

    for i=1:nstages
        meas_all(i,1) = xTrue_0(4*nstages + i);
    end
    meas_all(nstages+1,1) = xTrue_0(nstages);

    for k=1:Nsteps
        xkTrue = xTrue_all(:,k);
        xTemp  = xkTrue;
        for ssub=1:nSubSteps
            xdot_num = full(fPDE(xTemp, u_val, paramVal_num));
            xTemp    = xTemp + dt_sub*xdot_num;
            xTemp    = clampStateNumerically(xTemp, nstages);
        end
        xTrue_next = xTemp;
        xTrue_all(:,k+1) = xTrue_next;

        % measurement with noise
        for i=1:nstages
            Ttrue_i = xTrue_next(4*nstages + i);
            meas_all(i,k+1) = Ttrue_i + noise_std_T*randn;
        end
        xBfinal_true = xTrue_next(nstages);
        xB_sig = noise_frac_xB*xBfinal_true;
        if xB_sig<1e-10, xB_sig=1e-10; end
        meas_all(nstages+1,k+1) = xBfinal_true + xB_sig*randn;
    end

    %% 5) MHE Setup
    Nwin = 30;  % horizon length (30 steps, each 1 s => 30 s window)
    Xbig   = SX.sym('Xbig', Nx*(Nwin+1),1);  % stacked states over horizon
    Xprev  = SX.sym('Xprev', Nx,1);          % last known estimate, for prior
    Ubig   = SX.sym('Ubig',  6*(Nwin+1),1);  % inputs over horizon
    Meas   = SX.sym('Meas',  nMeas*(Nwin+1),1);
    Param  = SX.sym('Param', 15,1);

    cost_sym = SX(0);
    g_sym    = SX([]);

    % Weights
    wT         = 100;
    wB         = 500;
    wPrior     = 50;    % penalty on X^0 - Xprev
    wSmooth    = 5;     % penalty on consecutive states
    wNodePrior = 5;     % optional penalty for each X^i vs. open-loop guess

    %  pass in the "open-loop guess" for each node as. 
    XOL = SX.sym('XOL', Nx*(Nwin+1),1);  % open-loop guess for each node

    NxBlock = Nx;
    for iWin=0:Nwin
        % x_i = states at node i
        x_i = Xbig(iWin*NxBlock+1 : (iWin+1)*NxBlock);

        % Measurement offset
        measOff  = iWin*nMeas;
        T_meas_i = Meas(measOff + (1:nstages));
        xB_meas  = Meas(measOff + nstages+1);

        % PDE states
        T_sym   = x_i(4*nstages+1 : 5*nstages);
        xBfinal = x_i(nstages);

        % Measurement mismatch cost
        cost_sym = cost_sym ...
                   + wT*sumsqr(T_meas_i - T_sym) ...
                   + wB*( xB_meas - xBfinal )^2;

        % Node prior
        xOL_i = XOL(iWin*NxBlock+1 : (iWin+1)*NxBlock);
        cost_sym = cost_sym + wNodePrior*sumsqr(x_i - xOL_i);

        % composition sum constraint => xB+xH+xP+xT <=1
        xB_ = x_i(1:nstages);
        xH_ = x_i(nstages+1 : 2*nstages);
        xP_ = x_i(2*nstages+1 : 3*nstages);
        xT_ = x_i(3*nstages+1 : 4*nstages);
        for s=1:nstages
            sumComp = xB_(s) + xH_(s) + xP_(s) + xT_(s);
            g_sym = [g_sym; sumComp - 1];  % <=0
        end

        % PDE continuity constraints via multi-shooting
        if iWin < Nwin
            x_next = Xbig((iWin+1)*NxBlock+1 : (iWin+2)*NxBlock);
            u_i    = Ubig(iWin*6+1 : iWin*6+6);

            % Multi-step integration for 1 second using sub-steps
            xPred = multiStepIntegrationCasadi(x_i, u_i, Param, dt_sub, nSubSteps);

            % PDE eq => x_next - xPred= 0
            g_sym = [g_sym; x_next - xPred];

            % smoothing cost => penalize big jumps
            cost_sym = cost_sym + wSmooth*sumsqr(x_next - x_i);
        end
    end

    % prior on first node: keep near Xprev
    x0Block = Xbig(1:Nx);
    cost_sym = cost_sym + wPrior*sumsqr(x0Block - Xprev);

    % Build NLP
    % gather all "free" variables as Xbig, plus we pass Ubig, Meas, Param, XOL as "parameters"
    nlp = struct('x',Xbig, 'f',cost_sym, 'g',g_sym,...
                 'p',[Xprev; Ubig; Meas; Param; XOL]);
    opts=struct; opts.ipopt.print_level=0;
    MHEsolver = nlpsol('MHEsolver','ipopt', nlp, opts);

    % Constraints
    %  - nstages => comp sum <=0 for each node => total (Nwin+1)*nstages
    %  - PDE eq => Nx for each of (Nwin) transitions => Nx*Nwin
    g_lAll=[]; g_uAll=[];
    for iWin=0:Nwin
        g_lAll=[g_lAll; -inf*ones(nstages,1)]; 
        g_uAll=[g_uAll; zeros(nstages,1)];    % sumComp -1 <= 0
        if iWin<Nwin
            g_lAll=[g_lAll; zeros(Nx,1)];
            g_uAll=[g_uAll; zeros(Nx,1)];
        end
    end

    % Bounds for states
    lbX = -inf*ones(Nx*(Nwin+1),1);
    ubX =  inf*ones(Nx*(Nwin+1),1);
    for j=0:Nwin
        jOff = j*Nx;
        % xB,xH,xP,xT in [1e-10, 1]
        lbX(jOff+(1:4*nstages)) = 1e-10;
        ubX(jOff+(1:4*nstages)) = 1;
        % T => >=320
        lbX(jOff+(4*nstages+1 : 5*nstages)) = 320;
        % Theta => [0,1]
        lbX(jOff+(5*nstages+1 : 6*nstages)) = 0;
        ubX(jOff+(5*nstages+1 : 6*nstages)) = 1;
    end

    %% 6) Solve MHE step-by-step
    xEst_all = zeros(Nx, Nsteps+1);
    xEst_all(:,1) = xTrue_0;

    % store the horizon solution for warm start
    % PDE forward to build an open-loop guess each time.
    Xsol_horizon = repmat(xTrue_0, (Nwin+1),1);

    for k=1:Nsteps

        Xprev_k = xEst_all(:,k);

        % Gather horizon inputs & measurements
        U_window = [];
        M_window = [];
        for iWin=k:(k+Nwin)
            iUse = min(iWin, Nsteps);
            U_window = [U_window; u_val];
            M_window = [M_window; meas_all(:, iUse+1)];
        end

        %==== Build an open-loop guess for the horizon (XOL) ====
        XOL_guess = buildHorizonGuess(Xprev_k, Nwin, paramVal_num, ...
                                      dt_sub, nSubSteps, fPDE, u_val, nstages);

        %==== Also build the initial Xbig from that same open-loop guess ====
        % I tried a start from last horizon, with:
        %  Xinit = shiftHorizon(Xsol_horizon, Nx, Nwin);
        % But for this sim use the open-loop guess directly:
        Xinit = XOL_guess;  % same dimension Nx*(Nwin+1)

        % pack the param vector
        pVal_k = [Xprev_k; U_window; M_window; paramVal_num; XOL_guess];

        sol_k  = MHEsolver('x0',  Xinit, ...
                           'lbx', lbX, 'ubx', ubX, ...
                           'lbg', g_lAll, 'ubg', g_uAll, ...
                           'p',   pVal_k);

        X_sol  = full(sol_k.x);
        Xsol_horizon = X_sol;  % store for next iteration

        xEst_kplus1 = X_sol(1:Nx);  % first block in solution is new estimate
        xEst_all(:,k+1) = xEst_kplus1;

        if mod(k,600)==0 || k==Nsteps
            fprintf('Step %4d/%4d => T5_est=%.2f, xB5=%.3f, Th5=%.3f\n',...
                k, Nsteps, xEst_kplus1(4*nstages+5), ...
                xEst_kplus1(nstages), xEst_kplus1(5*nstages+5));
        end
    end

    %% 7) Plots
    figure('Name','MHECD_MultiHorizon_MultiShoot','Color','white');
    stageColors = lines(nstages);

    % xB
    subplot(2,3,1); hold on; grid on;
    plotVar(xTrue_all(1:nstages,:), xEst_all(1:nstages,:), timeVec, stageColors, 'x_B');

    % xH
    subplot(2,3,2); hold on; grid on;
    plotVar(xTrue_all(nstages+(1:nstages),:), xEst_all(nstages+(1:nstages),:), timeVec, stageColors, 'x_H');

    % xP
    subplot(2,3,3); hold on; grid on;
    plotVar(xTrue_all(2*nstages+(1:nstages),:), xEst_all(2*nstages+(1:nstages),:), timeVec, stageColors, 'x_P');

    % xT
    subplot(2,3,4); hold on; grid on;
    plotVar(xTrue_all(3*nstages+(1:nstages),:), xEst_all(3*nstages+(1:nstages),:), timeVec, stageColors, 'x_T');

    % T
    subplot(2,3,5); hold on; grid on;
    plotVar_T_only( xTrue_all(4*nstages+(1:nstages),:), ...
                    xEst_all(4*nstages+(1:nstages),:), ...
                    meas_all(1:nstages,:), ...
                    timeVec, stageColors );
    title('Temperature in Each Stage');

    % placeholder in subplot(2,3,6) (or leave blank)
    subplot(2,3,6); axis off;
    text(0.1,0.5,'\Theta is in separate figure','FontSize',10);

    % Now a separate figure for Catalyst (Theta)
    figure('Name','Catalyst_Theta_Only','Color','white');
    hold on; grid on;
    stride = 500;
    idxPlot= 1:stride:length(timeVec);
    for st=1:nstages
        ThTrue_st = xTrue_all(5*nstages + st, :);
        ThEst_st  = xEst_all(5*nstages + st, :);
        plot(timeVec(idxPlot), ThTrue_st(idxPlot), 'o-',...
            'LineWidth',1.2, 'MarkerSize',4, 'Color',stageColors(st,:),...
            'DisplayName', sprintf('\\Theta True (stage %d)', st));
        plot(timeVec(idxPlot), ThEst_st(idxPlot), 'x--',...
            'LineWidth',1.2, 'MarkerSize',4, 'Color',stageColors(st,:),...
            'HandleVisibility','off');
    end
    xlabel('Time [s]');
    ylabel('\Theta');
    legend('Location','best','FontSize',8);
    title('Catalyst Decay (Theta): True vs. MHE Estimate');
end

%% ========================================================================
%%  A) Multi-step integration in CasADi => 1 second with nSubSteps
function xNextSym = multiStepIntegrationCasadi(xCurSym, uSym, paramSym, dt_sub, nSubSteps)
    % multiStepIntegrationCasadi:
    %   xNext = xCur + sum_{sub=1..nSubSteps} [ dt_sub * fPDE(...) ],
    %   done in a loop of symbolic expressions
    import casadi.*
    xTemp = xCurSym;
    fPDEfun = @modelCD_pde_symbolic;  
    nstages = 5;  % be consistent with main file. Change if needed!

    for ss=1:nSubSteps
        xdot_sub = fPDEfun(xTemp, uSym, paramSym, nstages);
        xTemp    = xTemp + dt_sub*xdot_sub;
    end
    xNextSym = xTemp;
end

%% ========================================================================
%%  B) Build Open-Loop Guess for the Horizon
function XOL_guess = buildHorizonGuess(xStart, Nwin, paramVal, dt_sub, nSubSteps, fPDE, u_val, nstages)
    % Repeatedly integrate PDE forward 1 second at a time (with 10 sub-steps)
    % from xStart for Nwin steps. This yields Nx*(Nwin+1) states.
    Nx = length(xStart);
    XOL_guess = zeros(Nx*(Nwin+1),1);

    xNow = xStart;
    for i=0:Nwin
        XOL_guess(i*Nx+1 : (i+1)*Nx) = xNow;
        if i < Nwin
            % Integrate forward 1 second
            for ss=1:nSubSteps
                xdot_now = full(fPDE(xNow, u_val, paramVal));
                xNow     = xNow + dt_sub*xdot_now;
                xNow     = clampStateNumerically(xNow, nstages);
            end
        end
    end
end

%% ========================================================================
%%  C) PDE model (same as your script)
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

    % "padding" for T
    Tpad = SX.zeros(nstages+1,1);
    Tpad(1) = T_in;
    for j=1:nstages
        Tpad(j+1) = T_(j);
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
        advEn  = -flow*(P/(R*T_cur))/rhoCp * cMixVal*(T_cur - T_up)/dz;
        heatTerm= (negDeltaH/rhoCp)*(P/(R*T_cur))*rH_*(rhoB/epz);
        xdot(iT)= advEn + heatTerm;

        iTh= 5*nstages + i;
        xdot(iTh)= -rD_;
    end

    function val= clampExp(e)
        val= if_else(e>50, 50, e);
        val= if_else(val<-50, -50, val);
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

%% ========================================================================
%%  D) Numeric clamp 
function xC = clampStateNumerically(xRaw, nstages)
    xC = xRaw;
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

%% ========================================================================
%%  E)Plotting
function plotVar(xTrueMat, xEstMat, tvec, stageColors, labelStr)
    stride  = 500;
    idxPlot = 1:stride:length(tvec);
    nstages = size(xTrueMat,1);
    for i=1:nstages
        plot(tvec(idxPlot), xTrueMat(i, idxPlot), 'o-', ...
            'LineWidth',1,'MarkerSize',4,'Color',stageColors(i,:), ...
            'DisplayName', sprintf('%s True(stg %d)', labelStr, i));
        plot(tvec(idxPlot), xEstMat(i, idxPlot), 'x--', ...
            'LineWidth',1,'MarkerSize',4,'Color',stageColors(i,:),...
            'HandleVisibility','off');
    end
    xlabel('Time [s]');
    ylabel(labelStr);
    legend('Location','best','FontSize',8);
    grid on;
end

function plotVar_T_only(TtrueMat, TestMat, TmeasMat, tvec, stageColors)
    stride  = 500;
    idxPlot = 1:stride:length(tvec);
    nstages = size(TtrueMat,1);
    for i=1:nstages
        plot(tvec(idxPlot), TtrueMat(i, idxPlot), 'o-', ...
            'LineWidth',1,'MarkerSize',4,'Color',stageColors(i,:), ...
            'DisplayName', sprintf('T True (stg %d)', i));
        plot(tvec(idxPlot), TestMat(i, idxPlot), 'x--',...
            'LineWidth',1,'MarkerSize',4,'Color',stageColors(i,:),...
            'HandleVisibility','off');
        plot(tvec(idxPlot), TmeasMat(i, idxPlot), 's',...
            'LineWidth',1,'MarkerSize',4,'Color',stageColors(i,:),...
            'HandleVisibility','off');
    end
    xlabel('Time [s]');
    ylabel('Temperature [K]');
    legend('Location','best','FontSize',8);
    grid on;
end
