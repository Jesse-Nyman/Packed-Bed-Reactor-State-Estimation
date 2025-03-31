function xdot = modelCD(~, x, u, theta)
    import casadi.*

    % Unpack parameters
    nstages   = theta{1};
    k0        = theta{2};
    K0        = theta{3};
    Q         = theta{4};
    E         = theta{5};
    rhoCp     = theta{6};
    negDeltaH = theta{7};
    rhoB      = theta{8};
    epz       = theta{9};
    MT        = theta{10};
    k0d       = theta{11};
    Ed_       = theta{12};
    P         = theta{13};
    L         = theta{14};
    r         = theta{15}; %#ok<NASGU>
    DeB       = theta{16};

    R         = 8.314e3;    
    deltaZ    = L/nstages;
    epsilon   = 1e-12;

    % For convenience
    i_mid = 2:nstages-1;

    % CasADi-safe max
    fmax_ = @(X,Y) if_else(X>Y,X,Y);

    % Extract states
    xB  = x(1:nstages);
    xH  = x(nstages+1:2*nstages);
    xP  = x(2*nstages+1:3*nstages);
    xT_ = x(3*nstages+1:4*nstages);
    T_  = x(4*nstages+1:5*nstages);
    Th  = x(5*nstages+1:6*nstages);

    % Enforce positivity
    xB  = fmax_(xB, epsilon);
    xH  = fmax_(xH, epsilon);
    xP  = fmax_(xP, epsilon);
    xT_ = fmax_(xT_, epsilon);
    T_  = fmax_(T_, 1e-3);
    Th  = fmax_(Th, epsilon);

    % Inputs
    xB_in = u(1);
    xH_in = u(2);
    xP_in = u(3);
    xT_in = u(4);
    T_in  = u(5);
    u_t   = u(6);

    expQ_E = @(Tval) exp((-Q - E)./(R.*Tval));
    denom = @(Tval,xBval) fmax_(1 + K0.*exp(-Q./(R.*Tval)).*P.*xBval, epsilon);

    % xBdot
    xBdot_1 = DeB*(xB(2)-2*xB(1)+xB_in)/deltaZ^2 ...
        - u_t*(xB(1)-xB_in)/deltaZ ...
        - (rhoB/(P*epz))*R.*T_(1).*(P.*xB(1)).*(P.*xH(1))*k0*K0.*expQ_E(T_(1))./denom(T_(1),xB(1)).*Th(1);

    xBdot_mid = DeB*(xB(i_mid+1)-2*xB(i_mid)+xB(i_mid-1))/deltaZ^2 ...
        - u_t*(xB(i_mid)-xB(i_mid-1))/deltaZ ...
        - (rhoB/(P*epz))*R.*T_(i_mid).*(P.*xB(i_mid)).*(P.*xH(i_mid))*k0*K0.*expQ_E(T_(i_mid))./denom(T_(i_mid),xB(i_mid)).*Th(i_mid);

    xBdot_end = DeB*(xB(nstages-1)-2*xB(nstages)+xB(nstages-2))/deltaZ^2 ...
        - u_t*(xB(nstages)-xB(nstages-1))/(2*deltaZ) ...
        - (rhoB/(P*epz))*R.*T_(nstages).*(P.*xB(nstages)).*(P.*xH(nstages))*k0*K0.*expQ_E(T_(nstages))./denom(T_(nstages),xB(nstages)).*Th(nstages);

    % xHdot
    xHdot_1 = DeB*(xH(2)-2*xH(1)+xH_in)/deltaZ^2 ...
        - u_t*(xH(1)-xH_in)/deltaZ ...
        - (rhoB/(P*epz))*R.*T_(1).*(P.*xB(1)).*(P.*xH(1))*k0*K0.*expQ_E(T_(1))./denom(T_(1),xB(1)).*Th(1);

    xHdot_mid = DeB*(xH(i_mid+1)-2*xH(i_mid)+xH(i_mid-1))/deltaZ^2 ...
        - u_t*(xH(i_mid)-xH(i_mid-1))/deltaZ ...
        - (rhoB/(P*epz))*R.*T_(i_mid).*(P.*xB(i_mid)).*(P.*xH(i_mid))*k0*K0.*expQ_E(T_(i_mid))./denom(T_(i_mid),xB(i_mid)).*Th(i_mid);

    xHdot_end = DeB*(xH(nstages-1)-2*xH(nstages)+xH(nstages-2))/deltaZ^2 ...
        - u_t*(xH(nstages)-xH(nstages-1))/deltaZ ...
        - (rhoB/(P*epz))*R.*T_(nstages).*(P.*xB(nstages)).*(P.*xH(nstages))*k0*K0.*expQ_E(T_(nstages))./denom(T_(nstages),xB(nstages)).*Th(nstages);

    % xPdot
    xPdot_1 = DeB*(xP(2)-2*xP(1)+xP_in)/deltaZ^2 ...
        - u_t*(xP(1)-xP_in)/deltaZ ...
        + (rhoB/(P*epz))*R.*T_(1).*(P.*xB(1)).*(P.*xH(1))*k0*K0.*expQ_E(T_(1))./denom(T_(1),xB(1)).*Th(1);

    xPdot_mid = DeB*(xP(i_mid+1)-2*xP(i_mid)+xP(i_mid-1))/deltaZ^2 ...
        - u_t*(xP(i_mid)-xP(i_mid-1))/deltaZ ...
        + (rhoB/(P*epz))*R.*T_(i_mid).*(P.*xB(i_mid)).*(P.*xH(i_mid))*k0*K0.*expQ_E(T_(i_mid))./denom(T_(i_mid),xB(i_mid)).*Th(i_mid);

    xPdot_end = DeB*(xP(nstages-1)-2*xP(nstages)+xP(nstages-2))/deltaZ^2 ...
        - u_t*(xP(nstages)-xP(nstages-1))/deltaZ ...
        + (rhoB/(P*epz))*R.*T_(nstages).*(P.*xB(nstages)).*(P.*xH(nstages))*k0*K0.*expQ_E(T_(nstages))./denom(T_(nstages),xB(nstages)).*Th(nstages);

    % xTdot
    xTdot_1 = DeB*(xT_(2)-2*xT_(1)+xT_in)/deltaZ^2 ...
        - u_t*(xT_(1)-xT_in)/deltaZ ...
        + (rhoB*R/(epz*P)).*T_(1).*MT.*P.*xT_(1)*(-k0d).*exp(-Ed_/(R.*T_(1))).*(Th(1).^2);

    xTdot_mid = DeB*(xT_(i_mid+1)-2*xT_(i_mid)+xT_(i_mid-1))/deltaZ^2 ...
        - u_t*(xT_(i_mid)-xT_(i_mid-1))/deltaZ ...
        + (rhoB*R/(epz*P)).*T_(i_mid).*MT.*P.*xT_(i_mid)*(-k0d).*exp(-Ed_./(R.*T_(i_mid))).*(Th(i_mid).^2);

    xTdot_end = DeB*(xT_(nstages-1)-2*xT_(nstages)+xT_(nstages-2))/deltaZ^2 ...
        - u_t*(xT_(nstages)-xT_(nstages-1))/deltaZ ...
        + (rhoB*R/(epz*P)).*T_(nstages).*MT.*P.*xT_(nstages)*(-k0d).*exp(-Ed_./(R.*T_(nstages))).*(Th(nstages).^2);

    % Tdot
    heatTerm_1 = (2.902.*xH(1)+9.686.*xB(1))*1e4;
    Tdot_1 = -u_t*(P/(R*T_(1)))/rhoCp*heatTerm_1*(T_(1)-T_in)/deltaZ ...
        + negDeltaH/rhoCp*(P/(R*T_(1)))*(rhoB/(P*epz))*R.*T_(1).*(P.*xB(1)).*(P.*xH(1))*k0*K0.*expQ_E(T_(1))./denom(T_(1),xB(1)).*Th(1);

    heatTerm_mid = (2.902.*xH(i_mid)+9.686.*xB(i_mid))*1e4;
    Tdot_mid = -u_t*(P./(R.*T_(i_mid)))./rhoCp.*heatTerm_mid.*(T_(i_mid)-T_(i_mid-1))/deltaZ ...
        + negDeltaH./rhoCp.*(P./(R.*T_(i_mid))).*(rhoB/(P*epz))*R.*T_(i_mid).*(P.*xB(i_mid)).*(P.*xH(i_mid))*k0*K0.*expQ_E(T_(i_mid))./denom(T_(i_mid),xB(i_mid)).*Th(i_mid);

    heatTerm_end = (2.902.*xH(nstages)+9.686.*xB(nstages))*1e4;
    Tdot_end = -u_t*(P/(R*T_(nstages)))/rhoCp*heatTerm_end*(T_(nstages)-T_(nstages-1))/deltaZ ...
        + negDeltaH/rhoCp*(P/(R*T_(nstages)))*(rhoB/(P*epz))*R.*T_(nstages).*(P.*xB(nstages)).*(P.*xH(nstages))*k0*K0.*expQ_E(T_(nstages))./denom(T_(nstages),xB(nstages)).*Th(nstages);

    % Thdot
    Thdot_1 = -k0d*P.*xT_(1).*Th(1).*exp(-Ed_/(R.*T_(1)));

    Thdot_mid = -k0d*P.*xT_(i_mid).*Th(i_mid).*exp(-Ed_./(R.*T_(i_mid)));

    % Add Thdot_end for the last stage
    Thdot_end = -k0d*P.*xT_(nstages).*Th(nstages).*exp(-Ed_./(R.*T_(nstages)));

    % Combine all derivatives
    xdot = [xBdot_1; xBdot_mid; xBdot_end;
            xHdot_1; xHdot_mid; xHdot_end;
            xPdot_1; xPdot_mid; xPdot_end;
            xTdot_1; xTdot_mid; xTdot_end;
            Tdot_1; Tdot_mid; Tdot_end;
            Thdot_1; Thdot_mid; Thdot_end];
end
