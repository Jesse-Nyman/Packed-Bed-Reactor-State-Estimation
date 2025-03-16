function [Ut] = Concentration(nSteps, xB0t, xH0t, xP0t, xT0t, T0t, ut)

% CONCENTRATION Summary of this function goes here
% Detailed explanation goes here

   Ut   = [xB0t, xH0t, xP0t, xT0t, T0t, ut];                                % Put stuff together 
   
   prcnt = 5*10^-2;                                                         % Upper and lower limit of disturbance
                                                                            % Pass this as input argument
   
    p = 0.02;                                                               % Probabilty of success
                                                                            % Pass this as input argument
    n = geornd(p) + 1;                                                      % geornd gives number of failures we want number of tries until success
    if n > nSteps
        return
    end
    
        while n < nSteps
            r = geornd(p) + 1;
            
            if n + r > nSteps
                Ut(n:nSteps,1) = Ut(n,1);
                Ut(n:nSteps,4) = Ut(n,4);
                Ut(n:nSteps,2) = Ut(n,2);
                break
            end
            
            Ut(n:n+r,1) = unifrnd(xB0t(1)*(1-prcnt),xB0t(1)*(1+prcnt));    % Multiplications/divisions can be done outside unifrnd
            Ut(n:n+r,4) = unifrnd(xT0t(1)*(1-prcnt),xT0t(1)*(1+prcnt));    % Same
            
            Ut(n:n+r,2) = 1 - Ut(n,1) -Ut(n,4);
            n = n+r;                                                        % New interval n:nSteps

        end
end

