function [T0t] = temperature_disturbance(nSteps,T0t,Temp, prcnt,p)
     n = geornd(p) + 1; 
     while n < nSteps
        r = geornd(p) + 1;
        if n + r > nSteps
            T0t(n:nSteps) = unifrnd((1-prcnt),(1+prcnt)) * Temp;
            break
        end
        T0t(n:n+r) = unifrnd((1-prcnt),(1+prcnt)) * Temp;
        n = n+r;  
     end
end

