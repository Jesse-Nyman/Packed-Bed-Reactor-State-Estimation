function [Feed_C] = Concentration_Disturbances(nSteps, xB0t, xT0t, xP0t, prcnt1,prcnt2,prcnt3, p, Feed_C)
                                                             
    n = geornd(p) + 1;               %geornd gives number of failures we want number of tries until success
    
        while n < nSteps
            r = geornd(p) + 1;
            if n + r > nSteps
                Feed_C(n:nSteps,1) = unifrnd((1-prcnt1),(1+prcnt1))*xB0t;     %xB
                Feed_C(n:nSteps,3) = unifrnd((1-prcnt3),(1+prcnt3))*xP0t;     %xP
                Feed_C(n:nSteps,4) = unifrnd((1-prcnt2),(1+prcnt2))*xT0t;     %xT
                Feed_C(n:nSteps,2) = 1 - Feed_C(n:nSteps,1) -Feed_C(n:nSteps,4)  -Feed_C(n:nSteps,3);   %xH 
                break
            end
            Feed_C(n:n+r,1) = unifrnd((1-prcnt1),(1+prcnt1))*xB0t;     %xB
            Feed_C(n:n+r,3) = unifrnd((1-prcnt3),(1+prcnt3))*xP0t;     %xP
            Feed_C(n:n+r,4) = unifrnd((1-prcnt2),(1+prcnt2))*xT0t;     %xT
            Feed_C(n:n+r,2) = 1 - Feed_C(n:n+r,1) -Feed_C(n:n+r,4)  -Feed_C(n:n+r,3);   %xH    
            n = n+r;                                                        %new interval n:nSteps

        end
end
