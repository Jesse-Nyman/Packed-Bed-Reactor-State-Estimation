function ut = Feed_Disturbances(nSteps,Feed,A,epz,prcnt,p)
    ut = ones(1,nSteps)*Feed / (A*epz);
    n = geornd(p) + 1;                                           %geornd gives number of failures we want number of tries until success
    while n < nSteps
        r = geornd(p) + 1;
        if n + r > nSteps
            ut(n:nSteps) = unifrnd((1-prcnt),(1+prcnt)) * Feed/( A*epz);
            break
        end
        ut(n:n+r) = unifrnd((1-prcnt),(1+prcnt)) * Feed/( A*epz);
        n = n+r;  
    end
    
end