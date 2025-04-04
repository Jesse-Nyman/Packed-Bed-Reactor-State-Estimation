function [yout,tout,xout] = simulateSystem(varargin)
    % y_l = Simulate(F, g, x0, u1); 
    
    F = varargin{1}; 
    g = varargin{2}; 
    x0 = varargin{3};
    
    u = varargin{4};
    N = size(u,1);
        
    f = @(x,u) full(F(x,u));
    
    % aux
    Nx = size(x0,1); Ny = size(g(x0),1);
    tout = 1:N;
    
    xout = nan(Nx,N+1);   
    xout(:,1) = x0;
    
    yout = nan(Ny,N+1);   
    yout(:,1) = g(xout(:,1));
    
    % sim
    for k = 1:N
        xout(:,k+1) = f(xout(:,k), u(k,:));
        yout(:,k+1) = g(xout(:,k+1));
    end
    
    % remove the initial state/measurement
    xout = xout(:,2:end); 
    yout = yout(:,2:end);

end

 
