function y = modelCD_out(t, x, theta)

%      y = 1;

     nstages = theta{1};
     sensorPosition = theta{17};
     T = x(4*nstages+1 : 5*nstages);
 
     nSensors = numel(sensorPosition);
     Mu = zeros(nSensors,1);
     sigmasquare = 1;
     Sigma = sigmasquare * eye(nSensors);
     
     y = T(sensorPosition) + mvnrnd(Mu,Sigma,1)';
    
%     nstages = theta{1};
%     y = zeros(6, 1);
%     m = [0 0 0 0 0 0];
%     e = eye(6)*10^-3;
%     noise = mvnrnd(m,e,1);
%     for n = 1:6
%         y(n) = x(nstages*4 + theta{16+n})+noise(n); 
%     end

end