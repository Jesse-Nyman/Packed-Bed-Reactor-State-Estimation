function y = modelCD_out(~, x, theta)
    % modelCD_out now returns temperature at sensor positions without noise
    nstages = theta{1};
    sensorPosition = theta{17};

    T = x(4*nstages+1 : 5*nstages);
    y = T(sensorPosition);
end
