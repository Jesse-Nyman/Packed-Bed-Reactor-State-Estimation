function F = buildIntegrator(f, d, varargin)
    import casadi.*

    stepSize = 1;
    method = 'cvodes';

    if numel(varargin) > 0
        stepSize = varargin{1};
    end
    if numel(varargin) > 1
        method = varargin{2};
    end

    Nx = d(1);
    Nu = d(2);

    x = MX.sym('x', Nx);
    u = MX.sym('u', Nu);

    DAE = struct('x', x, 'p', u, 'ode', f(x,u));
    options_struct = struct('tf', stepSize, 'abstol',1e-10,'reltol',1e-10,'max_num_steps',1e5);

    F_integrator = integrator('F_integrator', method, DAE, options_struct);
    F_res = F_integrator('x0', x, 'p', u);

    F = Function('F', {x,u}, {F_res.xf}, {'x','u'}, {'x_next'});
end
