function [A,B,C,D] = linearizeSystem(x_bar, u_bar, theta)
% linearizeSystem.m
% This function linearizes the nonlinear model f(x,u), y=g(x) at the given 
% operating point (x_bar,u_bar). It returns the state-space matrices A,B,C,D.
%
% Inputs:
%   x_bar - Nx x 1 state vector at equilibrium/operating point
%   u_bar - Nu x 1 input vector at equilibrium/operating point
%   theta - cell array of parameters as in the original model
%
% Outputs:
%   A,B,C,D - state-space matrices linearizing:
%       dx/dt = A (x - x_bar) + B (u - u_bar)
%       y     = C (x - x_bar) + D (u - u_bar)

import casadi.*

% Unpack problem sizes
nstages = theta{1};
nstates = 6; % from original code: xB,xH,xP,xT,T,Activity each nstages long
Nx = nstates*nstages;
Nu = length(u_bar);

% Symbolics
x = MX.sym('x', Nx);
u = MX.sym('u', Nu);

% Nonlinear functions
f_func = @(x,u) modelCD(0, x, u, theta);
g_func = @(x)   modelCD_out(0, x, theta);

f_expr = f_func(x,u);
g_expr = g_func(x);

% Jacobians
Jf_x = jacobian(f_expr, x);
Jf_u = jacobian(f_expr, u);
Jg_x = jacobian(g_expr, x);
Jg_u = jacobian(g_expr, u);

% Create CasADi functions for evaluation
F_Jac_x = Function('F_Jac_x',{x,u},{Jf_x});
F_Jac_u = Function('F_Jac_u',{x,u},{Jf_u});
G_Jac_x = Function('G_Jac_x',{x},{Jg_x});
G_Jac_u = Function('G_Jac_u',{x,u},{Jg_u});

A_val = full(F_Jac_x(x_bar,u_bar));
B_val = full(F_Jac_u(x_bar,u_bar));
C_val = full(G_Jac_x(x_bar));
D_val = full(G_Jac_u(x_bar,u_bar));

A = A_val;
B = B_val;
C = C_val;
D = D_val;

end
