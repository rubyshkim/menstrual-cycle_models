function sol = ddeSolver(rhs, lag, Y_0, tspan, pars, solver)
% inputs:
%   rhs         the function containing the equations
%   lag         delay to be used with solver
%   Y_0         initial conditions
%   tspan       time vector
%   pars        parameters
%   solver      which solver to use (e.g., @dde23, ...)

sol = solver(rhs, lag, Y_0, tspan, pars);

end