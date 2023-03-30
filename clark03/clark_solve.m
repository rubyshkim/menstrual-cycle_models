function [T,sols,vars_i,pars,period,lag] = clark_solve(pars,numCycles,lhdata)

model = @clark_rhs;

% specify indices of variables
vars_i.lh = 2;
vars_i.fsh = 4;
vars_i.e2 = 14;
vars_i.p4 = 15;

IC = zeros(1,13); % set initial condition. program will recompute this.

% load lag
lag = [2 1]; % delays in IH and P4

if ~exist('lhdata','var')
    [T,sols,pars,period] = runSim(pars,lag,model,vars_i,IC,numCycles);
else
    [T,sols,pars,period] = runSim(pars,lag,model,vars_i,IC,numCycles,lhdata);
end

% parameters
e_0 = pars(31); e_1 = pars(32); e_2 = pars(33); 
e_3 = pars(34); p_1 = pars(35); p_2 = pars(36); 

% variables
SeF  = sols(6,:);  PrF = sols(7,:); 
Lut3 = sols(12,:); Lut4 = sols(13,:);

% auxiliary variables
		E2=e_0+e_1*SeF+e_2*PrF+e_3*Lut4;
		P4=p_1*Lut3+p_2*Lut4;

sols = [sols; E2; P4];

end