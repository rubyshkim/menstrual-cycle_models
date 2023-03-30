function [T,sols,vars_i,pars,period,lag] = wright_solve(pars,numCycles,lhdata)

model = @wright_rhs;

% specify indices of variables
vars_i.lh = 2;
vars_i.fsh = 4;
vars_i.e2 = 14;
vars_i.p4 = 15;

IC = zeros(1,13); % set initial condition. program will recompute this.

% load lag
lag = 1.5; % delay in IH 

if ~exist('lhdata','var')
    [T,sols,pars,period] = runSim(pars,lag,model,vars_i,IC,numCycles);
else
    [T,sols,pars,period] = runSim(pars,lag,model,vars_i,IC,numCycles,lhdata);
end

% parameters
e_0 = pars(31); e_1 = pars(32); e_2 = pars(33); 
e_3 = pars(34); p_0 = pars(35); p_1 = pars(36); 
p_2 = pars(37); p_dose = pars(44); e_dose = pars(45);

% variables
GrF  = sols(6,:);  DomF = sols(7,:); 
Lut3 = sols(12,:); Lut4 = sols(13,:);

% auxiliary variables
E2 = e_0+e_1*GrF+e_2*DomF+e_3*Lut4+e_dose;
P4 = p_0+p_1*Lut3+p_2*Lut4+p_dose;

sols = [sols; E2; P4];

end