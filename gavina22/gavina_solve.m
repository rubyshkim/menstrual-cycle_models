function [T,sols,vars_i,pars,period,lag] = gavina_solve(pars,numCycles,lhdata)

%LdeP Updated for Gavina 2022 December 29, 2022

model = @gavina_rhs;

% specify indices of variables
vars_i.lh = 2;
vars_i.fsh = 4;
vars_i.e2 = 14;
vars_i.p4 = 15;

IC = zeros(1,13); % set initial condition. program will recompute this.
%gavina_ICs; % set ICs according to Gavina

% load lag
lag = pars(17); %1.5; % delay in IH %LdeP pars(17) addition

if ~exist('lhdata','var')
    [T,sols,pars,period] = runSim(pars,lag,model,vars_i,IC,numCycles);
else
    [T,sols,pars,period] = runSim(pars,lag,model,vars_i,IC,numCycles,lhdata);
end

% parameters
e0 = pars(32); e1 = pars(33); e2 = pars(34); 
e3 = pars(35); p0 = pars(42); p1 = pars(36); 
p2 = pars(37); p_dose = pars(44); e_dose = pars(43);

% variables
GrF  = sols(6,:);  DomF = sols(7,:); 
Lut3 = sols(12,:); Lut4 = sols(13,:);

% auxiliary variables
E2 = e0+e1*GrF+e2*DomF+e3*Lut4+e_dose;
P4 = p0+p1*Lut3+p2*Lut4+p_dose;

sols = [sols; E2; P4];

end