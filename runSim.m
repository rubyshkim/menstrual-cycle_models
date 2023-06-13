function [T,sols,pars,period] = runSim(pars,lag,model,vars_i,IC,numCycles,lhdata)
% run simulations given pars vector and data

% compute day of LH peak for fitting purposes
if ~exist('lhdata','var')
    LH_PEAKTIME = 14;
else
    LHDATA = lhdata(:,2);
    LHTIME = lhdata(:,1);
    [~,index] = max(LHDATA);
    LH_PEAKTIME = LHTIME(index);
end

%% compute initial conditions

% transient
tstart=0; dt=0.1; tend=1000;
tspan=tstart:dt:tend;
sol = ddeSolver(@(t,y,lag) model(t,y,lag,pars), lag, IC, tspan, pars, @dde23);

% save last few cycles of solutions
sols = deval(sol,tstart:dt:tend);
tspan = tspan(end-29*3/dt:end);
sols = sols(:,end-29*3/dt:end);

% find LH peaks with min peak prominence 0.2 of max concentration
[~,locs] = findpeaks(sols(vars_i.lh,:),'MinPeakProminence',0.2*max(sols(vars_i.lh,:)));

% compute period
try
    period = tspan(locs(end))-tspan(locs(end-1));
catch
    period = NaN;
end

% compute initial condition to use
try 
    peakIndex = locs(end);
    Y_0 = sols(:,peakIndex-LH_PEAKTIME/dt); % start simulation LH_PEAKTIME days before LH peak
catch
    Y_0 = sols(:,end);
end

%% simulations

% time span
tstart=0; dt=0.1; 
if ~isnan(period)
    tend=period*numCycles;
else
    tend=29*numCycles;
end

T = tstart:dt:tend;

% solver
sol = ddeSolver(@(t,y,lag) model(t,y,lag,pars), lag, Y_0, T, pars, @dde23);

sols = deval(sol,tstart:dt:tend);

end