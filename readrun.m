function [T,sols,vars_i,pars,period,lag] = readrun(paramfile,numCycles,modelsolver,datafile)

pars = readcell(paramfile); 
pars = pars(2:end,2); % Assume row 1 is labels and column 2 is parameter values
pars = cell2mat(pars);

if ~exist('datafile','var')
    [T,sols,vars_i,pars,period,lag] = modelsolver(pars,numCycles);

else
    readlh = readcell(datafile,'Sheet','LH');
    lhdata = cell2mat(readlh(2:end,:));
    
    [T,sols,vars_i,pars,period,lag] = modelsolver(pars,numCycles,lhdata);
end
end